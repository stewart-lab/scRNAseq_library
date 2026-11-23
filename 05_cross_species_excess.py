#!/usr/bin/env python3
"""
04 -- other-species reads in species-called cells: is there more than ambient?

Takes the combined human+pig STARsolo matrix and the species calls from
03_classify_species.py and asks, separately for human reads in pig cells
("human_in_pig") and pig reads in human cells ("pig_in_human"), whether the
other-species signal looks like ambient RNA or like something cell-associated
(e.g. transfer by vesicles / nanotubes / fusion / phagocytosis).

Why the tests are framed the way they are
-----------------------------------------
rho (contamination fraction) from 03 is itself computed FROM the other-species
reads, so "total foreign UMIs vs rho-expected" is equal by construction and
tests nothing.  Everything here therefore tests PATTERN, not overall level:

1. Gene composition (rho-free).  Within a group of cells, is gene g's share of
   the foreign reads higher than its share of the same species' reads in the
   soup (empty droplets)?  Exact conditional test:
       O_g | O_g + K_g  ~  Binomial(O_g + K_g,  T_O / (T_O + T_K))
   O = foreign UMIs of gene g in the cells, K = UMIs of g in the soup, T = sums
   over the foreign species' genes.  Accounts for sampling noise in the soup
   profile too.  enrichment = (O_g/T_O) / (K_g/T_K).
2. Cell level (needs no clusters, better with them).  Ambient UMIs grow
   sub-linearly with cell size, so foreign UMIs per cell are modelled as
       y_c ~ NegBin(mu_c),  log mu_c = a + b * log(N_c)
   fitted on all cells of that species.  Per-cell upper-tail p-values flag
   outlier cells; per-cluster observed/expected flags clusters carrying more
   foreign RNA than other clusters of similar cell size.  This is RELATIVE to
   the species as a whole -- a cluster can only stand out against the others.
   With clusters, genes are also tested cluster-vs-rest (`q_vs_rest`): is g a
   larger share of the foreign reads in this cluster than in the other
   clusters of the same species?  Soup-composition and mismapping artifacts
   present in every cluster cancel, so this is the column to read for
   cell-type-specific transfer (mismapped genes can still differ if the host
   ortholog is cluster-specific -- check `mismap_like` in the pooled rows).
3. Ambient load per droplet.  foreign / (foreign species' soup share) at a
   common reference cell size, for each species.  Each direction uses its own
   tracer, so a large asymmetry is semi-independent evidence; it can also be
   technical (cell fragility, surface-bound debris differ between species).

Mismapping ceiling.  A host-species read that maps uniquely to a foreign gene
(assembly gap, missing/divergent ortholog) appears both in the cells and in
the soup, which is mostly host RNA.  That inflates its foreign-read share by a
fixed factor:
    ceiling = (1 / foreign fraction in cells) / (host soup share / foreign soup share)
(~6.6x for human genes in pig cells here).  Genes enriched near the ceiling
are consistent with being 100% mismapped host reads, and are flagged
`mismap_like` (enrichment >= --mismap-frac x ceiling, and not significantly
above the ceiling).  Enrichment well ABOVE
the ceiling (`q_above_ceiling`) cannot be explained by ambient or mismapping
alone.  Genes below
it are ambiguous: check reads in the BAM (nM tag, does the read also align to
the host genome?) before calling transfer.

Intronic check (--exon-solo-dir = exon-only Gene/raw of the same STARsolo
run; --solo-dir should then be GeneFull / GeneFull_Ex50pAS).  Ambient and
mismapped reads are mostly spliced; a fragment of the other species' nucleus
brings unspliced pre-mRNA.  Adds O_exonic / exonic_frac / soup_exonic_frac per
gene, foreign_exonic_frac + intronic_q + nuclear_like per cell (excess cells
whose foreign reads are more intronic than ordinary cells', one-sided
binomial, BH).  Per cluster: q_nuclear_enriched (more nuclear-like cells than
the rest of the species, Fisher) and q_intronic (are the cluster's ORDINARY,
non-nuclear-like cells' foreign reads more intronic than the rest's --
diffuse nuclear-type RNA; quasi-binomial, between-cell dispersion estimated
from ordinary cells).  03 run with the same option
calls these cells `nuclear_doublet`; 04 keeps them (by call_purity) so they
can be inspected, with `called_nuclear_doublet_03` marking them.

Caveats: soup (empty droplets) is enriched for cell-free RNA relative to the
ambient stuck to cells, so moderate (~2x) gene-level differences can be
technical, especially MT / short transcripts.  Check top hits for pig<->human
mismapping (ortholog present? reads with mismatches?) before calling transfer.
Use the UNIQUE matrix (matrix.mtx), not UniqueAndMult-EM.mtx: EM splits
ambiguous reads using each cell's own-species prior.

Clusters
--------
Optional --clusters CSV/TSV(s) with a barcode column and a cluster column, e.g.
from Seurat:  write.csv(data.frame(barcode=colnames(obj),
                                   cluster=obj$seurat_clusters), "pig_clusters.csv")
Cluster each species on ITS OWN genes only (e.g. the single-species
alignments); if other-species genes drive clustering, "cluster X is enriched
for foreign RNA" becomes circular.  Barcodes are matched on the longest
A/C/G/T/N run, so Seurat prefixes/suffixes (lane1_, -1) are tolerated; use
--barcode-filter to keep only one lane's rows from a merged object.

usage
-----
  python3 05_cross_species_excess.py \
      --solo-dir   OUT/SAMPLE/Solo.out/GeneFull/raw \
      --assignment OUT/SAMPLE/classified_cells/species_assignment.csv.gz \
      [--clusters pig_clusters.csv human_clusters.csv --cluster-col cluster] \
      --outdir     OUT/SAMPLE/cross_species
"""
import argparse, gzip, json, os, re, sys
import numpy as np
import pandas as pd
from scipy import stats

DIRECTIONS = {  # name: (cell species, foreign species)
    "human_in_pig": ("pig", "human"),
    "pig_in_human": ("human", "pig"),
}


# ---------------------------------------------------------------- io helpers
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _find(d, stem):
    for c in (f"{stem}.gz", stem):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{stem}[.gz] not found in {d}")


def bh(p):
    """Benjamini-Hochberg; NaN-aware."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        out = np.empty_like(pv)
        out[o] = np.minimum(r, 1)
        q[ok] = out
    return q


def norm_bc(b):
    runs = re.findall(r"[ACGTN]{8,}", str(b))
    return max(runs, key=len) if runs else str(b)


# ------------------------------------------------------------- clusters
def load_clusters(paths, bc_col, cl_col, bc_filter):
    frames = []
    for p in paths:
        t = pd.read_csv(p, sep=None, engine="python")
        b = bc_col if bc_col else t.columns[0]
        c = cl_col if cl_col in t.columns else t.columns[1]
        if b not in t.columns:
            sys.exit(f"{p}: barcode column '{b}' not found ({list(t.columns)})")
        t = t[[b, c]].rename(columns={b: "raw", c: "cluster"})
        if bc_filter:
            t = t[t["raw"].astype(str).str.contains(bc_filter, regex=True)]
        frames.append(t)
    t = pd.concat(frames, ignore_index=True)
    t["cb"] = t["raw"].map(norm_bc)
    t["cluster"] = t["cluster"].astype(str)
    dup = t["cb"].duplicated(keep=False)
    if dup.any():
        sys.stderr.write(f"warning: {dup.sum()} cluster rows share a barcode after "
                         "normalisation (merged lanes? use --barcode-filter); dropped\n")
        t = t[~dup]
    return t.set_index("cb")["cluster"]


# ------------------------------------------------------------- matrix
def grouped_gene_sums(solo_dir, assign_cb, grp_of_row, ngroups, row_masks=None,
                      chunksize=5_000_000):
    """Gene x group UMI sums streamed from matrix.mtx.  grp_of_row gives the
    group index (or -1) for each row of the assignment table.  With row_masks
    ({name: gene mask}), also returns per-assignment-row sums over each mask."""
    feats = pd.read_csv(_find(solo_dir, "features.tsv"), sep="\t", header=None,
                        usecols=[0, 1], names=["gene_id", "gene_name"])
    bcs = pd.read_csv(_find(solo_dir, "barcodes.tsv"), sep="\t", header=None)[0].values
    pos = pd.Index(assign_cb).get_indexer(bcs)  # matrix col -> assignment row
    grp_of_col = np.where(pos >= 0, grp_of_row[np.maximum(pos, 0)], -1)

    mtx = _find(solo_dir, "matrix.mtx")
    n_head = 0
    with _open(mtx) as fh:
        for line in fh:
            n_head += 1
            if not line.startswith("%"):
                ngene, nbc, nnz = (int(x) for x in line.split())
                break
    if ngene != len(feats) or nbc != len(bcs):
        sys.exit("matrix header does not match features/barcodes length")
    sys.stderr.write(f"matrix: {ngene} x {nbc}, {nnz} nonzero\n")

    acc = np.zeros(ngene * ngroups, np.int64)
    row_masks = row_masks or {}
    nrow = len(assign_cb)
    rows = {k: np.zeros(nrow, np.int64) for k in row_masks}
    reader = pd.read_csv(mtx, sep=r"\s+", skiprows=n_head, header=None,
                         names=["g", "c", "v"], dtype=np.int64,
                         chunksize=chunksize, engine="c")
    for ch in reader:
        g = ch["g"].values - 1
        c = ch["c"].values - 1
        v = ch["v"].values
        k = grp_of_col[c]
        sel = k >= 0
        if sel.any():
            acc += np.bincount(g[sel] * ngroups + k[sel], weights=v[sel],
                               minlength=ngene * ngroups).astype(np.int64)
        if row_masks:
            r = pos[c]
            for name, m in row_masks.items():
                s = (r >= 0) & m[g]
                if s.any():
                    rows[name] += np.bincount(r[s], weights=v[s],
                                              minlength=nrow).astype(np.int64)
    return feats, acc.reshape(ngene, ngroups), rows


# ------------------------------------------------------------- tests
def composition_test(O, K, min_obs, fold=1.0):
    """Is gene g's share of O more than `fold` x its share of K?
    Conditional binomial on O_g + K_g.  Returns E (at fold 1), enrichment,
    p/q up, p/q down (down only meaningful at fold 1)."""
    TO, TK = O.sum(), K.sum()
    n = O + K
    E = TO * K / TK
    enr = ((O + 0.5) / TO) / ((K + 0.5) / TK)
    p_f = fold * TO / (fold * TO + TK)
    up = O >= min_obs
    dn = E >= min_obs
    p_up = np.where(up, stats.binom.sf(O - 1, n, p_f), np.nan)
    p_dn = np.where(dn, stats.binom.cdf(O, n, TO / (TO + TK)), np.nan)
    return E, enr, p_up, bh(p_up), p_dn, bh(p_dn)


def fit_size_model(y, N, iters=50):
    """Poisson IRLS for log mu = a + b log N, then moment NB dispersion."""
    y = np.asarray(y, float)
    X = np.c_[np.ones(len(N)), np.log(N)]
    beta = np.linalg.lstsq(X, np.log(y + 0.5), rcond=None)[0]
    for _ in range(iters):
        eta = X @ beta
        mu = np.exp(eta)
        z = eta + (y - mu) / mu
        new = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (mu * z))
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    mu = np.exp(X @ beta)
    phi = max(0.0, float(((y - mu) ** 2 - mu).sum() / (mu ** 2).sum()))
    return beta, mu, phi


def nb_sf(y, mu, phi):
    if phi <= 0:
        return stats.poisson.sf(y - 1, mu)
    r = 1.0 / phi
    return stats.nbinom.sf(y - 1, r, r / (r + mu))


# ------------------------------------------------------------- plots
def plot_composition(gt, path, top=8):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 300,
                         "savefig.bbox": "tight"})
    dirs = [d for d in DIRECTIONS if ((gt.direction == d) & (gt.group == "all")).any()]
    fig, axes = plt.subplots(1, len(dirs), figsize=(3.4 * len(dirs), 3.2), squeeze=False)
    for ax, d in zip(axes[0], dirs):
        s = gt[(gt.direction == d) & (gt.group == "all") & ((gt.O > 0) | (gt.E > 0.5))]
        mis = s.mismap_like.astype(bool)
        abv = s.q_above_ceiling < 0.05
        # above the ceiling line as a point estimate, but not significantly
        # (usually few soup UMIs -> uncertain enrichment)
        ns_abv = mis & (s.enrichment > s.mismap_ceiling)
        mis_in = mis & ~ns_abv
        hit = (s.q_up < 0.05) & ~mis & ~abv
        ax.scatter(s.E + 0.5, s.O + 0.5, s=2, c="0.6", linewidths=0, rasterized=True)
        ax.scatter(s.E[mis_in] + 0.5, s.O[mis_in] + 0.5, s=4, c="#e08214", linewidths=0,
                   label=f"enriched, mismap-like (n={int(mis_in.sum())})")
        ax.scatter(s.E[ns_abv] + 0.5, s.O[ns_abv] + 0.5, s=8, facecolors="none",
                   edgecolors="#e08214", linewidths=0.6,
                   label=f"above ceiling line, not significant (n={int(ns_abv.sum())})")
        ax.scatter(s.E[hit] + 0.5, s.O[hit] + 0.5, s=4, c="#b2182b", linewidths=0,
                   label=f"enriched, below ceiling (n={int(hit.sum())})")
        ax.scatter(s.E[abv] + 0.5, s.O[abv] + 0.5, s=8, c="#542788", linewidths=0,
                   label=f"significantly above ceiling, q<0.05 (n={int(abv.sum())})")
        lim = [0.5, max(s.E.max(), s.O.max()) * 2]
        ax.plot(lim, lim, ls="--", lw=0.6, c="0.3")
        c0 = s.mismap_ceiling.iloc[0]
        ax.plot(lim, [v * c0 for v in lim], ls=":", lw=0.6, c="#e08214")
        hit = s.q_up < 0.05
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lim); ax.set_ylim(lim)
        for _, r in s[hit].nlargest(top, "O").iterrows():
            ax.annotate(r.gene_name.split("_", 1)[-1], (r.E + 0.5, r.O + 0.5),
                        fontsize=5, xytext=(2, 2), textcoords="offset points")
        ax.text(0.98, 0.02,
                f"dashed: soup expectation (1x)\ndotted: mismap ceiling ({c0:.1f}x)\n"
                "purple only if significantly above the ceiling;\n"
                "points above the dotted line with few soup UMIs\n"
                "are too uncertain to call",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=5, color="0.35")
        cs, fs = DIRECTIONS[d]
        ax.set_xlabel(f"expected {fs} UMIs from soup composition + 0.5")
        ax.set_ylabel(f"observed {fs} UMIs in {cs} cells + 0.5")
        ax.set_title(d.replace("_", " "))
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2),
                  fontsize=6, ncol=1)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def plot_cells(cells, fits, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dirs = [d for d in DIRECTIONS if d in fits]
    fig, axes = plt.subplots(1, len(dirs), figsize=(3.4 * len(dirs), 3.2), squeeze=False)
    for ax, d in zip(axes[0], dirs):
        s = cells[cells.direction == d]
        cl = s.cluster.astype(str)
        levels = cl.value_counts().index[:20]
        cmap = plt.get_cmap("tab20")
        for i, lv in enumerate(levels):
            m = cl == lv
            ax.scatter(s.N[m], s.foreign[m] + 0.5, s=1.5, color=cmap(i % 20),
                       linewidths=0, alpha=0.6, rasterized=True,
                       label=lv if len(levels) > 1 else None)
        a, b = fits[d]["a"], fits[d]["b"]
        x = np.linspace(np.log(s.N.min()), np.log(s.N.max()), 50)
        ax.plot(np.exp(x), np.exp(a + b * x), c="k", lw=0.8,
                label=f"fit: N^{b:.2f}")
        sig = s.q < 0.05
        nuc = s["nuclear_like"].astype(bool) if "nuclear_like" in s else pd.Series(False, index=s.index)
        ax.scatter(s.N[sig & ~nuc], s.foreign[sig & ~nuc] + 0.5, s=6, facecolors="none",
                   edgecolors="#b2182b", linewidths=0.5,
                   label=f"excess, q<0.05 (n={int((sig & ~nuc).sum())})")
        if nuc.any():
            ax.scatter(s.N[nuc], s.foreign[nuc] + 0.5, s=6, facecolors="none",
                       edgecolors="#e08214", linewidths=0.7,
                       label=f"excess + intronic (n={int(nuc.sum())})")
        ax.set_xscale("log"); ax.set_yscale("log")
        cs, fs = DIRECTIONS[d]
        ax.set_xlabel(f"{cs} cell UMIs (species-assignable)")
        ax.set_ylabel(f"{fs} UMIs per cell + 0.5")
        ax.set_title(d.replace("_", " "))
        ax.legend(frameon=False, fontsize=5, markerscale=3, ncol=2)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--solo-dir", required=True,
                    help="combined-reference STARsolo raw matrix dir (unique matrix.mtx)")
    ap.add_argument("--assignment", required=True,
                    help="species_assignment.csv.gz from 03_classify_species.py")
    ap.add_argument("--clusters", nargs="*", default=[],
                    help="CSV/TSV(s) with barcode + cluster columns (optional)")
    ap.add_argument("--barcode-col", default=None, help="default: first column")
    ap.add_argument("--cluster-col", default="cluster", help="default: 'cluster' or 2nd column")
    ap.add_argument("--barcode-filter", default=None,
                    help="regex; keep only cluster rows whose barcode matches (e.g. '^lane1_')")
    ap.add_argument("--human-prefix", default="GRCh38")
    ap.add_argument("--pig-prefix", default="Sscrofa11")
    ap.add_argument("--soup-max-umi", type=int, default=100)
    ap.add_argument("--min-obs", type=int, default=3,
                    help="min observed (enrichment) / expected (depletion) UMIs to test a gene")
    ap.add_argument("--min-cells", type=int, default=20, help="min cells to test a cluster")
    ap.add_argument("--mismap-frac", type=float, default=0.5,
                    help="flag enriched genes with enrichment >= this x mismap ceiling")
    ap.add_argument("--exon-solo-dir", default=None,
                    help="exon-only Gene/raw matrix of the same STARsolo run: adds "
                         "exonic fractions and the nuclear (intronic) check")
    ap.add_argument("--intronic-fdr", type=float, default=0.05)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    pfx = {"human": a.human_prefix, "pig": a.pig_prefix}

    sa = pd.read_csv(a.assignment)
    need = {"cb", "human", "pig", "total", "call_purity"}
    if need - set(sa.columns):
        sys.exit(f"assignment table lacks {need - set(sa.columns)} -- re-run 03")
    if (sa["total"] > 1e6).any():
        sys.exit("assignment has a barcode with >1e6 UMIs: produced by the old 03 "
                 "matrix-header bug. Re-run 03_classify_species.py first.")

    # ---- cluster labels
    sa["cluster"] = "all"
    if a.clusters:
        cl = load_clusters(a.clusters, a.barcode_col, a.cluster_col, a.barcode_filter)
        lab = sa["cb"].map(cl)
        called = sa["call_purity"].isin(["human", "pig"])
        n_match = int(lab[called].notna().sum())
        sys.stderr.write(f"clusters: {len(cl)} labelled barcodes, {n_match} match "
                         f"species-called cells ({n_match / called.sum():.1%} of called)\n")
        if n_match == 0:
            sys.exit("no cluster barcodes match -- check --barcode-col / --barcode-filter")
        sa["cluster"] = lab.fillna("unclustered")

    # ---- groups: 0 = soup, then (species, cluster) cells
    soup = (sa["total"] > 0) & (sa["total"] <= a.soup_max_umi)
    keys = [("soup", "soup")]
    grp = np.full(len(sa), -1, np.int64)
    grp[soup.values] = 0
    for sp in ("pig", "human"):
        m = (sa["call_purity"] == sp).values
        for c in sorted(sa.loc[m, "cluster"].unique()):
            keys.append((sp, c))
            grp[m & (sa["cluster"] == c).values] = len(keys) - 1
    feats, G, _ = grouped_gene_sums(a.solo_dir, sa["cb"].values, grp, len(keys))
    gid = feats["gene_id"].values.astype(str)
    is_sp = {sp: np.char.startswith(gid, pfx[sp] + "_") for sp in pfx}

    # exon-only counts: per gene x group, and per cell per species
    GX = None
    if a.exon_solo_dir:
        sys.stderr.write("reading exon-only matrix for the intronic check\n")
        fx, GX, rx = grouped_gene_sums(
            a.exon_solo_dir, sa["cb"].values, grp, len(keys),
            row_masks={sp: is_sp[sp] for sp in pfx})
        if not np.array_equal(fx["gene_id"].values.astype(str), gid):
            sys.exit("--exon-solo-dir features differ from --solo-dir; "
                     "use the Gene/raw matrix from the same STARsolo run")
        for sp in pfx:
            sa[sp + "_exonic_04"] = np.minimum(rx[sp], sa[sp].values)  # exonic <= total

    K_all = G[:, 0]
    f_soup = {sp: K_all[is_sp[sp]].sum() / K_all[is_sp["human"] | is_sp["pig"]].sum()
              for sp in pfx}
    sys.stderr.write(f"soup: {int(soup.sum())} barcodes, {int(K_all.sum())} UMIs, "
                     f"{f_soup['human']:.4f} human\n")

    # ---- 1. gene composition per direction x group (incl. pooled "all")
    gene_rows, summ = [], {"soup_barcodes": int(soup.sum()),
                           "soup_frac_human": round(float(f_soup["human"]), 4),
                           "directions": {}}
    for d, (cs, fs) in DIRECTIONS.items():
        cols = [i for i, (sp, _) in enumerate(keys) if sp == cs]
        if not cols:
            continue
        groups = {"all": G[:, cols].sum(1)}
        groups_x = {"all": GX[:, cols].sum(1)} if GX is not None else {}
        if a.clusters:
            for i in cols:
                if (sa["call_purity"].eq(cs) & sa["cluster"].eq(keys[i][1])).sum() >= a.min_cells:
                    groups[keys[i][1]] = G[:, i]
                    if GX is not None:
                        groups_x[keys[i][1]] = GX[:, i]
        fm = is_sp[fs]
        both = is_sp["human"] | is_sp["pig"]
        K = K_all[fm]
        ceilings = {}
        pooled = groups["all"][fm]
        for gname, col in groups.items():
            O = col[fm]
            E, enr, pu, qu, pd_, qd = composition_test(O, K, a.min_obs)
            x_cells = O.sum() / col[both].sum()
            ceil = (1 / x_cells) / (f_soup[cs] / f_soup[fs])
            ceilings[gname] = ceil
            # enrichment significantly above what pure mismapping could give
            _, _, _, q_ab, _, _ = composition_test(O, K, a.min_obs, fold=ceil)
            row = {
                "direction": d, "group": gname,
                "gene_id": gid[fm], "gene_name": feats["gene_name"].values[fm],
                "O": O, "K_soup": K, "E": np.round(E, 3),
                "enrichment": np.round(enr, 3), "p_up": pu, "q_up": qu,
                "p_down": pd_, "q_down": qd,
                "mismap_ceiling": round(float(ceil), 2),
                # near the ceiling but NOT significantly above it
                "mismap_like": (qu < 0.05) & (enr >= a.mismap_frac * ceil) & ~(q_ab < 0.05),
                "q_above_ceiling": q_ab}
            # cluster vs the other clusters of the same species: artifacts
            # shared by all clusters (soup composition, mismapping) cancel
            if gname != "all":
                rest = pooled - O
                _, enr_r, p_r, q_r, _, _ = composition_test(O, rest, a.min_obs)
                row.update({"O_rest": rest, "enrichment_vs_rest": np.round(enr_r, 3),
                            "p_vs_rest": p_r, "q_vs_rest": q_r})
            if GX is not None:
                Ox = np.minimum(groups_x[gname][fm], O)
                Kx = np.minimum(GX[:, 0][fm], K)
                with np.errstate(invalid="ignore", divide="ignore"):
                    row.update({"O_exonic": Ox,
                                "exonic_frac": np.round(np.where(O > 0, Ox / O, np.nan), 3),
                                "soup_exonic_frac": np.round(np.where(K > 0, Kx / K, np.nan), 3)})
            gene_rows.append(pd.DataFrame(row))
        # composition agreement (pooled): Spearman of O vs E for genes with E>=5
        g0 = gene_rows[-len(groups)]
        ok = g0.E >= 5
        rho_s = stats.spearmanr(g0.O[ok], g0.E[ok])[0] if ok.sum() > 10 else np.nan
        summ["directions"][d] = {
            "foreign_umis_in_cells": int(g0.O.sum()),
            "soup_umis_of_foreign_species": int(K.sum()),
            "spearman_obs_vs_soup_expected_E>=5": round(float(rho_s), 3),
            "n_genes_enriched_q05_pooled": int((g0.q_up < 0.05).sum()),
            "mismap_ceiling_pooled": round(float(ceilings["all"]), 2),
            "n_enriched_mismap_like": int(g0.mismap_like.sum()),
            "n_enriched_above_ceiling_q05": int((g0.q_above_ceiling < 0.05).sum()),
            "n_genes_depleted_q05_pooled": int((g0.q_down < 0.05).sum()),
        }
    gt = pd.concat(gene_rows, ignore_index=True)
    gt.sort_values(["direction", "group", "p_up"]).to_csv(
        os.path.join(a.outdir, "gene_excess.csv.gz"), index=False)

    # ---- 2. cell level: size model, per-cell + per-cluster excess
    cell_rows, clus_rows, fits = [], [], {}
    for d, (cs, fs) in DIRECTIONS.items():
        c = sa[sa["call_purity"] == cs]
        if len(c) < 20:
            continue
        y = c[fs].values
        N = c["total"].values.astype(float)
        beta, mu, phi = fit_size_model(y, N)
        p = nb_sf(y, mu, phi)
        q = bh(p)
        exc = q < 0.05
        fits[d] = {"a": float(beta[0]), "b": float(beta[1]), "phi": phi}
        crow = {"cb": c["cb"].values, "direction": d, "cluster": c["cluster"].values,
                "N": N.astype(int), "foreign": y, "expected": np.round(mu, 3),
                "ratio": np.round((y + 0.5) / (mu + 0.5), 3), "p": p, "q": q}
        if GX is not None:
            # nuclear check: are the excess cells' foreign reads more intronic
            # than the foreign reads of ordinary (ambient-only) cells?
            ex = c[fs + "_exonic_04"].values
            ref = ex[~exc].sum() / max(y[~exc].sum(), 1)
            pi = np.full(len(y), np.nan)
            sel = exc & (y > 0)
            pi[sel] = stats.binom.cdf(ex[sel], y[sel], ref)
            qi = bh(pi)
            with np.errstate(invalid="ignore", divide="ignore"):
                crow.update({"foreign_exonic_frac": np.round(np.where(y > 0, ex / y, np.nan), 4),
                             "intronic_q": qi,
                             "nuclear_like": sel & (qi < a.intronic_fdr)})
            # between-cell (extra-binomial) variability of the exonic fraction
            # among ordinary cells, for the quasi-binomial cluster test
            o_ = ~exc & (y > 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                pr = (ex[o_] - y[o_] * ref) ** 2 / (y[o_] * ref * (1 - ref))
            phi_x = max(1.0, float(pr.sum() / max(o_.sum() - 1, 1)))
            # cluster intronic test uses non-nuclear-like cells only; nuclear-
            # like cells are tested separately as a count (Fisher)
            nl = crow["nuclear_like"]
            keep_x = ~nl
            ex_tot, y_tot = ex[keep_x].sum(), y[keep_x].sum()
            n_nl, n_all = int(nl.sum()), len(y)
        if "nuclear_doublet" in c.columns:
            crow["called_nuclear_doublet_03"] = c["nuclear_doublet"].values
        cell_rows.append(pd.DataFrame(crow))
        for cname, idx in pd.Series(np.arange(len(c))).groupby(c["cluster"].values):
            ii = idx.values
            o, e = y[ii].sum(), mu[ii].sum()
            var = (mu[ii] + phi * mu[ii] ** 2).sum()
            z = (o - e) / np.sqrt(var) if var > 0 else np.nan
            if GX is not None:
                jj = ii[keep_x[ii]]          # this cluster's non-nuclear-like cells
                ox, oy = ex[jj].sum(), y[jj].sum()
                # cluster vs the rest of the species, quasi-binomial, one-sided:
                # are this cluster's (ordinary) foreign reads more intronic?
                r_rest = (ex_tot - ox) / max(y_tot - oy, 1)
                pz = np.nan
                if oy and len(ii) >= a.min_cells and 0 < r_rest < 1 and y_tot > oy:
                    zx = (ox - oy * r_rest) / np.sqrt(phi_x * oy * r_rest * (1 - r_rest))
                    pz = float(stats.norm.cdf(zx))
                k = int(nl[ii].sum())
                pf = (float(stats.fisher_exact([[k, len(ii) - k],
                                                [n_nl - k, n_all - len(ii) - (n_nl - k)]],
                                               alternative="greater")[1])
                      if len(ii) >= a.min_cells and n_all > len(ii) else np.nan)
                xrow = {"foreign_exonic_frac": round(float(ex[ii].sum() / o), 4) if o else np.nan,
                        "foreign_exonic_frac_ordinary": round(float(ox / oy), 4) if oy else np.nan,
                        "foreign_exonic_frac_rest": round(float(r_rest), 4),
                        "p_intronic": pz,
                        "n_nuclear_like": k,
                        "frac_cells_nuclear_like": round(k / len(ii), 4),
                        "p_nuclear_enriched": pf}
            else:
                xrow = {}
            clus_rows.append({**xrow,
                "direction": d, "cluster": cname, "n_cells": len(ii),
                "median_N": float(np.median(N[ii])),
                "foreign_obs": int(o), "foreign_exp": round(float(e), 1),
                "obs_over_exp": round(float(o / e), 3) if e else np.nan,
                "z": round(float(z), 2),
                "p_up": float(stats.norm.sf(z)) if len(ii) >= a.min_cells else np.nan,
                "rho_tracer": round(float(o / N[ii].sum() / f_soup[fs]), 4),
                "frac_cells_q05": round(float(exc[ii].mean()), 4)})
        # ambient load per droplet at a common reference size (see docstring #3)
        summ["directions"][d].update({
            "size_model_log_mu": f"{beta[0]:.3f} + {beta[1]:.3f}*log(N)",
            "nb_dispersion": round(phi, 4),
            "n_cells": int(len(c)),
            "n_cells_excess_q05": int(exc.sum()),
            "rho_tracer_pooled": round(float(y.sum() / N.sum() / f_soup[fs]), 4),
        })
        if GX is not None:
            summ["directions"][d].update({
                "foreign_exonic_frac_ordinary_cells": round(float(ref), 4),
                "foreign_exonic_frac_excess_cells":
                    round(float(ex[exc].sum() / y[exc].sum()), 4) if exc.any() else None,
                "n_cells_nuclear_like": int(crow["nuclear_like"].sum()),
            })
    cells = pd.concat(cell_rows, ignore_index=True)
    cells.to_csv(os.path.join(a.outdir, "cell_excess.csv.gz"), index=False)
    ct = pd.DataFrame(clus_rows)
    ct["q_up"] = np.nan
    xcols = [c for c in ("p_intronic", "p_nuclear_enriched") if c in ct]
    for c in xcols:
        ct["q" + c[1:]] = np.nan
    for d in ct.direction.unique():
        m = ct.direction == d
        ct.loc[m, "q_up"] = bh(ct.loc[m, "p_up"].values)
        for c in xcols:
            ct.loc[m, "q" + c[1:]] = bh(ct.loc[m, c].values)
    ct.to_csv(os.path.join(a.outdir, "cluster_excess.csv"), index=False)

    called = sa[sa["call_purity"].isin(["human", "pig"])]
    n_ref = float(np.median(called["total"]))
    summ["reference_cell_umis"] = n_ref
    for d, (cs, fs) in DIRECTIONS.items():
        if d in fits:
            amb = np.exp(fits[d]["a"] + fits[d]["b"] * np.log(n_ref)) / f_soup[fs]
            summ["directions"][d]["ambient_umis_per_droplet_at_ref"] = round(float(amb), 1)
    if all(d in fits for d in DIRECTIONS):
        x = summ["directions"]
        summ["ambient_ratio_pig_over_human_droplets"] = round(
            x["human_in_pig"]["ambient_umis_per_droplet_at_ref"]
            / x["pig_in_human"]["ambient_umis_per_droplet_at_ref"], 3)

    with open(os.path.join(a.outdir, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    plot_composition(gt, os.path.join(a.outdir, "gene_composition.png"))
    plot_cells(cells, fits, os.path.join(a.outdir, "cell_foreign_vs_size.png"))

    print(json.dumps(summ, indent=2))
    if a.clusters:
        print(ct.sort_values(["direction", "p_up"]).to_string(index=False))


if __name__ == "__main__":
    main()
