#!/usr/bin/env python3
"""
03 -- assign each barcode to pig or human from the combined-reference matrix.

Reads a STARsolo raw feature-barcode matrix, sums UMIs per barcode by species
prefix, and calls each barcode human / pig / mixed / unassigned.

Method notes
------------
* Assignment uses a purity threshold on species-assignable UMIs, not a simple
  majority.  A >50% rule forces every heterotypic doublet and every ambient-
  dominated empty droplet into one species or the other.
* Run this on the GeneFull matrix.  Gene-level counting depends on annotated
  3' UTRs, which are systematically truncated in Sscrofa11.1, so Gene-level
  totals understate pig and shift the ratio.  GeneFull (intronic reads
  included) is much less sensitive to that.
* Mitochondrial UMIs are reported separately.  Pig and human mtDNA are
  divergent enough to map uniquely almost everywhere and MT transcripts are a
  large fraction of a typical library, so the MT-only ratio is an independent,
  annotation-insensitive check on the whole-transcriptome call.
* Per-species purity cuts (--human-purity / --pig-purity).  Ambient RNA is
  dominated by the majority species, so the two arms are not symmetric: with
  ~10% contamination and a mostly-pig soup, clean human cells sit at ~88-90%
  human while clean pig cells sit at >95% pig.  A single symmetric 0.90 cut
  then throws most human cells into "mixed".  Look at the histogram in
  barnyard.png and set each cut below its arm's peak.
* The minority fraction inside confidently-assigned cells is a direct
  measurement of ambient (soup) contamination -- no modelling assumptions.
  It is converted to a SoupX contamination fraction (rho) by dividing by that
  species' share of the soup, measured from empty droplets
  (0 < UMIs <= --soup-max-umi, SoupX's default soupRange):
      rho_pig   = (human UMIs / all UMIs, pooled over pig cells)   / f_human_soup
      rho_human = (pig UMIs   / all UMIs, pooled over human cells) / (1 - f_human_soup)
  rho counts ALL ambient (same-species included); the other species just
  makes it visible.  Pass rho_pig / rho_human to SoupX setContaminationFraction()
  on the single-species alignments.  Pooled (UMI-weighted) fractions are used,
  matching SoupX's definition.  Ambient grows sub-linearly with cell size
  (~N^0.6 here), so small cells have a higher contamination fraction than
  large ones; rho is an average.  A noisy per-cell estimate is written
  to the `rho_cell` column (average it per cluster if wanted).
* Heterotypic doublets are observed directly; homotypic ones are invisible.
  total_rate ~ observed_hetero / (2 * f_human * f_pig).
* Reporter override: GRCh38_NLuc_P2A_tdTomato and GRCh38_NeoR are transgenes
  carried only by the human cells.  Any barcode passing --min-umi with
  >= --min-reporter-umi reporter UMIs is called human regardless of purity.
  The purity-only call is kept in `call_purity` and overrides are flagged in
  `reporter_override`.  Caveat: ambient reporter mRNA can land in pig droplets
  and human+pig doublets carry reporter reads too, so check how many pig/mixed
  calls were flipped (summary.json) and raise --min-reporter-umi if the
  overrides look like soup rather than real human cells.
* Nuclear doublets (needs --exon-solo-dir = the exon-only Gene/raw matrix of
  the same STARsolo run).  Some cells pass the purity cut but carry a small
  piece of the other species' NUCLEUS (nuclear fragment / chromatin-bound RNA
  from lysed cells, or fusion).  Signature: (1) more foreign UMIs than ambient
  predicts for the cell's size, and (2) those foreign reads are intronic
  (unspliced pre-mRNA), whereas ambient + mismapped reads are mostly exonic.
  (1) foreign UMIs ~ NegBin(mu), log mu = a + b log N, fitted per species;
      excess if BH q < --excess-fdr.
  (2) foreign exonic fraction (Gene / intron-inclusive matrix) tested one-sided
      against the pooled fraction of the non-excess cells (binomial);
      intronic if BH q (over excess cells) < --intronic-fdr.
  Cells with both get call = "nuclear_doublet" (call_purity is kept), are
  written to barcodes_nuclear_doublet.txt, and are left out of rho.  Takes
  precedence over the reporter override.  Excess cells whose foreign reads are
  NOT intronic stay in their species and are only flagged (foreign_excess_q).
  Run on the intron-inclusive matrix (GeneFull / GeneFull_Ex50pAS) for this.
usage
-----
  python3 03_classify_species.py \
      --solo-dir  OUT/combined/SAMPLE.Solo.out/GeneFull/raw \
      --gtf       REF/combined/combined.gtf \
      --human-prefix GRCh38 --pig-prefix Sscrofa11 \
      --min-umi 500 --human-purity 0.80 --pig-purity 0.90 \
      --reporter-genes GRCh38_NLuc_P2A_tdTomato GRCh38_NeoR --min-reporter-umi 1 \
      --exon-solo-dir OUT/combined/SAMPLE.Solo.out/Gene/raw \
      --outdir    OUT/species_call
"""
import argparse, gzip, json, os, sys
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------- io helpers
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _find(d, stem):
    for c in (f"{stem}.gz", stem):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{stem}[.gz] not found in {d}")


def mito_gene_ids(gtf, prefixes):
    """gene_ids on a mitochondrial contig, from the combined GTF."""
    wanted = set()
    for p in prefixes:
        wanted |= {f"{p}_chrM", f"{p}_MT", f"{p}_chrMT", f"{p}_M"}
    ids = set()
    with _open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[0] not in wanted:
                continue
            i = f[8].find('gene_id "')
            if i >= 0:
                ids.add(f[8][i + 9:f[8].index('"', i + 9)])
    return ids


# ------------------------------------------------------------- matrix sums
def _stream_barcode_sums(mtx, masks, nbc, chunksize=5_000_000):
    """Per-barcode UMI sums over each gene mask, streamed from a .mtx file."""
    n_head = 0  # comment lines + the dimension line, all skipped below
    with _open(mtx) as fh:
        for line in fh:
            n_head += 1
            if not line.startswith("%"):
                break
    acc = {k: np.zeros(nbc, np.int64) for k in masks}
    # skip the header *including* the dimension line; reading that line as data
    # adds nnz UMIs to (last gene, last barcode) -- a fake 60M-UMI "cell"
    reader = pd.read_csv(mtx, sep=r"\s+", skiprows=n_head, header=None,
                         names=["g", "c", "v"], dtype=np.int64,
                         chunksize=chunksize, engine="c")
    for chunk in reader:
        g = chunk["g"].values - 1
        c = chunk["c"].values - 1
        v = chunk["v"].values
        for key, m in masks.items():
            sel = m[g]
            if sel.any():
                acc[key] += np.bincount(c[sel], weights=v[sel],
                                        minlength=nbc).astype(np.int64)
    return acc


def species_sums(solo_dir, hpfx, ppfx, gtf=None, reporter_genes=(),
                 exon_dir=None, chunksize=5_000_000):
    """Per-barcode UMI sums split by species, streamed from matrix.mtx.
    With exon_dir, also human_exonic / pig_exonic from that (Gene) matrix."""
    features = pd.read_csv(_find(solo_dir, "features.tsv"), sep="\t",
                           header=None, usecols=[0, 1], names=["gene_id", "gene_name"])
    barcodes = pd.read_csv(_find(solo_dir, "barcodes.tsv"), sep="\t",
                           header=None, names=["cb"])["cb"].values

    gid = features["gene_id"].values.astype(str)
    is_h = np.char.startswith(gid, hpfx + "_")
    is_p = np.char.startswith(gid, ppfx + "_")
    if not (is_h.any() and is_p.any()):
        sys.exit(f"features.tsv has no {hpfx}_/{ppfx}_ gene_ids -- was the GTF prefixed?")
    unknown = int((~(is_h | is_p)).sum())
    if unknown:
        sys.stderr.write(f"warning: {unknown} features match neither prefix; ignored\n")

    if gtf:
        mt = mito_gene_ids(gtf, (hpfx, ppfx))
        in_mt = np.isin(gid, list(mt))
    else:
        in_mt = np.zeros(len(gid), bool)
    is_hmt, is_pmt = is_h & in_mt, is_p & in_mt

    # human-only reporter transgenes (matched on gene_id or gene_name)
    gname = features["gene_name"].values.astype(str)
    reporter_genes = list(reporter_genes)
    is_rep = np.isin(gid, reporter_genes) | np.isin(gname, reporter_genes)
    found = set(gid[is_rep]) | set(gname[is_rep])
    missing = [r for r in reporter_genes if r not in found]
    if missing:
        sys.stderr.write(f"warning: reporter gene(s) not in features.tsv: {missing}\n")
    sys.stderr.write(f"reporter features: {list(gid[is_rep])}\n")
    sys.stderr.write(f"features: {is_h.sum()} human ({is_hmt.sum()} MT), "
                     f"{is_p.sum()} pig ({is_pmt.sum()} MT)\n")

    mtx = _find(solo_dir, "matrix.mtx")
    with _open(mtx) as fh:
        for line in fh:
            if not line.startswith("%"):
                ngene, nbc, nnz = (int(x) for x in line.split())
                break
    if ngene != len(gid) or nbc != len(barcodes):
        sys.exit("matrix header does not match features/barcodes length")
    sys.stderr.write(f"matrix: {ngene} x {nbc}, {nnz} nonzero\n")

    masks = {"human": is_h, "pig": is_p, "human_mt": is_hmt, "pig_mt": is_pmt,
             "reporter": is_rep}
    acc = _stream_barcode_sums(mtx, masks, nbc, chunksize)

    if exon_dir:
        # must be the same STARsolo run: identical features and barcodes
        fx = pd.read_csv(_find(exon_dir, "features.tsv"), sep="\t", header=None,
                         usecols=[0])[0].values.astype(str)
        bx = pd.read_csv(_find(exon_dir, "barcodes.tsv"), sep="\t",
                         header=None)[0].values
        if not (np.array_equal(fx, gid) and np.array_equal(bx, barcodes)):
            sys.exit("--exon-solo-dir features/barcodes differ from --solo-dir; "
                     "use the Gene/raw matrix from the same STARsolo run")
        sys.stderr.write("reading exon-only matrix for the intronic check\n")
        ex = _stream_barcode_sums(_find(exon_dir, "matrix.mtx"),
                                  {"human_exonic": is_h, "pig_exonic": is_p},
                                  nbc, chunksize)
        acc.update(ex)

    df = pd.DataFrame(acc)
    df.insert(0, "cb", barcodes)
    return df


# ---------------------------------------------------------------- calling
def classify(df, min_umi, human_purity, pig_purity, min_reporter_umi=1):
    """human: frac_human >= human_purity; pig: frac_pig >= pig_purity."""
    df = df.copy()
    df["total"] = df["human"] + df["pig"]
    with np.errstate(invalid="ignore", divide="ignore"):
        df["frac_human"] = df["human"] / df["total"]
        mt_tot = df["human_mt"] + df["pig_mt"]
        df["frac_human_mt"] = np.where(mt_tot > 0, df["human_mt"] / mt_tot, np.nan)

    call = np.full(len(df), "unassigned", dtype=object)
    keep = df["total"].values >= min_umi
    fh = df["frac_human"].values
    is_h = keep & (fh >= human_purity)
    is_p = keep & (fh <= 1 - pig_purity)
    call[is_h] = "human"
    call[is_p] = "pig"
    call[keep & ~is_h & ~is_p] = "mixed"
    df["call_purity"] = call.copy()

    # reporter transgenes are human-only: any detected reporter => human
    has_rep = df["reporter"].values >= min_reporter_umi
    override = keep & has_rep & (call != "human")
    call[keep & has_rep] = "human"
    df["reporter_override"] = override
    df["call"] = call
    return df


# ------------------------------------------------------- nuclear doublets
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


def fit_size_model(y, N, iters=50):
    """Poisson IRLS for log mu = a + b log N, then moment NB dispersion.
    (Same model as 05_cross_species_excess.py.)"""
    y = np.asarray(y, float)
    X = np.c_[np.ones(len(N)), np.log(N)]
    beta = np.linalg.lstsq(X, np.log(y + 0.5), rcond=None)[0]
    for _ in range(iters):
        eta = X @ beta
        mu = np.exp(eta)
        new = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (mu * (eta + (y - mu) / mu)))
        done = np.max(np.abs(new - beta)) < 1e-8
        beta = new
        if done:
            break
    mu = np.exp(X @ beta)
    phi = max(0.0, float(((y - mu) ** 2 - mu).sum() / (mu ** 2).sum()))
    return beta, mu, phi


def nb_sf(y, mu, phi):
    if phi <= 0:
        return stats.poisson.sf(y - 1, mu)
    r = 1.0 / phi
    return stats.nbinom.sf(y - 1, r, r / (r + mu))


def call_nuclear_doublets(df, excess_fdr=0.05, intronic_fdr=0.05):
    """Flag cells with excess foreign UMIs (vs a cell-size model) and, if
    exon-only counts are present, call those whose foreign reads are
    significantly intronic `nuclear_doublet`.  Returns (df, stats dict)."""
    df = df.copy()
    has_exon = {"human_exonic", "pig_exonic"} <= set(df.columns)
    for col in ("foreign_expected", "foreign_excess_q", "foreign_exonic_frac",
                "intronic_q"):
        df[col] = np.nan
    df["nuclear_doublet"] = False
    out = {"nuclear_check": "on" if has_exon else "off (no --exon-solo-dir)"}
    for sp, fs in (("pig", "human"), ("human", "pig")):
        m = (df["call_purity"] == sp).values
        if m.sum() < 50:
            continue
        y = df.loc[m, fs].values
        N = df.loc[m, "total"].values.astype(float)
        _, mu, phi = fit_size_model(y, N)
        q = bh(nb_sf(y, mu, phi))
        exc = q < excess_fdr
        idx = np.where(m)[0]
        df.iloc[idx, df.columns.get_loc("foreign_expected")] = np.round(mu, 2)
        df.iloc[idx, df.columns.get_loc("foreign_excess_q")] = q
        o = {f"n_{sp}_cells_foreign_excess": int(exc.sum())}
        if has_exon:
            ex = np.minimum(df.loc[m, fs + "_exonic"].values, y)  # exonic <= total
            with np.errstate(invalid="ignore", divide="ignore"):
                df.iloc[idx, df.columns.get_loc("foreign_exonic_frac")] = np.round(ex / y, 4)
            ref = ex[~exc].sum() / max(y[~exc].sum(), 1)
            p = np.full(len(y), np.nan)
            sel = exc & (y > 0)
            p[sel] = stats.binom.cdf(ex[sel], y[sel], ref)
            qi = bh(p)
            nuc = sel & (qi < intronic_fdr)
            df.iloc[idx, df.columns.get_loc("intronic_q")] = qi
            df.iloc[idx[nuc], df.columns.get_loc("nuclear_doublet")] = True
            with np.errstate(invalid="ignore", divide="ignore"):
                o.update({
                    f"{fs}_exonic_frac_ref_in_{sp}_cells": round(float(ref), 4),
                    f"{fs}_exonic_frac_in_excess_{sp}_cells":
                        round(float(ex[exc].sum() / y[exc].sum()), 4) if exc.any() else None,
                    f"n_{sp}_cells_nuclear_doublet": int(nuc.sum())})
        out.update(o)
    df.loc[df["nuclear_doublet"], "call"] = "nuclear_doublet"
    return df, out


def estimate_rho(df, soup_max_umi=100):
    """SoupX contamination fraction per species, using the other species as a
    tracer.  Adds a per-cell `rho_cell` column; returns (df, dict)."""
    soup = df[(df["total"] > 0) & (df["total"] <= soup_max_umi)]
    soup_umi = int(soup["total"].sum())
    f_hs = soup["human"].sum() / soup_umi if soup_umi else np.nan

    # nuclear doublets carry real foreign nuclear RNA, not ambient: leave out
    nd = df["nuclear_doublet"] if "nuclear_doublet" in df else False
    pc = (df["call_purity"] == "pig") & ~nd
    hc = (df["call_purity"] == "human") & ~nd
    with np.errstate(invalid="ignore", divide="ignore"):
        x_p = df.loc[pc, "human"].sum() / df.loc[pc, "total"].sum() if pc.any() else np.nan
        x_h = df.loc[hc, "pig"].sum() / df.loc[hc, "total"].sum() if hc.any() else np.nan
        rho_p = x_p / f_hs
        rho_h = x_h / (1 - f_hs)
        rc = np.full(len(df), np.nan)
        rc[pc.values] = df.loc[pc, "frac_human"].values / f_hs
        rc[hc.values] = (1 - df.loc[hc, "frac_human"].values) / (1 - f_hs)
    df = df.copy()
    df["rho_cell"] = rc

    if f_hs == f_hs and min(f_hs, 1 - f_hs) < 0.01:
        sys.stderr.write("warning: soup is >99% one species; rho for the cells of "
                         "that species is poorly determined\n")
    r = lambda v, n=4: round(float(v), n) if v == v else None
    return df, {
        "soup_max_umi": soup_max_umi,
        "n_soup_barcodes": int(len(soup)), "soup_umis": soup_umi,
        "soup_frac_human": r(f_hs),
        "pooled_human_frac_in_pig_cells": r(x_p),
        "pooled_pig_frac_in_human_cells": r(x_h),
        "rho_pig": r(rho_p, 3), "rho_human": r(rho_h, 3),
        # where a clean cell of each species is expected to sit, given rho
        "expected_frac_human_clean_human_cell": r(1 - rho_h * (1 - f_hs), 3),
        "expected_frac_pig_clean_pig_cell": r(1 - rho_p * f_hs, 3),
    }


def summarise(df, min_umi, human_purity, pig_purity, min_reporter_umi=1):
    n_h = int((df["call"] == "human").sum())
    n_p = int((df["call"] == "pig").sum())
    n_m = int((df["call"] == "mixed").sum())
    n_nd = int((df["call"] == "nuclear_doublet").sum())
    n_cells = n_h + n_p + n_m + n_nd

    # ambient: minority fraction within purity-assigned barcodes (reporter
    # overrides and nuclear doublets excluded)
    nd = df["call"] == "nuclear_doublet"
    amb_h = (1 - df.loc[(df["call_purity"] == "human") & ~nd, "frac_human"]).median()
    amb_p = df.loc[(df["call_purity"] == "pig") & ~nd, "frac_human"].median()

    # reporter detection / overrides
    cells = df["call_purity"] != "unassigned"
    rep_pos = df["reporter"] >= min_reporter_umi
    ov = df["reporter_override"]
    n_ov_pig = int((ov & (df["call_purity"] == "pig")).sum())
    n_ov_mixed = int((ov & (df["call_purity"] == "mixed")).sum())
    ph = df["call_purity"] == "human"
    rep_frac_h = float(rep_pos[ph].mean()) if ph.any() else np.nan

    # total doublet rate from the observed heterotypic rate
    f_h = n_h / (n_h + n_p) if (n_h + n_p) else np.nan
    f_p = 1 - f_h
    obs = n_m / n_cells if n_cells else np.nan
    tot = obs / (2 * f_h * f_p) if (n_cells and f_h and f_p) else np.nan

    # MT-based call as an independent check (on the purity call)
    sub = df[df["call_purity"].isin(["human", "pig"]) & df["frac_human_mt"].notna()]
    if len(sub):
        mt_call = np.where(sub["frac_human_mt"] >= 0.5, "human", "pig")
        concord = float((mt_call == sub["call_purity"].values).mean())
    else:
        concord = None

    return {
        "min_umi": min_umi,
        "human_purity": human_purity, "pig_purity": pig_purity,
        "n_human": n_h, "n_pig": n_p, "n_mixed": n_m,
        "n_nuclear_doublet": n_nd, "n_cells": n_cells,
        "n_unassigned_barcodes": int((df["call"] == "unassigned").sum()),
        "species_fraction_human": round(f_h, 4) if f_h == f_h else None,
        "median_pig_ambient_in_human_cells": round(float(amb_h), 4) if amb_h == amb_h else None,
        "median_human_ambient_in_pig_cells": round(float(amb_p), 4) if amb_p == amb_p else None,
        "observed_heterotypic_doublet_rate": round(float(obs), 4) if obs == obs else None,
        "estimated_total_doublet_rate": round(float(tot), 4) if tot == tot else None,
        "mt_call_concordance": round(concord, 4) if concord is not None else None,
        "min_reporter_umi": min_reporter_umi,
        "n_reporter_positive_cells": int((cells & rep_pos).sum()),
        "frac_purity_human_cells_reporter_positive":
            round(rep_frac_h, 4) if rep_frac_h == rep_frac_h else None,
        "n_reporter_override_from_pig": n_ov_pig,
        "n_reporter_override_from_mixed": n_ov_mixed,
    }


# ---------------------------------------------------------------- figure
def barnyard_plot(df, human_purity, pig_purity, min_umi, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 300, "savefig.bbox": "tight",
    })
    col = {"human": "#2166ac", "pig": "#b2182b", "mixed": "#762a83",
           "nuclear_doublet": "#e08214", "unassigned": "#cccccc"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.1),
                                   gridspec_kw={"width_ratios": [1.25, 1]})

    for k in ("unassigned", "human", "pig", "mixed", "nuclear_doublet"):
        s = df[df["call"] == k]
        if not len(s):
            continue
        ax1.scatter(s["human"] + 1, s["pig"] + 1, s=2 if k != "unassigned" else 1,
                    c=col[k], alpha=0.5 if k != "unassigned" else 0.25,
                    linewidths=0, rasterized=True, label=f"{k} (n={len(s)})")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("human UMIs per barcode + 1")
    ax1.set_ylabel("pig UMIs per barcode + 1")
    ax1.set_title("Two clean arms separate the species")
    # limits from a high quantile, not the max: one extreme droplet otherwise
    # leaves most of the panel empty
    hi = max(1e3, 5 * max(np.percentile(df["human"], 99.9),
                          np.percentile(df["pig"], 99.9)))
    n_off = int(((df["human"] > hi) | (df["pig"] > hi)).sum())
    r = np.array([1, hi])
    # human cut: pig/human = (1-h)/h ; pig cut: pig/human = p/(1-p)
    ax1.plot(r, r * (1 - human_purity) / human_purity, lw=0.6, ls="--", c="0.35", zorder=0)
    ax1.plot(r, r * pig_purity / (1 - pig_purity), lw=0.6, ls="--", c="0.35", zorder=0)
    ax1.set_xlim(1, hi); ax1.set_ylim(1, hi)
    ax1.legend(frameon=False, loc="upper left", markerscale=4,
               handletextpad=0.3, borderpad=0.1)
    note = f"dashed: human >= {human_purity:.0%}, pig >= {pig_purity:.0%}"
    if n_off:
        note += f"\n{n_off} barcode(s) beyond axis"
    ax1.text(0.97, 0.03, note, transform=ax1.transAxes, fontsize=6,
             color="0.35", ha="right", va="bottom")

    cells = df[df["total"] >= min_umi]
    ax2.hist(cells["frac_human"], bins=60, color="0.45", linewidth=0)
    ax2.set_yscale("log")
    ax2.axvline(1 - pig_purity, lw=0.6, ls="--", c="0.35")
    ax2.axvline(human_purity, lw=0.6, ls="--", c="0.35")
    ax2.set_xlabel("fraction of species-assignable UMIs that are human")
    ax2.set_ylabel("barcodes")
    ax2.set_title("Assignment is bimodal, not graded")
    ax2.margins(x=0.02)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo-dir", required=True,
                    help="STARsolo raw matrix dir, e.g. .../Solo.out/GeneFull/raw")
    ap.add_argument("--gtf", default=None, help="combined GTF (for MT genes)")
    ap.add_argument("--human-prefix", default="GRCh38")
    ap.add_argument("--pig-prefix", default="Sscrofa11")
    ap.add_argument("--min-umi", type=int, default=500)
    ap.add_argument("--purity", type=float, default=0.90,
                    help="default for both cuts if --human-purity/--pig-purity not given")
    ap.add_argument("--human-purity", type=float, default=None,
                    help="min fraction human UMIs to call human (default: --purity)")
    ap.add_argument("--pig-purity", type=float, default=None,
                    help="min fraction pig UMIs to call pig (default: --purity)")
    ap.add_argument("--soup-max-umi", type=int, default=100,
                    help="barcodes with 0 < UMIs <= this define the soup (SoupX default 100)")
    ap.add_argument("--reporter-genes", nargs="*",
                    default=["GRCh38_NLuc_P2A_tdTomato", "GRCh38_NeoR"],
                    help="human-only reporter genes; detection forces a human call "
                         "(pass with no values to disable)")
    ap.add_argument("--min-reporter-umi", type=int, default=1,
                    help="summed reporter UMIs needed to force a human call")
    ap.add_argument("--exon-solo-dir", default=None,
                    help="exon-only Gene/raw matrix of the same STARsolo run; "
                         "enables the intronic check and nuclear_doublet calls")
    ap.add_argument("--excess-fdr", type=float, default=0.05,
                    help="BH q for 'more foreign UMIs than ambient predicts'")
    ap.add_argument("--intronic-fdr", type=float, default=0.05,
                    help="BH q for 'foreign reads more intronic than reference'")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    hp = a.purity if a.human_purity is None else a.human_purity
    pp = a.purity if a.pig_purity is None else a.pig_purity
    if not (0.5 <= hp <= 1 and 0.5 <= pp <= 1):
        sys.exit("purity cuts must be between 0.5 and 1")

    os.makedirs(a.outdir, exist_ok=True)
    df = species_sums(a.solo_dir, a.human_prefix, a.pig_prefix, a.gtf,
                      a.reporter_genes, a.exon_solo_dir)
    df = classify(df, a.min_umi, hp, pp, a.min_reporter_umi)
    df, nuc = call_nuclear_doublets(df, a.excess_fdr, a.intronic_fdr)
    df, rho = estimate_rho(df, a.soup_max_umi)
    s = summarise(df, a.min_umi, hp, pp, a.min_reporter_umi)
    s.update(rho)
    s.update(nuc)

    # a cut above where clean cells are expected to sit discards real cells
    for sp, cut, key in (("human", hp, "expected_frac_human_clean_human_cell"),
                         ("pig", pp, "expected_frac_pig_clean_pig_cell")):
        exp = s.get(key)
        if exp is not None and cut > exp - 0.02:
            sys.stderr.write(
                f"warning: --{sp}-purity {cut} is at/above where clean {sp} cells "
                f"are expected (~{exp}); many {sp} cells will be called mixed. "
                f"Check the barnyard histogram and lower the cut.\n")

    df.to_csv(os.path.join(a.outdir, "species_assignment.csv.gz"), index=False)
    for k in ("human", "pig", "mixed", "nuclear_doublet"):
        df.loc[df["call"] == k, "cb"].to_csv(
            os.path.join(a.outdir, f"barcodes_{k}.txt"), index=False, header=False)
    with open(os.path.join(a.outdir, "summary.json"), "w") as fh:
        json.dump(s, fh, indent=2)
    barnyard_plot(df, hp, pp, a.min_umi, os.path.join(a.outdir, "barnyard.png"))

    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
