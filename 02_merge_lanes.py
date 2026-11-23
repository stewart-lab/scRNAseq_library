#!/usr/bin/env python3
"""
02 -- merge STARsolo outputs of several sequencing lanes of ONE 10x library into
a single sample, without re-aligning.

Lanes of one library (FASTQs ..._S3_L001_..., ..._S3_L002_...) hold the SAME
cells: the barcoded cDNA library is split across flow-cell lanes, so every cell
barcode appears in every lane.  Treating lanes as separate samples counts each
cell twice downstream.  The proper fix is one STARsolo run over all lanes
(COMBINE_LANES: true in config.json); this script is the no-realignment shortcut.

What it does, per Solo.out feature (Gene, GeneFull, GeneFull_Ex50pAS, ...):
  raw/       lane matrices summed barcode-by-gene (identical features and
             barcode whitelist are required and checked)
  filtered/  the summed raw matrix restricted to the UNION of the lanes'
             STARsolo-filtered barcodes

Caveat: a molecule (cell barcode + UMI + gene) sequenced in more than one lane
is counted once per lane.  On the ChimeraTLS library this inflated counts by
~3.2%, nearly uniformly across genes (2.9-3.2% from low to high expression)
-- harmless for species calls, rho, SoupX, clustering and DE, but UMI totals
(and so --min-umi in 03) read ~3% high.  It grows with sequencing saturation;
an exact count needs a joint STARsolo run.

Output layout matches a COMBINE_LANES alignment, <out>/Solo.out/<feature>/
{raw,filtered}/, so 03/04/05 and the Seurat pipeline read it unchanged.  For
run_seuratv5.sh, put the output in shared_mount/<NAME>/ (the default: the lane
folder name without its _lane<i> suffix) and set "COMBINE_LANES": true for that
sample.  BAMs are not merged.

usage
-----
  python3 02_merge_lanes.py \
      --lanes shared_mount/ChimeraTLS_20H80P_mapped2pigonly_lane1 \
              shared_mount/ChimeraTLS_20H80P_mapped2pigonly_lane2 \
      [--features Gene GeneFull_Ex50pAS] [--em] [--out shared_mount/ChimeraTLS_20H80P_mapped2pigonly]
"""
import argparse, gzip, json, os, re, shutil, sys, time
import numpy as np
import pandas as pd
from scipy import sparse

COUNT_FEATURES = ("Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron")


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _find(d, stem):
    for c in (stem, f"{stem}.gz"):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    return None


def read_lines(path):
    with _open(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


def read_mtx(path, chunksize=10_000_000):
    """Read a MatrixMarket coordinate file into a CSC matrix (genes x barcodes)."""
    n_head, field = 0, "integer"
    with _open(path) as fh:
        for line in fh:
            n_head += 1
            if line.startswith("%%MatrixMarket"):
                field = line.split()[3]
            if not line.startswith("%"):
                ngene, nbc, nnz = (int(x) for x in line.split())
                break
    vdtype = np.float64 if field == "real" else np.int64
    g, c, v = [], [], []
    # skiprows covers the header AND the dimension line (not data)
    for ch in pd.read_csv(path, sep=r"\s+", skiprows=n_head, header=None,
                          names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": vdtype},
                          chunksize=chunksize, engine="c"):
        g.append(ch["g"].values - 1)
        c.append(ch["c"].values - 1)
        v.append(ch["v"].values)
    g, c, v = (np.concatenate(x) if x else np.array([], dt)
               for x, dt in ((g, np.int32), (c, np.int32), (v, vdtype)))
    if len(v) != nnz:
        sys.exit(f"{path}: header says {nnz} entries, read {len(v)}")
    m = sparse.coo_matrix((v, (g, c)), shape=(ngene, nbc)).tocsc()
    m.sum_duplicates()
    return m, field


def write_mtx(path, m, field, chunk=5_000_000):
    """Write a CSC matrix as MatrixMarket coordinate, column-major like STARsolo."""
    m = m.tocsc()
    m.sort_indices()
    m.eliminate_zeros()
    cols = np.repeat(np.arange(m.shape[1], dtype=np.int64) + 1, np.diff(m.indptr))
    rows = m.indices.astype(np.int64) + 1
    vals = m.data
    with open(path, "w") as fh:
        fh.write(f"%%MatrixMarket matrix coordinate {field} general\n%\n")
        fh.write(f"{m.shape[0]} {m.shape[1]} {m.nnz}\n")
        for s in range(0, m.nnz, chunk):
            e = s + chunk
            df = pd.DataFrame({"g": rows[s:e], "c": cols[s:e], "v": vals[s:e]})
            df.to_csv(fh, sep=" ", header=False, index=False,
                      float_format="%.6g" if field == "real" else None)


def merge_feature(lanes, feature, out_root, with_em):
    dirs = [os.path.join(l, "Solo.out", feature) for l in lanes]
    raws = [os.path.join(d, "raw") for d in dirs]
    for r in raws:
        if not _find(r, "matrix.mtx"):
            sys.exit(f"missing {r}/matrix.mtx")
    feats = [read_lines(_find(r, "features.tsv")) for r in raws]
    bcs = [read_lines(_find(r, "barcodes.tsv")) for r in raws]
    for i in range(1, len(lanes)):
        if feats[i] != feats[0]:
            sys.exit(f"{feature}: features.tsv differs between {lanes[0]} and {lanes[i]} "
                     "(different index/GTF?) -- cannot sum")
        if bcs[i] != bcs[0]:
            sys.exit(f"{feature}: raw barcodes.tsv differs between {lanes[0]} and {lanes[i]} "
                     "(different whitelist/chemistry?) -- cannot sum")

    out_raw = os.path.join(out_root, "Solo.out", feature, "raw")
    out_filt = os.path.join(out_root, "Solo.out", feature, "filtered")
    os.makedirs(out_raw, exist_ok=True)
    os.makedirs(out_filt, exist_ok=True)

    info = {"feature": feature}
    mats = ["matrix.mtx"] + (["UniqueAndMult-EM.mtx"] if with_em else [])
    total = None
    for mname in mats:
        t0 = time.time()
        acc, field, per_lane = None, None, []
        for r in raws:
            p = _find(r, mname)
            if p is None:
                sys.exit(f"missing {r}/{mname}")
            m, field = read_mtx(p)
            per_lane.append(float(m.sum()))
            acc = m if acc is None else acc + m
            del m
        write_mtx(os.path.join(out_raw, mname), acc, field)
        info[mname] = {"umis_per_lane": per_lane, "umis_merged": float(acc.sum()),
                       "nnz_merged": int(acc.nnz)}
        sys.stderr.write(f"  {feature}/{mname}: {acc.shape[0]} x {acc.shape[1]}, "
                         f"{acc.nnz} nonzero, {time.time() - t0:.0f}s\n")
        if mname == "matrix.mtx":
            total = acc
        else:
            del acc
    for fn, lines in (("features.tsv", feats[0]), ("barcodes.tsv", bcs[0])):
        with open(os.path.join(out_raw, fn), "w") as fh:
            fh.write("\n".join(lines) + "\n")

    # filtered = union of the lanes' STARsolo cell calls, columns of the summed raw
    union = set()
    for d in dirs:
        fb = _find(os.path.join(d, "filtered"), "barcodes.tsv")
        if fb:
            lane_cells = set(read_lines(fb))
            info.setdefault("filtered_cells_per_lane", []).append(len(lane_cells))
            union |= lane_cells
    if union:
        keep = np.array([b in union for b in bcs[0]])
        idx = np.where(keep)[0]
        write_mtx(os.path.join(out_filt, "matrix.mtx"), total[:, idx], "integer")
        with open(os.path.join(out_filt, "barcodes.tsv"), "w") as fh:
            fh.write("\n".join(bcs[0][i] for i in idx) + "\n")
        with open(os.path.join(out_filt, "features.tsv"), "w") as fh:
            fh.write("\n".join(feats[0]) + "\n")
        info["filtered_cells_union"] = int(len(idx))
    else:
        sys.stderr.write(f"  {feature}: no filtered/ barcodes in the lanes; filtered/ left empty\n")
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--lanes", nargs="+", required=True,
                    help="STARsolo output folders (each containing Solo.out/) of one library")
    ap.add_argument("--features", nargs="*", default=None,
                    help=f"Solo.out features to merge (default: all of {', '.join(COUNT_FEATURES)} present)")
    ap.add_argument("--em", action="store_true",
                    help="also sum UniqueAndMult-EM.mtx (multimapper matrices)")
    ap.add_argument("--out", default=None,
                    help="output folder (default: first lane folder without its _lane<i> suffix)")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output folder")
    a = ap.parse_args()

    lanes = [os.path.normpath(l) for l in a.lanes]
    if len(lanes) < 2:
        sys.exit("give at least two lane folders")
    for l in lanes:
        if not os.path.isdir(os.path.join(l, "Solo.out")):
            sys.exit(f"{l}: no Solo.out/ folder")
    out = a.out
    if out is None:
        base = os.path.basename(lanes[0])
        stem = re.sub(r"_lane\d+$", "", base)
        if stem == base:
            sys.exit("cannot derive --out: first lane folder has no _lane<i> suffix")
        out = os.path.join(os.path.dirname(lanes[0]), stem)
    out = os.path.normpath(out)
    if out in lanes:
        sys.exit("--out must differ from the lane folders")
    if os.path.exists(os.path.join(out, "Solo.out")) and not a.force:
        sys.exit(f"{out}/Solo.out exists; use --force to overwrite")

    feats = a.features
    if not feats:
        feats = [f for f in COUNT_FEATURES
                 if all(os.path.isdir(os.path.join(l, "Solo.out", f, "raw")) for l in lanes)]
        if not feats:
            sys.exit("no count features found in all lanes")
    bad = [f for f in feats if f not in COUNT_FEATURES]
    if bad:
        sys.exit(f"only gene-count features can be summed, not {bad}")

    sys.stderr.write(f"merging {len(lanes)} lanes -> {out}\n  features: {', '.join(feats)}\n")
    infos = [merge_feature(lanes, f, out, a.em) for f in feats]
    meta = {"lanes": lanes, "out": out, "features": infos,
            "note": "lane matrices summed; molecules read in several lanes are counted once "
                    "per lane (~3% over-count at ~35% sequencing saturation)",
            "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(os.path.join(out, "merge_lanes_info.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for i in infos:
        m = i["matrix.mtx"]
        sys.stderr.write(f"{i['feature']}: UMIs per lane {[int(x) for x in m['umis_per_lane']]} "
                         f"-> merged {int(m['umis_merged'])}; filtered cells per lane "
                         f"{i.get('filtered_cells_per_lane')} -> union {i.get('filtered_cells_union')}\n")
    sys.stderr.write(f"wrote {out}/Solo.out/ and merge_lanes_info.json\n")


if __name__ == "__main__":
    main()
