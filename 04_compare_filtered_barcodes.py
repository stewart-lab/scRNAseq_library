#!/usr/bin/env python3
"""
04 -- compare species calls (03_classify_species.py) to STARsolo's filtered
(cell-called) barcodes.

For each call category (human / pig / mixed / unassigned) reports how many
barcodes are also in the filtered set, and writes:

  overlap_summary.tsv              per-category counts + UMI medians
  filtered_not_called.tsv          filtered barcodes that are NOT human/pig/mixed,
                                   with their species UMI counts (mostly
                                   < --min-umi, i.e. 'unassigned')
  filtered_not_called.txt          same, barcodes only
  called_not_filtered_<call>.txt   human/pig/mixed barcodes missing from filtered
  filtered_species_calls.tsv       every filtered barcode with its species call

usage
-----
  python3 04_compare_filtered_barcodes.py \
      --classified-dir OUT/classified_cells_0.9 \
      --filtered-dir   OUT/Solo.out/GeneFull/filtered \
      --outdir         OUT/classified_cells_0.9/filtered_comparison
"""
import argparse, gzip, json, os, sys
import pandas as pd

CALLS = ("human", "pig", "mixed")


def read_barcodes(path):
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        path += ".gz"
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [ln.strip().split("\t")[0] for ln in fh if ln.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classified-dir", required=True,
                    help="03_classify_species.py --outdir")
    ap.add_argument("--filtered-dir", required=True,
                    help="STARsolo filtered matrix dir, e.g. .../Solo.out/GeneFull/filtered")
    ap.add_argument("--outdir", default=None,
                    help="default: <classified-dir>/filtered_comparison")
    a = ap.parse_args()
    outdir = a.outdir or os.path.join(a.classified_dir, "filtered_comparison")
    os.makedirs(outdir, exist_ok=True)

    df = pd.read_csv(os.path.join(a.classified_dir, "species_assignment.csv.gz"))
    filt = pd.Index(read_barcodes(os.path.join(a.filtered_dir, "barcodes.tsv")))
    if filt.has_duplicates:
        sys.exit("duplicate barcodes in filtered barcodes.tsv")

    min_umi = None
    sj = os.path.join(a.classified_dir, "summary.json")
    if os.path.exists(sj):
        min_umi = json.load(open(sj)).get("min_umi")

    df["in_filtered"] = df["cb"].isin(filt)
    missing = filt.difference(df["cb"])
    if len(missing):
        print(f"WARNING: {len(missing)} filtered barcodes absent from "
              "species_assignment (different whitelist / matrix?)", file=sys.stderr)

    # per-category overlap
    rows = []
    for k in CALLS + ("unassigned",):
        s = df[df["call"] == k]
        inf = s[s["in_filtered"]]
        rows.append({
            "call": k,
            "n_called": len(s),
            "n_in_filtered": len(inf),
            "pct_called_in_filtered": round(100 * len(inf) / len(s), 2) if len(s) else None,
            "pct_of_filtered": round(100 * len(inf) / len(filt), 2),
            "median_total_umi_in_filtered": inf["total"].median() if len(inf) else None,
            "median_total_umi_not_filtered": s.loc[~s["in_filtered"], "total"].median()
                                             if (~s["in_filtered"]).any() else None,
        })
    rows.append({"call": "absent_from_assignment", "n_called": None,
                 "n_in_filtered": len(missing),
                 "pct_of_filtered": round(100 * len(missing) / len(filt), 2)})
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(outdir, "overlap_summary.tsv"), sep="\t", index=False)

    # filtered barcodes not called human/pig/mixed
    nc = df[df["in_filtered"] & ~df["call"].isin(CALLS)]
    nc = nc.drop(columns="in_filtered").sort_values("total", ascending=False)
    if len(missing):
        nc = pd.concat([nc, pd.DataFrame({"cb": missing, "call": "absent_from_assignment"})])
    nc.to_csv(os.path.join(outdir, "filtered_not_called.tsv"), sep="\t", index=False)
    nc["cb"].to_csv(os.path.join(outdir, "filtered_not_called.txt"),
                    index=False, header=False)

    # called barcodes missing from filtered
    for k in CALLS:
        df.loc[(df["call"] == k) & ~df["in_filtered"], "cb"].to_csv(
            os.path.join(outdir, f"called_not_filtered_{k}.txt"),
            index=False, header=False)

    # every filtered barcode with its call
    df[df["in_filtered"]].drop(columns="in_filtered").to_csv(
        os.path.join(outdir, "filtered_species_calls.tsv"), sep="\t", index=False)

    # console report
    pd.set_option("display.width", 200)
    print(f"filtered barcodes: {len(filt)}   called human/pig/mixed: "
          f"{int(df['call'].isin(CALLS).sum())}")
    print(summ.to_string(index=False))
    if len(nc) and "total" in nc:
        t = nc["total"].dropna()
        print(f"\nfiltered-but-not-called: n={len(nc)}  total UMI "
              f"median={t.median():.0f}  max={t.max():.0f}"
              + (f"  (n >= min_umi {min_umi}: {int((t >= min_umi).sum())})"
                 if min_umi is not None else ""))
        q = t.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(0)
        print("  quantiles:", dict(q))
    print(f"\nwrote {outdir}")


if __name__ == "__main__":
    main()
