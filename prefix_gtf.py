#!/usr/bin/env python3
"""Prefix seqnames and every identifier attribute in a GTF.

STAR counts on gene_id, but features.tsv carries gene_name -- if you prefix only
the contigs, human and pig GAPDH both appear as "GAPDH" and Seurat/scanpy
make.unique() them into GAPDH / GAPDH.1, after which you cannot tell which is
which.  So we prefix the IDs and the names.

usage: prefix_gtf.py <in.gtf[.gz]> <prefix> <out.gtf>
"""
import gzip, re, sys

ATTRS = ("gene_id", "transcript_id", "exon_id", "protein_id", "ccds_id",
         "gene_name", "transcript_name", "havana_gene", "havana_transcript")

def main(src, pfx, dst):
    pat = re.compile(r'(\b(?:%s)\s+")' % "|".join(ATTRS))
    op = gzip.open if src.endswith(".gz") else open
    n = 0
    with op(src, "rt") as fh, open(dst, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            f[0] = pfx + "_" + f[0]
            f[8] = pat.sub(lambda m: m.group(1) + pfx + "_", f[8])
            out.write("\t".join(f) + "\n")
            n += 1
    sys.stderr.write(f"{dst}: {n} feature lines, prefix {pfx}_\n")

if __name__ == "__main__":
    main(*sys.argv[1:4])
