#!/usr/bin/env python3
"""Extend gene 3' ends in a GTF using read coverage from 3' scRNA-seq alignments.

Non-model annotations (e.g. Ensembl pig) often stop genes at the stop codon with
no 3' UTR.  10x 3' reads pile up just upstream of the poly(A) site, so they land
past the annotated end and STARsolo discards them as intergenic (pig FOXC2 loses
~85-90% of its reads this way).

For each gene we search downstream of its 3' end, on the gene's strand, and
extend it to where sense-strand coverage stops.  The search is capped at
--max-ext and stops --buffer bp before the next gene or exon of another gene on
the same strand, so no new same-strand overlaps are created (the SOX2/SOX2-OT
problem).  Walking out from the 3' end in --bin sized bins, a bin counts as
covered if it has >= max(--min-reads, --rel x the highest bin seen so far); the
walk ends after more than --max-gap bp of uncovered bins.  The new end is --pad
bp past the last covered bin (still inside the window), because STARsolo Gene
only counts reads that fit entirely within a transcript.

--refseq-gff adds a second source: NCBI RefSeq 3' ends (same assembly, e.g.
GCF_000003025.6_Sscrofa11.1_genomic.gff for Sscrofa11.1).  Each Ensembl gene is
matched to an overlapping RefSeq gene on the same strand (same gene name
preferred, else the largest overlap) and the RefSeq gene's 3'-most end becomes a
candidate.  It is clipped at the same neighbor boundary as the coverage search,
but a RefSeq end beyond --max-ext is dropped rather than clipped, since it is
usually a read-through model or a RefSeq gene spanning several Ensembl genes
(report limit "refseq_too_far").  The new end is whichever candidate reaches
further.  RefSeq carries 3' UTRs for many genes
Ensembl leaves CDS-only, including genes too lowly expressed in your data for
coverage to find; coverage catches genes RefSeq also stops short (pig SOX2 is
CDS-only in both, and RefSeq's FOXC2 UTR misses ~1/3 of the reads).  Give
--refseq-gff without alignments for RefSeq-only extension.  RefSeq sequence
names are translated to chromosome names (1..18, X, Y, MT) from the GFF's region
lines; genes on unplaced scaffolds are left to coverage.  Only 3' ends are taken
from RefSeq; no RefSeq genes or transcripts are added and Ensembl IDs are kept.

Only gene, transcript, exon, three_prime_utr and UTR lines that end at the
gene's 3'-most coordinate are moved; CDS / stop_codon are left alone, as are
isoforms that end further upstream.

Reads: primary, mapped, uniquely mapped (NH:i:1) unless --include-multimappers.
Strand follows the read (STARsolo --soloStrand Forward: R2 is sense).  Coverage
is summed over all inputs, so --min-reads applies to the total; raise it if you
pass many lanes.  Inputs can be SAM or BAM (BAM needs samtools on PATH, e.g.
`conda activate star_env`).

Scanning full SAMs is slow, so --coverage-out saves the binned counts and
--coverage-in reuses them to retune --min-reads/--rel/--max-gap/--buffer
without rescanning (same GTF, same --bin, and --max-ext no larger than the scan).

Report columns: refseq_3p is the matched RefSeq gene's 3' end before clipping
(NA if unmatched), coverage_3p the coverage-method end (NA if none), source
which one set new_3p ("both" when they agree within --agree bp), and cov_in_ext
the summed per-bin read count over the coverage extension; a read touching
several bins counts in each, so it overstates the read number.

example (pig, both lanes of one sample, plus RefSeq; ~20-40 min per 67 GB SAM):
  python3 extend_3prime_gtf.py \\
      --gtf Sus_scrofa.Sscrofa11.1.109.chr.filtered.gtf \\
      --refseq-gff ../ncbi/GCF_000003025.6_Sscrofa11.1_genomic.gff \\
      --out Sus_scrofa.Sscrofa11.1.109.chr.filtered.3pext.gtf \\
      --report 3pext_report.tsv --coverage-out 3pext_coverage.tsv.gz \\
      shared_mount/pTLS_D4_Mgel_F20CH5_lane{1,2}/Aligned.out.sam

usage: extend_3prime_gtf.py --gtf in.gtf --out out.gtf --report report.tsv
                            [--refseq-gff refseq.gff]
                            [--coverage-out cov.tsv.gz | --coverage-in cov.tsv.gz]
                            [aligned.sam|bam ...]
"""
import argparse, bisect, gzip, re, subprocess, sys
from collections import defaultdict

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
NH = re.compile(r"\tNH:i:(\d+)")
EXTEND_FEATURES = {"gene", "transcript", "exon", "three_prime_utr", "UTR"}

def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else ""

def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def read_gtf(path):
    genes, utr3 = {}, set()
    bounds = defaultdict(lambda: ([], []))   # (chrom, strand) -> (starts, ends) of genes and exons
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon", "three_prime_utr"):
                continue
            gid = attr(f[8], "gene_id")
            if f[2] == "three_prime_utr":
                utr3.add(gid)
                continue
            s, e = int(f[3]), int(f[4])
            bounds[(f[0], f[6])][0].append(s)
            bounds[(f[0], f[6])][1].append(e)
            if f[2] == "gene":
                genes[gid] = dict(chrom=f[0], strand=f[6], start=s, end=e,
                                  name=attr(f[8], "gene_name"), biotype=attr(f[8], "gene_biotype"))
    for starts, ends in bounds.values():
        starts.sort()
        ends.sort()
    return genes, utr3, bounds

def read_refseq(path):
    """(chrom, strand) -> [(start, end, name)] for RefSeq genes on assembled chromosomes."""
    chrom, genes = {}, defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "region":
                m = re.search(r"(?:^|;)chromosome=([^;]+)", f[8])
                if "genome=mitochondrion" in f[8]:
                    chrom[f[0]] = "MT"
                elif m and "genome=chromosome" in f[8]:
                    chrom[f[0]] = m.group(1)
            elif f[2] == "gene" and f[0] in chrom:
                m = re.search(r"(?:^|;)Name=([^;]+)", f[8])
                genes[(chrom[f[0]], f[6])].append((int(f[3]), int(f[4]), m.group(1) if m else ""))
    return genes

def match_refseq(genes, refseq):
    """gene_id -> RefSeq 3' end of the overlapping same-strand RefSeq gene,
    preferring one with the same gene name, else the largest overlap."""
    hits = {}
    for gid, g in genes.items():
        best = None
        for s, e, name in refseq.get((g["chrom"], g["strand"]), ()):
            ov = min(e, g["end"]) - max(s, g["start"]) + 1
            if ov <= 0:
                continue
            score = (bool(g["name"]) and name == g["name"], ov)
            if best is None or score > best[0]:
                best = (score, e if g["strand"] == "+" else s)
        if best:
            hits[gid] = best[1]
    return hits

def build_windows(genes, bounds, max_ext, buffer):
    """gene_id -> (lo, hi, limit): the genomic interval searched past the 3' end.
    The gene's own features all lie inside [start, end], so the nearest boundary
    beyond the 3' end always belongs to another gene."""
    win = {}
    for gid, g in genes.items():
        starts, ends = bounds[(g["chrom"], g["strand"])]
        if g["strand"] == "+":
            hi, why = g["end"] + max_ext, "max_ext"
            i = bisect.bisect_right(starts, g["end"])
            if i < len(starts) and starts[i] - buffer - 1 < hi:
                hi, why = starts[i] - buffer - 1, "neighbor"
            if hi > g["end"]:
                win[gid] = (g["end"] + 1, hi, why)
        elif g["strand"] == "-":
            lo, why = max(1, g["start"] - max_ext), "max_ext"
            i = bisect.bisect_left(ends, g["start"]) - 1
            if i >= 0 and ends[i] + buffer + 1 > lo:
                lo, why = ends[i] + buffer + 1, "neighbor"
            if lo < g["start"]:
                win[gid] = (lo, g["start"] - 1, why)
    return win

def wanted_bins(genes, win, binsize):
    wanted = defaultdict(set)
    for gid, (lo, hi, _) in win.items():
        g = genes[gid]
        wanted[(g["chrom"], g["strand"])].update(range(lo // binsize, hi // binsize + 1))
    return wanted

def alignments(path):
    if path.endswith(".bam"):
        proc = subprocess.Popen(["samtools", "view", "-F", "0x904", path],
                                stdout=subprocess.PIPE, universal_newlines=True)
        yield from proc.stdout
        if proc.wait() != 0:
            sys.exit(f"samtools view failed on {path}")
    else:
        with opener(path) as fh:
            yield from fh

def scan(paths, wanted, binsize, unique_only):
    cov = defaultdict(int)   # (chrom, strand, bin) -> reads
    for path in paths:
        n = kept = 0
        for line in alignments(path):
            if line[0] == "@":
                continue
            n += 1
            if n % 20000000 == 0:
                sys.stderr.write(f"  {path}: {n:,} alignments, {kept:,} counted\n")
            f = line.split("\t", 6)
            flag = int(f[1])
            if flag & 0x904:   # unmapped, secondary, supplementary
                continue
            strand = "-" if flag & 16 else "+"
            w = wanted.get((f[2], strand))
            if w is None:
                continue
            pos, cigar = int(f[3]), f[5]
            # cheap pre-filter; spliced reads can start far upstream of a window
            if "N" not in cigar and pos // binsize not in w and (pos + 150) // binsize not in w:
                continue
            if unique_only:
                m = NH.search(line)
                if m and m.group(1) != "1":
                    continue
            bins, ref = set(), pos
            for length, op in CIGAR.findall(cigar):
                length = int(length)
                if op in "M=XD":
                    bins.update(range(ref // binsize, (ref + length - 1) // binsize + 1))
                    ref += length
                elif op == "N":
                    ref += length
            hit = False
            for b in bins:
                if b in w:
                    cov[(f[2], strand, b)] += 1
                    hit = True
            kept += hit
        sys.stderr.write(f"  {path}: {n:,} alignments, {kept:,} counted\n")
    return cov

def save_coverage(cov, path, binsize, max_ext):
    with gzip.open(path, "wt") as out:
        out.write(f"#bin={binsize}\tmax_ext={max_ext}\n")
        for (chrom, strand, b), c in sorted(cov.items()):
            out.write(f"{chrom}\t{strand}\t{b}\t{c}\n")

def load_coverage(path, binsize, max_ext):
    cov = defaultdict(int)
    with gzip.open(path, "rt") as fh:
        meta = dict(kv.split("=") for kv in fh.readline().lstrip("#").split())
        if int(meta["bin"]) != binsize:
            sys.exit(f"{path} was scanned with --bin {meta['bin']}, not {binsize}")
        if max_ext > int(meta["max_ext"]):
            sys.exit(f"{path} was scanned with --max-ext {meta['max_ext']}; rescan to search further")
        for line in fh:
            chrom, strand, b, c = line.rstrip("\n").split("\t")
            cov[(chrom, strand, int(b))] = int(c)
    return cov

def extend(genes, win, cov, binsize, min_reads, rel, max_gap, pad):
    """gene_id -> (new 3' coordinate, coverage in extension, what stopped it)"""
    ext = {}
    for gid, (lo, hi, why) in win.items():
        g = genes[gid]
        plus = g["strand"] == "+"
        first, last_bin, step = (lo // binsize, hi // binsize, 1) if plus else (hi // binsize, lo // binsize, -1)
        last = None
        gap = reads = pending = peak = 0
        for b in range(first, last_bin + step, step):
            c = cov.get((g["chrom"], g["strand"], b), 0)
            if c >= max(min_reads, rel * peak):
                last, gap, peak = b, 0, max(peak, c)
                reads += pending + c
                pending = 0
            else:
                gap += binsize
                pending += c
                if gap > max_gap:
                    break
        if last is None:
            continue
        new = min(hi, (last + 1) * binsize - 1 + pad) if plus else max(lo, last * binsize - pad)
        ext[gid] = (new, reads, why if new == (hi if plus else lo) else "coverage")
    return ext

def combine(genes, win, cov_ext, refseq_3p, agree, max_ext):
    """gene_id -> (new 3' coordinate, coverage in extension, what stopped it, source),
    taking whichever of the coverage and RefSeq ends reaches further downstream.
    A RefSeq end past a neighbor is clipped there, but one past --max-ext is
    dropped (usually a read-through model or a RefSeq gene spanning several
    Ensembl genes), and returned in the second value."""
    ext, too_far = {}, set()
    for gid, (lo, hi, why) in win.items():
        g = genes[gid]
        plus = g["strand"] == "+"
        further = (lambda x, y: x > y) if plus else (lambda x, y: x < y)
        old = g["end"] if plus else g["start"]
        cands = []
        if gid in cov_ext:
            new, reads, limit = cov_ext[gid]
            cands.append((new, reads, limit, "coverage"))
        r = refseq_3p.get(gid)
        if r is not None and further(r, old):
            clipped = min(r, hi) if plus else max(r, lo)
            if abs(r - old) > max_ext:
                too_far.add(gid)
            elif clipped == r:
                cands.append((r, cov_ext.get(gid, (0, 0))[1], "refseq", "refseq"))
            else:
                cands.append((clipped, cov_ext.get(gid, (0, 0))[1], why, "refseq"))
        if not cands:
            continue
        best = max(cands, key=lambda c: c[0]) if plus else min(cands, key=lambda c: c[0])
        if len(cands) == 2 and abs(cands[0][0] - cands[1][0]) <= agree:
            best = best[:3] + ("both",)
        ext[gid] = best
    return ext, too_far

def write_gtf(src, dst, genes, ext):
    moved = 0
    with opener(src) as fh, open(dst, "w") as out:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(f) >= 9 and f[2] in EXTEND_FEATURES:
                gid = attr(f[8], "gene_id")
                if gid in ext:
                    g = genes[gid]
                    if g["strand"] == "+" and int(f[4]) == g["end"]:
                        f[4] = str(ext[gid][0])
                    elif g["strand"] == "-" and int(f[3]) == g["start"]:
                        f[3] = str(ext[gid][0])
                    else:
                        out.write(line)
                        continue
                    line = "\t".join(f) + "\n"
                    moved += 1
            out.write(line)
    return moved

def write_report(path, genes, win, ext, cov_ext, refseq_3p, too_far, utr3):
    with open(path, "w") as out:
        out.write("gene_id\tgene_name\tbiotype\tchrom\tstrand\thad_3utr\told_3p\trefseq_3p\tcoverage_3p\t"
                  "new_3p\text_bp\tsource\tcov_in_ext\tlimit\twindow_bp\n")
        for gid, g in sorted(genes.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start"])):
            old = g["end"] if g["strand"] == "+" else g["start"]
            lo, hi, why = win.get(gid, (0, -1, "no_room"))
            none_why = "refseq_too_far" if gid in too_far else "no_evidence" if gid in win else why
            new, reads, limit, source = ext.get(gid, (old, 0, none_why, "none"))
            if gid in too_far and gid in ext:
                limit += ";refseq_too_far"
            out.write(f"{gid}\t{g['name']}\t{g['biotype']}\t{g['chrom']}\t{g['strand']}\t"
                      f"{'yes' if gid in utr3 else 'no'}\t{old}\t{refseq_3p.get(gid, 'NA')}\t"
                      f"{cov_ext[gid][0] if gid in cov_ext else 'NA'}\t{new}\t{abs(new - old)}\t"
                      f"{source}\t{reads}\t{limit}\t{hi - lo + 1}\n")

def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("alignments", nargs="*", help="SAM/BAM files (omit with --coverage-in)")
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True, help="extended GTF")
    p.add_argument("--report", required=True, help="per-gene TSV")
    p.add_argument("--coverage-out", help="save binned coverage (.tsv.gz)")
    p.add_argument("--coverage-in", help="reuse binned coverage instead of scanning")
    p.add_argument("--refseq-gff", help="NCBI RefSeq GFF3 for the same assembly; its 3' ends are a second source")
    p.add_argument("--agree", type=int, default=100, help="bp within which RefSeq and coverage ends count as agreeing [%(default)s]")
    p.add_argument("--max-ext", type=int, default=5000, help="max bp added past the 3' end [%(default)s]")
    p.add_argument("--buffer", type=int, default=100, help="bp kept clear before the next same-strand gene/exon [%(default)s]")
    p.add_argument("--bin", type=int, default=50, help="bin size in bp [%(default)s]")
    p.add_argument("--min-reads", type=int, default=10, help="min reads for a bin to count as covered [%(default)s]")
    p.add_argument("--rel", type=float, default=0.05, help="...and min fraction of the highest bin so far [%(default)s]")
    p.add_argument("--max-gap", type=int, default=300, help="bp of uncovered bins tolerated before stopping [%(default)s]")
    p.add_argument("--pad", type=int, default=100, help="bp added past the last covered bin [%(default)s]")
    p.add_argument("--include-multimappers", action="store_true")
    a = p.parse_args()
    if a.alignments and a.coverage_in:
        p.error("give alignment files or --coverage-in, not both")
    if not (a.alignments or a.coverage_in or a.refseq_gff):
        p.error("give alignment files, --coverage-in and/or --refseq-gff")

    genes, utr3, bounds = read_gtf(a.gtf)
    win = build_windows(genes, bounds, a.max_ext, a.buffer)
    sys.stderr.write(f"{len(genes):,} genes, {len(win):,} with room to extend\n")
    refseq_3p = {}
    if a.refseq_gff:
        refseq_3p = match_refseq(genes, read_refseq(a.refseq_gff))
        sys.stderr.write(f"{len(refseq_3p):,} genes matched to a RefSeq gene\n")
    cov_ext = {}
    if a.alignments or a.coverage_in:
        if a.coverage_in:
            cov = load_coverage(a.coverage_in, a.bin, a.max_ext)
        else:
            cov = scan(a.alignments, wanted_bins(genes, win, a.bin), a.bin, not a.include_multimappers)
            if a.coverage_out:
                save_coverage(cov, a.coverage_out, a.bin, a.max_ext)
        cov_ext = extend(genes, win, cov, a.bin, a.min_reads, a.rel, a.max_gap, a.pad)
    ext, too_far = combine(genes, win, cov_ext, refseq_3p, a.agree, a.max_ext)
    moved = write_gtf(a.gtf, a.out, genes, ext)
    write_report(a.report, genes, win, ext, cov_ext, refseq_3p, too_far, utr3)
    if a.refseq_gff:
        sys.stderr.write(f"{len(too_far):,} RefSeq ends dropped for lying beyond --max-ext\n")

    lens = sorted(abs(e[0] - (genes[g]["end"] if genes[g]["strand"] == "+" else genes[g]["start"])) for g, e in ext.items())
    limits, sources = defaultdict(int), defaultdict(int)
    for e in ext.values():
        limits[e[2]] += 1
        sources[e[3]] += 1
    sys.stderr.write(f"{len(ext):,} genes extended ({moved:,} GTF lines moved)")
    if lens:
        sys.stderr.write(f"; median {lens[len(lens) // 2]} bp, 90th pct {lens[int(0.9 * len(lens))]} bp; "
                         f"source: {dict(sources)}; stopped by: {dict(limits)}")
    sys.stderr.write(f"\nwrote {a.out} and {a.report}\n")

if __name__ == "__main__":
    main()
