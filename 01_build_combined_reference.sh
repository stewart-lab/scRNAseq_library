#!/usr/bin/env bash
# ============================================================================
# 01 -- build the combined human (GRCh38 + reporters) + pig (Sscrofa11.1)
#       reference and its STAR index, as used for the ChimeraTLS analysis.
#
# Reproduces /w5home/bmoore/genomes/human-pig_chimera/ as built 2026-11-12
# (index star_solo/: 71,329 genes = 39,548 GRCh38_ + 31,781 Sscrofa11_).
# See README_chimeric_analysis.md for the whole pipeline.
#
# Inputs (already on disk, nothing is downloaded):
#   human FASTA  Ensembl GRCh38 dna_sm primary assembly + the reporter contig
#                GRCh38_MEOX1_reporter (NLuc-P2A-tdTomato, SV40 pA, FRT-NeoR)
#   human GTF    Ensembl 108, filtered, + reporter genes NLuc_P2A_tdTomato, NeoR
#   pig FASTA    Ensembl Sscrofa11.1 primary assembly, all sequences
#                (613: chromosomes 1-18, X, Y, MT + unplaced scaffolds)
#   pig GTF      Ensembl 109 Sscrofa11.1 .chr, filtered, with 3' ends extended
#                by extend_3prime_gtf.py (step 0; report: 3pext_report.tsv)
#
# The pig FASTA and GTF MUST be the same assembly with the same sequence
# names.  An earlier build paired the Berkshire_pig_v1 contig assembly
# (LUXW... contigs) with the Sscrofa11.1 GTF: STAR silently dropped every pig
# exon and the matrices had no pig genes.  The checks at the end catch this.
#
# Single-species indexes (for the per-species re-alignment) are separate:
#   human  /w5home/bmoore/genomes/human_genome_38/star_solo_wtdTomato/
#   pig    built via run_aligner.sh's genome-build option.
#
# RAM: combined genome ~5.7 Gb; with --genomeSAsparseD 2 genomeGenerate needs
#      ~40-60 GB (limit set to 90 GB).  Disk: ~60 GB for FASTA/GTF + index.
# Time: a few hours at 24 threads.  Needs STAR on PATH (conda activate star_env).
#
# usage:  bash 01_build_combined_reference.sh          (refuses to overwrite)
#         FORCE=1 bash 01_build_combined_reference.sh  (rebuild in place)
# ============================================================================
#SBATCH --job-name=star_ref
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=12:00:00
set -euo pipefail

# ---- inputs ----------------------------------------------------------------
HUMAN_DIR=/w5home/bmoore/genomes/human_genome_38
PIG_DIR=/w5home/bmoore/genomes/pig_genome/ensemble
HUMAN_FA=$HUMAN_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly_tdTomato_NanoLuc.fa
HUMAN_GTF=$HUMAN_DIR/Homo_sapiens.GRCh38.108.filtered.gtf
PIG_FA=$PIG_DIR/Sus_scrofa.Sscrofa11.1.dna.primary_assembly.allchr.fa
PIG_GTF=$PIG_DIR/Sus_scrofa.Sscrofa11.1.109.chr.filtered.3pext.gtf
PREFIX_GTF=/w5home/bmoore/scRNAseq_library/prefix_gtf.py

# ---- outputs ---------------------------------------------------------------
OUT=/w5home/bmoore/genomes/human-pig_chimera
INDEX=star_solo
THREADS=24
HPFX=GRCh38        # species prefixes: every contig, gene_id and gene_name
PPFX=Sscrofa11     # carries one, so the species is readable everywhere

for f in "$HUMAN_FA" "$HUMAN_GTF" "$PIG_FA" "$PIG_GTF" "$PREFIX_GTF"; do
  [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done
if [ -e "$OUT/$INDEX/SA" ] && [ "${FORCE:-0}" != 1 ]; then
  echo "ERROR: $OUT/$INDEX already exists; set FORCE=1 to rebuild" >&2; exit 1
fi
mkdir -p "$OUT"/{human,pig}
cd "$OUT"

# ---- 1. prefix FASTA headers -----------------------------------------------
# `print $1` keeps only the sequence name: STAR uses just the first token, and
# keeping the description would leave the unprefixed name in the header.
prefix_fa() { awk -v p="$2" '/^>/{sub(/^>/,">"p"_"); print $1; next}{print}' "$1"; }
prefix_fa "$HUMAN_FA" "$HPFX" > human/${HPFX}.fa
prefix_fa "$PIG_FA"   "$PPFX" > pig/${PPFX}.fa

# ---- 2. prefix GTFs (seqnames + gene/transcript/exon IDs and names) --------
# No biotype filtering here: GENCODE/Ensembl-human and Ensembl-pig use
# different biotype vocabularies, so a `protein_coding` filter would remove a
# different fraction of each species and bias the species ratio.
python3 "$PREFIX_GTF" "$HUMAN_GTF" "$HPFX" human/${HPFX}.108.gtf
python3 "$PREFIX_GTF" "$PIG_GTF"   "$PPFX" pig/${PPFX}_v1.109_3pext.gtf

# ---- 3. concatenate --------------------------------------------------------
cat human/${HPFX}.fa pig/${PPFX}.fa > combined.fa
cat human/${HPFX}.108.gtf pig/${PPFX}_v1.109_3pext.gtf > combined.gtf

if grep '^>' combined.fa | sort | uniq -d | grep -q .; then
  echo "ERROR: duplicate contig names in combined.fa" >&2; exit 1
fi
# every GTF seqname must exist in the FASTA (the Berkshire/Sscrofa11.1 trap)
missing=$(comm -13 <(grep '^>' combined.fa | sed 's/^>//' | sort -u) \
                   <(grep -v '^#' combined.gtf | cut -f1 | sort -u) | head -5)
if [ -n "$missing" ]; then
  echo "ERROR: GTF seqnames absent from combined.fa, e.g.:" >&2
  echo "$missing" >&2; exit 1
fi

# ---- 4. combined STAR index ------------------------------------------------
# sjdbOverhang 149 = 150 bp reads - 1; SAsparseD 2 roughly halves index RAM
# at a modest mapping-speed cost.
STAR --runMode genomeGenerate \
     --runThreadN $THREADS \
     --genomeDir $INDEX \
     --genomeFastaFiles combined.fa \
     --sjdbGTFfile combined.gtf \
     --sjdbOverhang 149 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 90000000000 \
     --outFileNamePrefix _gg_

# ---- 5. checks -------------------------------------------------------------
nf=$(grep -c 'not found in Genome fasta' $INDEX/Log.out || true)
if [ "$nf" -gt 0 ]; then
  echo "ERROR: STAR dropped $nf GTF lines whose chromosome is not in the FASTA" >&2
  exit 1
fi
echo "genes in index by species (expect ~39,548 $HPFX and ~31,781 $PPFX):"
tail -n +2 $INDEX/geneInfo.tab | cut -f1 | sed 's/_.*//' | sort | uniq -c
grep -q "${HPFX}_NLuc_P2A_tdTomato" $INDEX/geneInfo.tab \
  && echo "reporter genes present" \
  || echo "WARNING: reporter gene ${HPFX}_NLuc_P2A_tdTomato not in index" >&2
