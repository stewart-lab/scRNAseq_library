#!/usr/bin/env bash
# ============================================================================
# 01 -- download GRCh38 + Sscrofa11.1, build one combined STAR index plus the
#       two single-species indexes used for re-quantification in step 04.
#
# RAM: the combined genome is ~5.6 Gb.  genomeGenerate wants ~60 GB with
#      default settings; --genomeSAsparseD 2 brings it to ~40 GB at a modest
#      mapping-speed cost (3 -> ~32 GB, slower still).  This step does not fit
#      on a 16 GB machine -- run it on a cluster node or a large cloud VM.
# Disk: ~250 GB for downloads + indexes.
# Time: 3-6 h at 16 threads.
# ============================================================================
#SBATCH --job-name=star_ref
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00


mkdir -p human-pig_chimera/{human,pig}
cd human-pig_chimera/

#GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}
#ENSEMBL=https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}
# already downloaded

# ---- 1. prefix contig names ----------------------------------------------
# GENCODE uses chr1.., Ensembl uses 1.. -- after prefixing both are unique, but
# prefix anyway so the species is readable off any contig/gene/feature name.
# export prefix
export HPFX=GRCh38
export PPFX=Sscrofa11

cat Homo_sapiens.GRCh38.dna_sm.primary_assembly_tdTomato_NanoLuc.fa \
| awk -v p="$HPFX" '/^>/{sub(/^>/,">"p"_"); print $1; next}{print}' \
> "${HPFX}.fa"
mv human_genome_38/GRCh38.fa human-pig_chimera/human/

cat Sus_scrofa_berkshire.Berkshire_pig_v1.dna.toplevel.fa \
| awk -v p="$PPFX" '/^>/{sub(/^>/,">"p"_"); print $1; next}{print}' \
> "${PPFX}.fa"
mv Sscrofa11.fa  ../../human-pig_chimera/pig/
# NB: `print $1` drops the FASTA description field -- STAR only uses the first
# token anyway, and keeping it would leave the unprefixed name in the header.

# ---- 3. prefix GTFs -------------------------------------------------------
python /w5home/bmoore/scRNAseq_library/prefix_gtf.py \
Homo_sapiens.GRCh38.108.filtered.gtf \
"$HPFX" "/w5home/bmoore/genomes/human-pig_chimera/${HPFX}.gtf"
mv GRCh38.gtf human/GRCh38.108.gtf

python /w5home/bmoore/scRNAseq_library/prefix_gtf.py \
Sus_scrofa_berkshire.Berkshire_pig_v1.109.filtered.gtf \
"$PPFX" "/w5home/bmoore/genomes/human-pig_chimera/${PPFX}.gtf"
mv Sscrofa11.gtf pig/Sscrofa11_v1.109.gtf

# No biotype filtering.  GENCODE and Ensembl-pig use different biotype
# vocabularies, so a naive `protein_coding` filter removes a different fraction
# from each species and biases the species ratio you are about to measure.

# ---- 4. concatenate -------------------------------------------------------
cat human/GRCh38.fa pig/Sscrofa11.fa > combined.fa
cat human/GRCh38.108.gtf pig/Sscrofa11_v1.109.gtf > combined.gtf

# sanity: no duplicate contig names
grep '^>' "combined.fa" | sort | uniq -d | grep . && \
  { echo "ERROR: duplicate contig names" >&2; exit 1; } || true

# ---- 5. combined index (used for species assignment only) -----------------
STAR --runMode genomeGenerate \
     --runThreadN 24 \
     --genomeDir star_solo \
     --genomeFastaFiles combined.fa \
     --sjdbGTFfile combined.gtf \
     --sjdbOverhang 149 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 90000000000 \
     --outFileNamePrefix _gg_


