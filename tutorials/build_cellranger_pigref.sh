# build sc genome index from scratch with CellRanger

# download CellRanger
# check: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
# can also install via conda: conda install -c hcc cellranger

 wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1685069631&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODUwNjk2MzF9fX1dfQ__&Signature=JzoThaT8sPx~KK0EKr1EEzir21yg7hBuxUUq99DJ5WlZ2OeNkiIq5jkKF2TeBnoJdwbmNH047COtOQHIT2oqVAo40DR3b6HDhYUiwbG4ho4vl2a3y-T0syLz4PUZpwlgAn7C0PeqxifibE3wytWrAHKGfkfzrGQjm6EozcloXFlN4c9UQnmL3VgFgj~V62A3wSVGX6fGYWrga-hsry~v4sCI3R3t2cAMf9M9p1lJbc6DgyyPPM1cburCcwdUQv~LAqVSNl6UppZHXuK1P-HPQYx92ajB~IbYz-3Z6d1KDkbBYLU3h56YAg5EZe3hfFgmqkHv7Ami-e0qZbDwBNGAHg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
--2023-05-25 09:54:46--  https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1685069631&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODUwNjk2MzF9fX1dfQ__&Signature=JzoThaT8sPx~KK0EKr1EEzir21yg7hBuxUUq99DJ5WlZ2OeNkiIq5jkKF2TeBnoJdwbmNH047COtOQHIT2oqVAo40DR3b6HDhYUiwbG4ho4vl2a3y-T0syLz4PUZpwlgAn7C0PeqxifibE3wytWrAHKGfkfzrGQjm6EozcloXFlN4c9UQnmL3VgFgj~V62A3wSVGX6fGYWrga-hsry~v4sCI3R3t2cAMf9M9p1lJbc6DgyyPPM1cburCcwdUQv~LAqVSNl6UppZHXuK1P-HPQYx92ajB~IbYz-3Z6d1KDkbBYLU3h56YAg5EZe3hfFgmqkHv7Ami-e0qZbDwBNGAHg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA
# untar- self contained, so no need to install anything else
tar -xzvf cellranger-7.1.0.tar.gz

# download human reference
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# untar
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

# export path
export PATH=/w5home/bmoore/scRNAseq/opt/cellranger-7.1.0:$PATH

# test
cellranger testrun --id=tiny
# output: tiny/tiny.mri.tgz

# use cell ranger to build pig reference

# get pig reference from ensemble (also in ncbi but ensemble is recommended by cell ranger)
wget https://ftp.ensembl.org/pub/release-109/fasta/sus_scrofa_berkshire/dna/Sus_scrofa_berkshire.Berkshire_pig_v1.dna.toplevel.fa.gz
gunzip Sus_scrofa_berkshire.Berkshire_pig_v1.dna.toplevel.fa.gz
# get gtf from ensemble
wget https://ftp.ensembl.org/pub/release-109/gtf/sus_scrofa_berkshire/Sus_scrofa_berkshire.Berkshire_pig_v1.109.gtf.gz
gunzip Sus_scrofa_berkshire.Berkshire_pig_v1.109.gtf.gz
# filter gtf file
# cellranger mkgtf input.gtf output.gtf --attribute=key:allowable_value
cellranger mkgtf Sus_scrofa_berkshire.Berkshire_pig_v1.109.gtf Sus_scrofa_berkshire.Berkshire_pig_v1.109.filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lncRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene
# Index the FASTA and GTF files with mkref
# cellranger mkref --genome=output_genome --fasta=input.fa --genes=input.gtf
cellranger mkref --genome=Sus_scrofa_genome --fasta=Sus_scrofa_berkshire.Berkshire_pig_v1.dna.toplevel.fa \
--genes=Sus_scrofa_berkshire.Berkshire_pig_v1.109.filtered.gtf --ref-version v1 --nthreads 8