# export path
export PATH=/w5home/bmoore/scRNAseq/opt/cellranger-7.1.0:$PATH
# filter gtf
cellranger mkgtf Homo_sapiens.GRCh38.108.gtf Homo_sapiens.GRCh38.108.filtered.gtf --attribute=gene_biotype:protein_coding \
>                    --attribute=gene_biotype:lncRNA \
>                    --attribute=gene_biotype:antisense \
>                    --attribute=gene_biotype:IG_LV_gene \
>                    --attribute=gene_biotype:IG_V_gene \
>                    --attribute=gene_biotype:IG_V_pseudogene \
>                    --attribute=gene_biotype:IG_D_gene \
>                    --attribute=gene_biotype:IG_J_gene \
>                    --attribute=gene_biotype:IG_J_pseudogene \
>                    --attribute=gene_biotype:IG_C_gene \
>                    --attribute=gene_biotype:IG_C_pseudogene \
>                    --attribute=gene_biotype:TR_V_gene \
>                    --attribute=gene_biotype:TR_V_pseudogene \
>                    --attribute=gene_biotype:TR_D_gene \
>                    --attribute=gene_biotype:TR_J_gene \
>                    --attribute=gene_biotype:TR_J_pseudogene \
>                    --attribute=gene_biotype:TR_C_gene
# activate star env
conda activate star_env
# generate genome index with star using filtered gtf file from cellranger
 STAR --runMode genomeGenerate --runThreadN 8 --genomeDir star_solo_index/ \
 --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile \
 Homo_sapiens.GRCh38.108.filtered.gtf --genomeSAsparseD 3 --limitGenomeGenerateRAM 101410123018
 # get fastq files
 # sra does not have barcodes
 # download bam file:
 wget https://sra-pub-src-2.s3.amazonaws.com/SRR10759479/ORG_205_possorted_genome_bam.bam.1
 # convert to fastq
  ~/scRNAseq/opt/cellranger-7.1.0/lib/bin/bamtofastq ORG_205_possorted_genome_bam.bam.1 bam2fastq_out/
 # run star solo 
 # --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read
 # for reh data, R1 is barcodes, R2 is cDNA reads
 nohup STAR --genomeDir ~/human_genome_38/star_solo_index/ \
 --readFilesIn bamtofastq_S1_L001_R2_001.fastq.gz,bamtofastq_S1_L002_R2_001.fastq.gz,bamtofastq_S1_L003_R2_001.fastq.gz,bamtofastq_S1_L004_R2_001.fastq.gz \
 bamtofastq_S1_L001_R1_001.fastq.gz,bamtofastq_S1_L002_R1_001.fastq.gz,bamtofastq_S1_L003_R1_001.fastq.gz,bamtofastq_S1_L004_R1_001.fastq.gz \
 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/3M-february-2018.txt \
 --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --readFilesCommand zcat \
 --soloUMIdedup 1MM_CR --outFileNamePrefix ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_reh/ --runThreadN 8 &
 # lane 1
  nohup STAR --genomeDir ~/human_genome_38/star_solo_index/ \
 --readFilesIn bamtofastq_S1_L001_R2_001.fastq.gz bamtofastq_S1_L001_R1_001.fastq.gz \
 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/3M-february-2018.txt \
 --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --readFilesCommand zcat \
 --soloUMIdedup 1MM_CR --outFileNamePrefix ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_reh_lane1/ --runThreadN 8 > nohup.out1 &
 # lane 2
   nohup STAR --genomeDir ~/human_genome_38/star_solo_index/ \
 --readFilesIn bamtofastq_S1_L002_R2_001.fastq.gz bamtofastq_S1_L002_R1_001.fastq.gz \
 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/3M-february-2018.txt \
 --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --readFilesCommand zcat \
 --soloUMIdedup 1MM_CR --outFileNamePrefix ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_reh_lane2/ --runThreadN 8 > nohup.out2 &
 # lane 3
    nohup STAR --genomeDir ~/human_genome_38/star_solo_index/ \
 --readFilesIn bamtofastq_S1_L003_R2_001.fastq.gz bamtofastq_S1_L003_R1_001.fastq.gz \
 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/3M-february-2018.txt \
 --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --readFilesCommand zcat \
 --soloUMIdedup 1MM_CR --outFileNamePrefix ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_reh_lane3/ --runThreadN 8 > nohup.out3 &
 # lane 4
 nohup STAR --genomeDir ~/human_genome_38/star_solo_index/ \
 --readFilesIn bamtofastq_S1_L004_R2_001.fastq.gz bamtofastq_S1_L004_R1_001.fastq.gz \
 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/3M-february-2018.txt \
 --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --readFilesCommand zcat \
 --soloUMIdedup 1MM_CR --outFileNamePrefix ~/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_reh_lane4/ --runThreadN 8 > nohup.out4 &
