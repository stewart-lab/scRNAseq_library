# star solo
conda activate star_env

# generate genome index with star using filtered gtf file from cellranger
 STAR --runMode genomeGenerate --runThreadN 8 --genomeDir Sus_scrofa_genome_forstar/ \
 --genomeFastaFiles Sus_scrofa_berkshire.Berkshire_pig_v1.dna.toplevel.fa --sjdbGTFfile \
 Sus_scrofa_berkshire.Berkshire_pig_v1.109.filtered.gtf --genomeSAsparseD 3 --limitGenomeGenerateRAM 101410123018
 
# run star-solo 
#/path/to/STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  [...other parameters...] --soloType ... --soloCBwhitelist ...
STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome/star/ --readFilesIn DATA/Mini_Fastqs/1_S1_L001_R2_001_mini.fastq.gz \
DATA/Mini_Fastqs/1_S1_L001_R1_001_mini.fastq.gz --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt
## this command didn't work because uses the index made from cell ranger. given the index
## generated above, should work
## Also note the white list and --soloUMIlen flag are for 10x Chromium v3.

# rerun star-solo
STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Mini_Fastqs/1_S1_L001_R2_001_mini.fastq.gz \
DATA/Mini_Fastqs/1_S1_L001_R1_001_mini.fastq.gz --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt
# add in decompression commmand:
--readFilesCommand zcat

# more complicated command with lanes
STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Mini_Fastqs/1_S1_L001_R2_001_mini.fastq.gz,DATA/Mini_Fastqs/1_S1_L002_R2_001_mini.fastq.gz DATA/Mini_Fastqs/1_S1_L001_R1_001_mini.fastq.gz,DATA/Mini_Fastqs/1_S1_L002_R1_001_mini.fastq.gz --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --readFilesCommand zcat --outFileNamePrefix output2/

# empty drops algorithm
STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Mini_Fastqs/1_S1_L001_R2_001_mini.fastq.gz,DATA/Mini_Fastqs/1_S1_L002_R2_001_mini.fastq.gz DATA/Mini_Fastqs/1_S1_L001_R1_001_mini.fastq.gz,DATA/Mini_Fastqs/1_S1_L002_R1_001_mini.fastq.gz --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --readFilesCommand zcat --outFileNamePrefix output3/

# sample 1 both lanes
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/1_S1_L001_R2_001.fastq.gz,DATA/Orig_Fastqs/1_S1_L002_R2_001.fastq.gz DATA/Orig_Fastqs/1_S1_L001_R1_001.fastq.gz,DATA/Orig_Fastqs/1_S1_L002_R1_001.fastq.gz --soloUMIlen 12 --soloType 
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --readFilesCommand zcat --
outFileNamePrefix output_S1/ --runThreadN 8 &

# sample 2 both lanes
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/2_S2_L001_R2_001.fastq.gz,DATA/Orig_Fastqs/2_S2_L002_R2_001.fastq.gz DATA/Orig_Fastqs/2_S2_L001_R1_001.fastq.gz,DATA/Orig_Fastqs/2_S2_L002_R1_001.fastq.gz --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --readFilesCommand zcat --outFileNamePrefix output_S2/ --runThreadN 8 &

# add multimapping EM, UMI demultiplexing
--soloMultiMappers EM 
--soloUMIdedup 1MM_CR
# sample 1 both lanes
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/1_S1_L001_R2_001.fastq.gz,DATA/Orig_Fastqs/1_S1_L002_R2_001.fastq.gz DATA/Orig_Fastqs/1_S1_L001_R1_001.fastq.gz,DATA/Orig_Fastqs/1_S1_L002_R1_001.fastq.gz --soloUMIlen 12 --soloType \
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloUMIdedup 1MM_CR --readFilesCommand zcat --outFileNamePrefix output_S1_mm/ --runThreadN 8 &
# lane 1
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/1_S1_L001_R2_001.fastq.gz DATA/Orig_Fastqs/1_S1_L001_R1_001.fastq.gz --soloUMIlen 12 --soloType \
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloUMIdedup 1MM_CR --readFilesCommand zcat --outFileNamePrefix output_S1-1_mm/ --runThreadN 8 &
# lane 2
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/1_S1_L002_R2_001.fastq.gz DATA/Orig_Fastqs/1_S1_L002_R1_001.fastq.gz --soloUMIlen 12 --soloType \
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloUMIdedup 1MM_CR --readFilesCommand zcat --outFileNamePrefix output_S1-2_mm/ --runThreadN 8 &
# sample 2 both lanes- #might need to redo this one
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/2_S2_L001_R2_001.fastq.gz,DATA/Orig_Fastqs/2_S2_L002_R2_001.fastq.gz DATA/Orig_Fastqs/2_S2_L001_R1_001.fastq.gz,DATA/Orig_Fastqs/2_S2_L002_R1_001.fastq.gz --soloUMIlen 12 --soloType \
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloUMIdedup 1MM_CR --readFilesCommand zcat --outFileNamePrefix output_S2_mm/ --runThreadN 8 &
# lane 1
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/2_S2_L001_R2_001.fastq.gz DATA/Orig_Fastqs/2_S2_L001_R1_001.fastq.gz --soloUMIlen 12 --soloType \
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloUMIdedup 1MM_CR --readFilesCommand zcat --outFileNamePrefix output_S2-1_mm/ --runThreadN 8 &
# lane 2
nohup STAR --genomeDir ~/pig_genome/ensemble/Sus_scrofa_genome_forstar/ --readFilesIn DATA/Orig_Fastqs/2_S2_L002_R2_001.fastq.gz DATA/Orig_Fastqs/2_S2_L002_R1_001.fastq.gz --soloUMIlen 12 --soloType \
CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloUMIdedup 1MM_CR --readFilesCommand zcat --outFileNamePrefix output_S2-2_mm/ --runThreadN 8 &