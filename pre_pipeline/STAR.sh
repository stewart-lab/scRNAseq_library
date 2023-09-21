#!/bin/bash

# Import variables from config.json
CONFIG_FILE="/config.json"

NUM_LANES=$(jq -r '.fastq_alignment.NUM_LANES' $CONFIG_FILE)
GENOME_DIR="/usr/local/bin/star_solo_index/"
SAMPLE_NAME=$(jq -r '.fastq_alignment.SAMPLE_NAME' $CONFIG_FILE)
CHEMISTRY_VERSION=$(jq -r '.fastq_alignment.CHEMISTRY_VERSION' $CONFIG_FILE)
SOLO_TYPE=$(jq -r '.fastq_alignment.SOLO_TYPE' $CONFIG_FILE)
SOLO_FEATURES=$(jq -r '.fastq_alignment.SOLO_FEATURES' $CONFIG_FILE)
SOLO_CELL_FILTER=$(jq -r '.fastq_alignment.SOLO_CELL_FILTER' $CONFIG_FILE)
SOLO_MULTI_MAPPERS=$(jq -r '.fastq_alignment.SOLO_MULTI_MAPPERS' $CONFIG_FILE)
READ_FILES_COMMAND=$(jq -r '.fastq_alignment.READ_FILES_COMMAND' $CONFIG_FILE)
SOLO_UMI_DEDUP=$(jq -r '.fastq_alignment.SOLO_UMI_DEDUP' $CONFIG_FILE)
RUN_THREAD_N=$(jq -r '.fastq_alignment.RUN_THREAD_N' $CONFIG_FILE)

# Define read files and solo options based on chemistry version
if [ "$CHEMISTRY_VERSION" == "V2" ]; then
  SOLO_CB_WHITELIST="/usr/local/bin/737K-august-2016.txt"  # V2 whitelist
  SOLO_UMI_LEN=10  # UMI length for V2
elif [ "$CHEMISTRY_VERSION" == "V3" ]; then
  SOLO_CB_WHITELIST="/usr/local/bin/3M-february-2018.txt"  # V3 whitelist
  SOLO_UMI_LEN=12  # UMI length for V3
else
  echo "Invalid chemistry version. Exiting."
  exit 1
fi

# Loop through the number of lanes to run STAR command
for i in $(seq 1 $NUM_LANES); do
  cDNA_LANE=$(jq -r ".fastq_alignment.cDNA_LANE${i}" $CONFIG_FILE)
  BARCODE_LANE=$(jq -r ".fastq_alignment.BARCODE_LANE${i}" $CONFIG_FILE)
  
  nohup STAR --genomeDir $GENOME_DIR \
  --readFilesIn "$cDNA_LANE" "$BARCODE_LANE" \
  --soloUMIlen $SOLO_UMI_LEN --soloType $SOLO_TYPE \
  --soloCBwhitelist $SOLO_CB_WHITELIST \
  --soloFeatures $SOLO_FEATURES \
  --soloCellFilter $SOLO_CELL_FILTER --soloMultiMappers $SOLO_MULTI_MAPPERS \
  --readFilesCommand $READ_FILES_COMMAND --soloUMIdedup $SOLO_UMI_DEDUP \
  --outFileNamePrefix "/shared_volume/${SAMPLE_NAME}_lane${i}/" \
  --runThreadN $RUN_THREAD_N > "/shared_volume/nohup_lane${i}.out" &
done


