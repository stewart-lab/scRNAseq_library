#!/bin/bash

# Import variables from config.json
CONFIG_FILE="/config.json"

NUM_LANES=$(jq -r '.fastq_alignment.NUM_LANES' $CONFIG_FILE)
GENOME_DIR="/usr/local/bin/star_solo_index/"
SAMPLE_NAME=$(jq -r '.fastq_alignment.SAMPLE_NAME' $CONFIG_FILE)
CHEMISTRY_VERSION=$(jq -r '.fastq_alignment.CHEMISTRY_VERSION' $CONFIG_FILE)
SOLO_TYPE=$(jq -r '.fastq_alignment.SOLO_TYPE' $CONFIG_FILE)
SOLO_MULTI_MAPPERS=$(jq -r '.fastq_alignment.SOLO_MULTI_MAPPERS' $CONFIG_FILE)
SOLO_UMI_DEDUP=$(jq -r '.fastq_alignment.SOLO_UMI_DEDUP' $CONFIG_FILE)
RUN_THREAD_N=$(jq -r '.fastq_alignment.RUN_THREAD_N' $CONFIG_FILE)
soloBarcodeMate=$(jq -r '.fastq_alignment.soloBarcodeMate' $CONFIG_FILE)
clip5pNbases=$(jq -r '.fastq_alignment.clip5pNbases' $CONFIG_FILE)
isBarcodeFollowedbyReads=$(jq -r '.fastq_alignment.isBarcodeFollowedbyReads' $CONFIG_FILE)
READ_FILES_COMMAND=$(jq -r '.fastq_alignment.READ_FILES_COMMAND' $CONFIG_FILE)
SOLO_FEATURES="Gene GeneFull SJ Velocyto"
SOLO_CELL_FILTER="EmptyDrops_CR"


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



for i in $(seq 1 $NUM_LANES); do
  cDNA_LANE=$(jq -r ".fastq_alignment.cDNA_LANE${i}" $CONFIG_FILE)
  BARCODE_LANE=$(jq -r ".fastq_alignment.BARCODE_LANE${i}" $CONFIG_FILE)

  STAR_OPTIONS=""
  
  if [ "$isBarcodeFollowedbyReads" = true ]; then
    soloBarcodeMate=2
    soloCBstart=1
    soloCBlen=16
    soloUMIstart=17
    STAR_OPTIONS="--soloBarcodeMate $soloBarcodeMate --clip5pNbases $clip5pNbases --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart"
  fi

  nohup STAR --genomeDir $GENOME_DIR \
  --readFilesIn "$cDNA_LANE" "$BARCODE_LANE" \
  --readFilesCommand $READ_FILES_COMMAND \
  --soloUMIlen $SOLO_UMI_LEN --soloType $SOLO_TYPE \
  --soloCBwhitelist $SOLO_CB_WHITELIST \
  --soloFeatures $SOLO_FEATURES \
  --soloCellFilter $SOLO_CELL_FILTER --soloMultiMappers $SOLO_MULTI_MAPPERS \
  --soloUMIdedup $SOLO_UMI_DEDUP \
  $STAR_OPTIONS \
  --outFileNamePrefix "/shared_volume/${SAMPLE_NAME}_lane${i}/" \
  --runThreadN $RUN_THREAD_N > "/shared_volume/nohup_lane${i}.out" &
done


