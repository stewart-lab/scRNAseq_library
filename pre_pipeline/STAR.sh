#!/bin/bash

# Path to the configuration file
CONFIG_FILE="/config.json"

declare -A GENOME_DIRS
GENOME_DIRS["human"]="/usr/local/bin/human_genome/star_solo_index"
GENOME_DIRS["pig"]="/usr/local/bin/pig_genome/Sus_scrofa_genome_forstar_MT"

# Define SOLO_FEATURES and SOLO_CELL_FILTER
SOLO_FEATURES="Gene GeneFull SJ Velocyto"
SOLO_CELL_FILTER="EmptyDrops_CR"

# Fetch all sample names from the config file
NUMBER_OF_SAMPLES=$(jq -r '.Number_of_samples' $CONFIG_FILE)

for ((SAMPLE_IDX=1; SAMPLE_IDX<=NUMBER_OF_SAMPLES; SAMPLE_IDX++)); do
  # Loop over each sample based on NUMBER_OF_SAMPLES
  SAMPLE_KEY="SAMPLE_$SAMPLE_IDX"  # Construct the key to access sample configurations

  # Fetch the species and determine the genome directory
  SPECIES=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.species" $CONFIG_FILE)
  GENOME_DIR=${GENOME_DIRS[$SPECIES]}

  # Check if genome directory is set, else skip the sample with a warning
  if [ -z "$GENOME_DIR" ]; then
    echo "Warning: Genome directory not set for species '$SPECIES' in sample $SAMPLE_KEY. Skipping."
    continue
  fi
  # Fetch configuration for each sample using the constructed key
  SAMPLE_NAME=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.NAME" $CONFIG_FILE)  # Unique name for the sample
  NUM_LANES=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.NUM_LANES" $CONFIG_FILE)
  CHEMISTRY_VERSION=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.CHEMISTRY_VERSION" $CONFIG_FILE)
  SOLO_TYPE=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.SOLO_TYPE" $CONFIG_FILE)
  SOLO_MULTI_MAPPERS=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.SOLO_MULTI_MAPPERS" $CONFIG_FILE)
  SOLO_UMI_DEDUP=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.SOLO_UMI_DEDUP" $CONFIG_FILE)
  RUN_THREAD_N=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.RUN_THREAD_N" $CONFIG_FILE)
  READ_FILES_COMMAND=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.READ_FILES_COMMAND" $CONFIG_FILE)
  clip5pNbases=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.clip5pNbases" $CONFIG_FILE)
  isBarcodeFollowedbyReads=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.isBarcodeFollowedbyReads" $CONFIG_FILE)

  # Determine SOLO_CB_WHITELIST and SOLO_UMI_LEN based on CHEMISTRY_VERSION
  if [ "$CHEMISTRY_VERSION" == "V2" ]; then
    SOLO_CB_WHITELIST="/usr/local/bin/737K-august-2016.txt"
    SOLO_UMI_LEN=10
  elif [ "$CHEMISTRY_VERSION" == "V3" ]; then
    SOLO_CB_WHITELIST="/usr/local/bin/3M-february-2018.txt"
    SOLO_UMI_LEN=12
  else
    echo "Invalid chemistry version for $SAMPLE_NAME ($SAMPLE_KEY). Exiting."
    continue
  fi

  # Process each lane for the current sample
  for i in $(seq 1 $NUM_LANES); do
    cDNA_LANE=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.cDNA_LANE${i}" $CONFIG_FILE)
    BARCODE_LANE=$(jq -r ".fastq_alignment.${SAMPLE_KEY}.BARCODE_LANE${i}" $CONFIG_FILE)

    # Define STAR options based on the configuration
    STAR_OPTIONS=""
    if [ "$isBarcodeFollowedbyReads" = true ]; then
      soloBarcodeMate=2
      soloCBstart=1
      soloCBlen=16
      soloUMIstart=17
      STAR_OPTIONS="--soloBarcodeMate $soloBarcodeMate --clip5pNbases $clip5pNbases --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart"
    fi

    # Run STAR command with options for the current sample and lane
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
    --runThreadN $RUN_THREAD_N > "/shared_volume/nohup_${SAMPLE_NAME}_lane${i}.out" &
  done
done



