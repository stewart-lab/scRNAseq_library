#!/bin/bash

# Path to the configuration file
CONFIG_FILE="/config.json"

umask 000

# Ensure the shared volume is writable
if [ ! -w /shared_mount ]; then
    echo "Cannot write to /shared_volume. Please check permissions."
    exit 1
fi

# Define SOLO_FEATURES and SOLO_CELL_FILTER
SOLO_FEATURES="Gene GeneFull SJ Velocyto"
SOLO_CELL_FILTER="EmptyDrops_CR"

# define genome dir
GENOME_DIR="/genome_dir/"

# Python function to extract JSON values
get_json_value() {
    local key="$1"
    python3 -c "import json, sys; data = json.load(open('$CONFIG_FILE')); print(data.get('$key', ''))"
}

get_json_nested_value() {
    local path="$1"
    python3 -c "import json, sys; data = json.load(open('$CONFIG_FILE')); keys = '$path'.split('.'); value = data; [value := value.get(k, '') for k in keys]; print(value)"
}

# Fetch all sample names from the config file
NUMBER_OF_SAMPLES=$(get_json_value "Number_of_samples")

for ((SAMPLE_IDX=1; SAMPLE_IDX<=NUMBER_OF_SAMPLES; SAMPLE_IDX++)); do
  # Loop over each sample based on NUMBER_OF_SAMPLES
  SAMPLE_KEY="SAMPLE_$SAMPLE_IDX"  # Construct the key to access sample configurations

  # Fetch the species
  SPECIES=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.species")
  
    # Fetch configuration for each sample using the constructed key
  SAMPLE_NAME=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.NAME")  # Unique name for the sample
  NUM_LANES=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.NUM_LANES")
  CHEMISTRY_VERSION=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.CHEMISTRY_VERSION")
  SOLO_TYPE=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.SOLO_TYPE")
  SOLO_MULTI_MAPPERS=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.SOLO_MULTI_MAPPERS")
  SOLO_UMI_DEDUP=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.SOLO_UMI_DEDUP")
  RUN_THREAD_N=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.RUN_THREAD_N")
  READ_FILES_COMMAND=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.READ_FILES_COMMAND")
  clip5pNbases=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.clip5pNbases")
  isBarcodeFollowedbyReads=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.isBarcodeFollowedbyReads")
  soloStrand=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.soloStrand")

  # Determine SOLO_CB_WHITELIST and SOLO_UMI_LEN based on CHEMISTRY_VERSION
  if [ "$CHEMISTRY_VERSION" == "V2" ]; then
    SOLO_CB_WHITELIST="data/737K-august-2016.txt"
    SOLO_UMI_LEN=10
  elif [ "$CHEMISTRY_VERSION" == "V3" ]; then
    SOLO_CB_WHITELIST="data/3M-february-2018.txt"
    SOLO_UMI_LEN=12
  elif [ "$CHEMISTRY_VERSION" == "V4" ]; then
    if [ "$FIVE_PRIME" = true ]; then
      SOLO_CB_WHITELIST="data/3M-5pgex-jan-2023.txt"
      SOLO_UMI_LEN=12
    else
      SOLO_CB_WHITELIST="data/3M-3pgex-may-2023_TRU.txt"
      SOLO_UMI_LEN=12
    fi
  else
    echo "Invalid chemistry version for $SAMPLE_NAME ($SAMPLE_KEY). Exiting."
    continue
  fi

  # Process each lane for the current sample
  for i in $(seq 1 $NUM_LANES); do
    cDNA_LANE=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.cDNA_LANE${i}")
    BARCODE_LANE=$(get_json_nested_value "fastq_alignment.${SAMPLE_KEY}.BARCODE_LANE${i}")

    # Define STAR options based on the configuration
    STAR_OPTIONS=""

    if [ "$isBarcodeFollowedbyReads" = true ]; then
      soloBarcodeMate=2
      soloCBstart=1
      soloCBlen=16
      soloUMIstart=17
      STAR_OPTIONS+=" --soloBarcodeMate $soloBarcodeMate --clip5pNbases $clip5pNbases --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart"
    elif [ -n "$clip5pNbases" ]; then
      soloCBstart=1
      soloCBlen=16
      soloUMIstart=17
      STAR_OPTIONS+=" --clip5pNbases $clip5pNbases --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart --soloBarcodeReadLength 0"
    fi

    if [ -n "$soloStrand" ]; then
      STAR_OPTIONS+=" --soloStrand $soloStrand"
    fi
    # Create output directory in shared volume
    mkdir -p /shared_mount/${SAMPLE_NAME}_lane${i}

    # Run STAR command with options for the current sample and lane
    STAR --genomeDir $GENOME_DIR \
    --readFilesIn "$cDNA_LANE" "$BARCODE_LANE" \
    --readFilesCommand $READ_FILES_COMMAND \
    --soloUMIlen $SOLO_UMI_LEN --soloType $SOLO_TYPE \
    --soloCBwhitelist $SOLO_CB_WHITELIST \
    --soloFeatures $SOLO_FEATURES \
    --soloCellFilter $SOLO_CELL_FILTER --soloMultiMappers $SOLO_MULTI_MAPPERS \
    --soloUMIdedup $SOLO_UMI_DEDUP \
    $STAR_OPTIONS \
    --outFileNamePrefix "/shared_mount/${SAMPLE_NAME}_lane${i}/" \
    --runThreadN $RUN_THREAD_N \
    --runDirPerm All_RWX # Add this line to set file permissions

    # Log the completion of each lane
    echo "Completed processing ${SAMPLE_NAME}_lane${i}" >> /shared_mount/alignment_log.txt
  done
done

echo "All samples processed successfully" >> /shared_mount/alignment_log.txt
