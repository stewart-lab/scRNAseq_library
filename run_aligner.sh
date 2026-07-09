#!/bin/bash

echo "Step 1: Importing DATA_DIR from config.json"
CONFIG_FILE="./config.json"
DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

GENOME_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['GENOME_DIR'])")
echo "GENOME_DIR imported as $GENOME_DIR"
GENOME_INDEX=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['GENOME_INDEX_DIR'])")
GENOME_INDEX_DIR=$GENOME_DIR$GENOME_INDEX
echo "GENOME_INDEX_DIR imported as $GENOME_INDEX_DIR"

echo "Step 2: Confirming data alignment"
read -p "Would you like to run alignment? [y/N]: " confirm

if [[ "$confirm" =~ ^[Yy]$ ]]; then
  echo "Step 2.1: Removing and recreating the SHARED_MOUNT"
  
  SHARED_MOUNT="./shared_mount"
  # rm -rf "$SHARED_MOUNT"
  # mkdir -p "$SHARED_MOUNT"
  chmod 777 "$SHARED_MOUNT"

  read -p "Are you running parse pipeline? [y/N]: " PARSE
  if [[ "$PARSE" =~ ^[Yy]$ ]]; then
    echo "Step 2.2: Building Parse Docker image for alignment"
    docker build -t parse-biosciences ./parse_pipeline
    
      read -p "Do you need to build a genome index? [y/N]: " confirm
      if [[ "$confirm" =~ ^[Yy]$ ]]; then
        echo "Building genome index"
        docker run -it \
          -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
          -v "$(realpath "$GENOME_DIR"):/genome_dir" \
          parse-biosciences /bin/bash -c "conda run -n base /bin/bash -c 'cd /src && ./parse_genome_build.sh'"
      fi
      echo "Step 2.3: Running Docker container for alignment"
      docker run -d \
          -v "$(realpath "$DATA_DIR"):/data:ro" \
          -v "$(realpath "$SHARED_MOUNT"):/shared_mount" \
          -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
          -v "$(realpath "$GENOME_INDEX_DIR"):/genome_index_dir" \
          parse-biosciences /bin/bash -c "conda run -n base /bin/bash -c 'cd /src && ./parse.sh'"
  else
    echo "Step 2.2: Building Docker image for alignment"
    docker build -t scaligner_v2_with_genomes_and_jq ./pre_pipeline

    echo "Step 2.3: Running Docker container for alignment"
    docker run -d \
      -v "$(realpath "$DATA_DIR"):/data:ro" \
      -v "$(realpath "$SHARED_MOUNT"):/shared_mount" \
      -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
      scaligner_v2_with_genomes_and_jq /bin/bash -c "conda run -n star_env /bin/bash -c 'cd /src && ./STAR.sh'"
  fi
fi

echo "Step 3: Updating config.json"
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json