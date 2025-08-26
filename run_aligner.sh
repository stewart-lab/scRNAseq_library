#!/bin/bash

echo "Step 1: Importing DATA_DIR from config.json"
CONFIG_FILE="./config.json"
DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

echo "Step 2: Confirming data alignment"
read -p "Would you like to run alignment? [y/N]: " confirm

if [[ "$confirm" =~ ^[Yy]$ ]]; then
  echo "Step 2.1: Removing and recreating the SHARED_MOUNT"
  
  SHARED_MOUNT="./shared_mount"
  rm -rf "$SHARED_MOUNT"
  mkdir -p "$SHARED_MOUNT"
  chmod 777 "$SHARED_MOUNT"

  echo "Step 2.2: Building Docker image for alignment"
  docker build -t scaligner_v2_with_genomes_and_jq ./pre_pipeline

  echo "Step 2.3: Running Docker container for alignment"
  docker run -d \
    -v "$(realpath "$DATA_DIR"):/data:ro" \
    -v "$(realpath "$SHARED_MOUNT"):/shared_mount" \
    -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
    scaligner_v2_with_genomes_and_jq /bin/bash -c "conda run -n star_env /bin/bash -c 'cd /src && ./STAR.sh'"
fi

echo "Step 3: Updating config.json"
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json