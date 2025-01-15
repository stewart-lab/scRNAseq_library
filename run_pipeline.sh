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
  docker run -it \
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

echo "Step 4: Preparing for Docker container run"
cd output
output_dir=$(pwd)
cd ..
shared_mount_dir="$(pwd)/shared_mount"

echo "Step 5: Building Docker image for the main container"
docker build -t seuratv5 ./sc_pipeline

echo "Step 6: Running the main Docker container"
chmod 777 "$output_dir"
chmod 777 ./sc_pipeline/src/config.json

# Run pipeline with simplified commands
docker run -it \
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_mount_dir",target=/scRNA-seq/shared_mount \
  --mount type=bind,source="$(realpath ./sc_pipeline/src/config.json)",target=/scRNA-seq/src/config.json \
  --entrypoint=/bin/sh \
<<<<<<< HEAD
  seuratv5 -c "python3 /scRNA-seq/get_data.py && /opt/conda/envs/scrnaseq/bin/Rscript /scRNA-seq/script.R"
=======
  seuratv5 -c "python3 /scRNA-seq/get_data.py $DATA_FLAG" && \

# Second command for script.R with full path to Rscript
docker run -it \
  --memory="64g" \
  --memory-swap="64g" \
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_mount_dir",target=/scRNA-seq/shared_mount \
  --mount type=bind,source="$(realpath ./sc_pipeline/src/config.json)",target=/scRNA-seq/src/config.json \
  --entrypoint=/bin/sh \
  seuratv5 -c "/opt/conda/envs/scrnaseq/bin/Rscript /scRNA-seq/script.R"







>>>>>>> 02d3a70 (wip)
