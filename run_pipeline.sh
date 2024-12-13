#!/bin/bash

echo "Step 1: Importing DATA_DIR from config.json"
CONFIG_FILE="./config.json"
DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

echo "Step 2: Confirming data alignment"
read -p "Have you loaded new data or would you like to realign? [y/N]: " confirm

if [[ "$confirm" =~ ^[Yy]$ ]]; then
  echo "Step 2.1: Removing and recreating the SHARED_MOUNT"
  
  SHARED_MOUNT="./shared_mount"

  rm -rf "$SHARED_MOUNT"
  mkdir -p "$SHARED_MOUNT"
  chmod 777 "$SHARED_MOUNT"

  echo "Step 2.2: Building Docker image for alignment"
  # Build the Docker image from the pre_pipeline directory
  docker build -t scaligner_v2_with_genomes_and_jq ./pre_pipeline

  echo "Step 2.3: Running Docker container for alignment"
  docker run -it \
    -v "$(realpath "$DATA_DIR"):/data:ro" \
    -v "$(realpath "$SHARED_MOUNT"):/shared_mount" \
    -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
    scaligner_v2_with_genomes_and_jq /bin/bash -c "conda run -n star_env /bin/bash -c 'cd /src && ./STAR.sh'"

else
  echo "Skipped scaligner."
  
  echo "Step 2.4: Asking for data flag"
  # Ask the user which flag to use for get_data.py
  read -p "If you'd like to load a stored experiment, select 'data'. If you have aligned FASTQs loaded and changed pipeline parameters, select 'fastq' [data/fastq]: " data_flag

  if [ "$data_flag" == "data" ]; then
    read -p "Specify data type to download and extract [REH/GAMM_S1/GAMM_S2]: " data_type
    DATA_FLAG="--data $data_type"
  else
    DATA_FLAG="--fastq"
  fi
fi

echo "Step 3: Updating config.json"
# Check if src/config.json already exists and remove it
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi

# Copy the top-level config.json to src/config.json
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json

echo "Step 4: Preparing for Docker container run"
# Change to the output directory
cd output

# Get the absolute path to the output directory
output_dir=$(pwd)

# Go back to the previous directory
cd ..

# Get the absolute path to the shared_mount directory
shared_mount_dir="$(pwd)/shared_mount"

echo "Step 5: Building Docker image for the main container"
# Build the Docker image for the main container
docker build -t seuratv5 ./sc_pipeline

echo "Step 6: Running the main Docker container"


chmod 777 "$output_dir"
chmod 777 ./sc_pipeline/src/config.json
# Run the Docker container with the correct volume mappings
# First command for get_data.py
# First command for get_data.py
docker run -it \
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_mount_dir",target=/scRNA-seq/shared_mount \
  --mount type=bind,source="$(realpath ./sc_pipeline/src/config.json)",target=/scRNA-seq/src/config.json \
  --entrypoint=/bin/sh \
  seuratv5 -c "python3 /scRNA-seq/get_data.py $DATA_FLAG" && \

# Second command for script.R with full path to Rscript
docker run -it \
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_mount_dir",target=/scRNA-seq/shared_mount \
  --mount type=bind,source="$(realpath ./sc_pipeline/src/config.json)",target=/scRNA-seq/src/config.json \
  --entrypoint=/bin/sh \
  seuratv5 -c "/opt/conda/envs/scrnaseq/bin/Rscript /scRNA-seq/script.R"







