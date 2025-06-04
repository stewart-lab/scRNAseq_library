#!/bin/bash
CONFIG_FILE="./config.json"
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
DATA_FLAG="--fastq"  # Default to --fastq

echo "Step 5: Building Docker image for the main container"
docker build -t seuratv5 ./sc_pipeline

echo "Step 6: Running the main Docker container"
chmod 777 "$output_dir"
chmod 777 ./sc_pipeline/src/config.json

# Run pipeline with simplified commands #--memory-swap="64g" \
docker run -it \
  --memory="100g" \
  --memory-swap="150g" \
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_mount_dir",target=/scRNA-seq/shared_mount \
  --mount type=bind,source="$(realpath ./sc_pipeline/src/config.json)",target=/scRNA-seq/src/config.json \
  --entrypoint=/bin/sh \
  seuratv5 -c "python3 /scRNA-seq/get_data.py $DATA_FLAG && /opt/conda/envs/scrnaseq/bin/Rscript /scRNA-seq/script.R" \
