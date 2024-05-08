#!/bin/bash

echo "Step 1: Importing DATA_DIR from config.json"
# Import DATA_DIR from config.json
CONFIG_FILE="./config.json"
DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

echo "Step 2: Confirming data alignment"
# Confirm with the user if they have loaded new data or would like to realign
read -p "Have you loaded new data or would you like to realign? [y/N]: " confirm

data_flag="fastq"  # Default to --fastq

if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
  echo "Step 2.1: Removing and recreating the SHARED_VOLUME"
  # Define the shared volume path
  SHARED_VOLUME="./shared_volume"

  rm -rf "$SHARED_VOLUME"

  # Recreate the SHARED_VOLUME directory
  mkdir -p "$SHARED_VOLUME"

  # Your next steps here
  DATA_FLAG="--fastq"

  echo "Step 2.2: Building Docker image for alignment"
  # Build the Docker image from the pre_pipeline directory
  docker build --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) -t align_fastqs ./pre_pipeline

  echo "Step 2.3: Running Docker container for alignment"
  # Run the Docker container and execute the STAR.sh script
docker run -it \
    --user=$(id -u):$(id -g) \
    --mount type=bind,source="$DATA_DIR",target=/data \
    --mount type=bind,source="$SHARED_VOLUME",target=/shared_volume \
    --mount type=bind,source="$CONFIG_FILE",target=/config.json \
    align_fastqs /bin/bash -c "conda run -n star_env /bin/bash -c 'cd /src && ./STAR.sh'"
else
  echo "Skipped align_fastqs."
  
  echo "Step 2.4: Asking for data flag"
  # Ask the user which flag to use for get_data.py
  read -p "If you'd like to load a stored experiment select data. If you have aligned FASTQs loaded and changed pipeline parameters, select fastq [data/fastq]: " data_flag
  
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
cp $CONFIG_FILE ./sc_pipeline/src/config.json

# change the ownership of the sc_pipeline directory to the current user
chown -R $(id -u):$(id -g) ./sc_pipeline/

echo "Step 4: Preparing for Docker container run"
# Change to the output directory
cd output

# Get the absolute path to the current directory
output_dir=$(pwd)

# Go back to the previous directory
cd ..

# Get the absolute path to the shared_volume directory
shared_volume_dir=$(pwd)/shared_volume
echo "Step 5: Building Docker image for the next container"
# Build the Docker image for the next container
docker build --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) -t seuratv5 ./sc_pipeline

echo "Step 6: Running the main Docker container"

# Run the Docker container with the correct volume mapping
docker run -it\
  --user=$(id -u):$(id -g)\
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_volume_dir",target=/scRNA-seq/shared_volume \
  --mount type=bind,source="./sc_pipeline/src/config.json",target=/scRNA-seq/src/config.json \
  seuratv5 /bin/bash -c "python /scRNA-seq/get_data.py $DATA_FLAG && /root/miniconda/bin/conda run -n scrnaseq Rscript /scRNA-seq/script.R"
  
echo "Pipeline completed successfully."













