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

    # Optionally rebuild a STAR index from GENOME_FASTA + GENOME_GTF into GENOME_DIR/genome_build/
    USE_CUSTOM_INDEX=""
    BUILD_FAILED=""
    read -p "Do you need to build a genome index from GENOME_FASTA and GENOME_GTF? [y/N]: " build_confirm
    if [[ "$build_confirm" =~ ^[Yy]$ ]]; then
      GENOME_BUILD_DIR="${GENOME_DIR%/}/genome_build"
      if [ -d "$GENOME_BUILD_DIR" ] && [ -n "$(ls -A "$GENOME_BUILD_DIR")" ]; then
        echo "An existing genome build was found at $GENOME_BUILD_DIR"
        read -p "Rename it to keep it? (otherwise it will be overwritten) [Y/n]: " rename_confirm
        if [[ ! "$rename_confirm" =~ ^[Nn]$ ]]; then
          default_name="genome_build_$(date +%Y%m%d_%H%M%S)"
          read -p "New name [$default_name]: " new_name
          new_name=${new_name:-$default_name}
          while [ -e "${GENOME_DIR%/}/$new_name" ]; do
            read -p "${GENOME_DIR%/}/$new_name already exists, choose another name: " new_name
          done
          mv "$GENOME_BUILD_DIR" "${GENOME_DIR%/}/$new_name"
          echo "Old genome build moved to ${GENOME_DIR%/}/$new_name"
        fi
      fi

      echo "Building genome index in $GENOME_BUILD_DIR (runs in the foreground, can take an hour or more)"
      if docker run --rm \
        -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
        -v "$(realpath "$GENOME_DIR"):/genome_dir" \
        -v "$(realpath ./pre_pipeline/STAR_genomebuild.sh):/src/STAR_genomebuild.sh:ro" \
        scaligner_v2_with_genomes_and_jq /bin/bash -c "conda run --no-capture-output -n star_env /bin/bash /src/STAR_genomebuild.sh" \
        && [ -f "$GENOME_BUILD_DIR/SA" ]; then
        GENOME_INDEX_DIR="$GENOME_BUILD_DIR/"
        USE_CUSTOM_INDEX="y"
        echo "Genome build finished. For future runs, set \"GENOME_INDEX_DIR\": \"genome_build/\" in config.json to reuse it."
      else
        echo "Genome build failed, skipping alignment. Check the STAR output above."
        BUILD_FAILED="y"
      fi
    else
      read -p "Align to GENOME_INDEX_DIR ($GENOME_INDEX_DIR) instead of the built-in genomes? [y/N]: " USE_CUSTOM_INDEX
    fi

    if [ -z "$BUILD_FAILED" ]; then
      CUSTOM_INDEX_MOUNT=()
      if [[ "$USE_CUSTOM_INDEX" =~ ^[Yy]$ ]]; then
        echo "Using custom genome index: $GENOME_INDEX_DIR"
        CUSTOM_INDEX_MOUNT=(-v "$(realpath "$GENOME_INDEX_DIR"):/genome_dir:ro")
      fi

      echo "Step 2.3: Running Docker container for alignment"
      docker run -d \
        -v "$(realpath "$DATA_DIR"):/data:ro" \
        -v "$(realpath "$SHARED_MOUNT"):/shared_mount" \
        -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
        "${CUSTOM_INDEX_MOUNT[@]}" \
        scaligner_v2_with_genomes_and_jq /bin/bash -c "conda run -n star_env /bin/bash -c 'cd /src && ./STAR.sh'"
    fi
  fi
fi

echo "Step 3: Updating config.json"
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json