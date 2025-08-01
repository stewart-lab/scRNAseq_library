#!/bin/bash

# Path to the configuration file
CONFIG_FILE="/config.json"

umask 000

# define genome dir
GENOME_DIR="/genome_dir"
# mk output dir for genome build
mkdir -p $GENOME_DIR/genome_build/
GENOME_BUILD_DIR=$GENOME_DIR/genome_build/

# Ensure the shared volume is writable
if [ ! -w $GENOME_BUILD_DIR ]; then
    echo "Cannot write to /genome_build/. Please check permissions."
    exit 1
fi

# Python function to extract JSON values
get_json_value() {
    local key="$1"
    python3 -c "import json, sys; data = json.load(open('$CONFIG_FILE')); print(data.get('$key', ''))"
}

# define genome fasta
GENOME_FASTA=$(get_json_value "GENOME_FASTA")
GENOME_FASTA_DIR=$GENOME_DIR/$GENOME_FASTA
# define genome gtf
GENOME_GTF=$(get_json_value "GENOME_GTF")
GENOME_GTF_DIR=$GENOME_DIR/$GENOME_GTF
# define run thread n
RUN_THREAD_N=$(get_json_value "RUN_THREAD_N")

# Run STAR command with options for the current sample and lane
STAR --runMode genomeGenerate --genomeDir $GENOME_BUILD_DIR \
    --genomeFastaFiles $GENOME_FASTA_DIR \
    --sjdbGTFfile $GENOME_GTF_DIR    \
    --sjdbOverhang 149 \
    --runThreadN $RUN_THREAD_N \