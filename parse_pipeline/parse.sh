#!/bin/bash

# Path to the configuration file
CONFIG_FILE="/config.json"

umask 000

# Ensure the shared volume is writable
if [ ! -w /shared_mount ]; then
    echo "Cannot write to /shared_volume. Please check permissions."
    exit 1
fi

# Fetch all run names from the config file
NUMBER_OF_RUNS=$(jq -r '.fastq_alignment_parse.Number_of_runs' $CONFIG_FILE)

# create output
# Create output directory in shared volume
timestamp=$(date +%Y%m%d_%H%M%S)
OUT_DIR=$(jq -r ".fastq_alignment_parse.OUT_DIR" $CONFIG_FILE)_${timestamp}
mkdir -p /shared_mount/$OUT_DIR

for ((RUN_IDX=1; RUN_IDX<=NUMBER_OF_RUNS; RUN_IDX++)); do
  # Loop over each run based on NUMBER_OF_RUNS
  RUN_KEY="RUN_$RUN_IDX"  # Construct the key to access run configurations
  NAME=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.NAME" $CONFIG_FILE)
  MODE=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.MODE" $CONFIG_FILE)
  CHEM=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.CHEM" $CONFIG_FILE)
  cDNA_LANE=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.cDNA_LANE" $CONFIG_FILE)
  BARCODE_LANE=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.BARCODE_LANE" $CONFIG_FILE)
  KIT=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.KIT" $CONFIG_FILE)  
  SAMPLE_LIST=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.SAMPLE_LIST" $CONFIG_FILE)
  NTHREADS=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.RUN_THREAD_N" $CONFIG_FILE)
  KIT_SCORE_SKIP=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.KIT_SCORE_SKIP" $CONFIG_FILE)
  START_TIMEOUT=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.START_TIMEOUT" $CONFIG_FILE)
  
  mkdir -p /shared_mount/$OUT_DIR/$NAME

  split-pipe \
    --mode $MODE \
    --chemistry $CHEM \
    --kit $KIT \
    --genome_dir /genome_index_dir \
    --fq1 $cDNA_LANE \
    --output_dir /shared_mount/$OUT_DIR/$NAME \
    --samp_list /data/$SAMPLE_LIST \
    --nthreads $NTHREADS --kit_score_skip --start_timeout $START_TIMEOUT

    # Log the completion of each run
    echo "Completed processing $NAME" >> /shared_mount/alignment_log.txt

done

# combine runs
echo "Combining runs"
# fetch runs to combine from config file
NUMBER_OF_RUNS=$(jq -r '.fastq_alignment_parse.Number_of_runs' $CONFIG_FILE)
SUB_LIBS=()
for ((RUN_IDX=1; RUN_IDX<=NUMBER_OF_RUNS; RUN_IDX++)); do
    RUN_KEY="RUN_$RUN_IDX"  # Construct the key to access run configurations
    NAME=$(jq -r ".fastq_alignment_parse.${RUN_KEY}.NAME" $CONFIG_FILE)
    # add each run's output dir to the array for --sublibraries
    SUB_LIBS+=("/shared_mount/$OUT_DIR/$NAME")

done

mkdir -p /shared_mount/$OUT_DIR/combined

split-pipe \
--mode comb \
--output_dir /shared_mount/$OUT_DIR/combined \
--sublibraries "${SUB_LIBS[@]}"

echo "finished parse alignments and combining runs" >> /shared_mount/alignment_log.txt