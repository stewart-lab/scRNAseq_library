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
  read -p "Do you need to build a genome index? [y/N]: " confirm
  if [[ "$confirm" =~ ^[Yy]$ ]]; then
    echo "Building genome index"
    apptainer run \
      --bind "$(realpath "$DATA_DIR"):/data:ro" \
      --bind "$(realpath "$CONFIG_FILE"):/config.json:ro" \
      --bind "$(realpath "$GENOME_DIR"):/genome_dir" \
      sc_aligner_v2_no_genomes_latest.sif \
      conda run -n star_env /bin/bash -c './pre_pipeline/STAR_genomebuild.sh'
  fi

  echo "Step 2.1: Removing and recreating the SHARED_MOUNT"
  
  SHARED_MOUNT="./shared_mount"
  rm -rf "$SHARED_MOUNT"
  mkdir -p "$SHARED_MOUNT"
  chmod 777 "$SHARED_MOUNT"

  echo "Step 2.2: Choose execution mode"
  read -p "Run in detached tmux session? [y/N]: " tmux_mode

  if [[ "$tmux_mode" =~ ^[Yy]$ ]]; then
    echo "Running Apptainer container for alignment in detached tmux session"
    
    # Create a tmux session with a descriptive name
    SESSION_NAME="alignment_$(date +%Y%m%d_%H%M%S)"
    
    # Create new detached tmux session and run the apptainer command
    tmux new-session -d -s "$SESSION_NAME" "
      cd $(pwd) && \
      apptainer run \
        --bind \"$(realpath "$DATA_DIR"):/data:ro\" \
        --bind \"$(realpath "$SHARED_MOUNT"):/shared_mount\" \
        --bind \"$(realpath "$CONFIG_FILE"):/config.json:ro\" \
        --bind \"$(realpath "$GENOME_INDEX_DIR"):/genome_dir\" \
        sc_aligner_v2_no_genomes_latest.sif \
        conda run -n star_env /bin/bash -c './pre_pipeline/STAR_genomedir.sh'
    "
    
    echo "Apptainer process started in tmux session: $SESSION_NAME"
    echo ""
    echo "=== TMUX COMMANDS ==="
    echo "To attach to the session: tmux attach-session -t $SESSION_NAME"
    echo "To list all sessions: tmux list-sessions"
    echo "To detach from session (when attached): Ctrl+B, then D"
    echo "To kill the session: tmux kill-session -t $SESSION_NAME"
    echo ""
    echo "=== MONITORING ==="
    echo "To monitor progress in real-time: tmux attach-session -t $SESSION_NAME"
    echo "To check if session is still running: tmux list-sessions | grep $SESSION_NAME"
  else
    echo "Running Apptainer container for alignment in foreground mode"
    apptainer run \
      --bind "$(realpath "$DATA_DIR"):/data:ro" \
      --bind "$(realpath "$SHARED_MOUNT"):/shared_mount" \
      --bind "$(realpath "$CONFIG_FILE"):/config.json:ro" \
      --bind "$(realpath "$GENOME_INDEX_DIR"):/genome_dir" \
      sc_aligner_v2_no_genomes_latest.sif \
      conda run -n star_env /bin/bash -c './pre_pipeline/STAR_genomedir.sh'
  fi
fi

echo "Step 3: Updating config.json"
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json