#!/bin/bash
CONFIG_FILE="./config.json"
echo "Step 3: Updating config.json"
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json

echo "Step 4: Preparing for Apptainer run"
cd output
output_dir=$(pwd)
cd ..
shared_mount_dir="$(pwd)/shared_mount"
chmod 777 "$output_dir"
chmod 777 ./sc_pipeline/src/config.json

if [[ "$tmux_mode" =~ ^[Yy]$ ]]; then
    echo "Running Apptainer container for alignment in detached tmux session"
    
    # Create a tmux session with a descriptive name
    SESSION_NAME="alignment_$(date +%Y%m%d_%H%M%S)"
    
    # Create new detached tmux session and run the apptainer command
    tmux new-session -d -s "$SESSION_NAME" "
      cd $(pwd) && \
      apptainer exec --cleanenv \
      --env LANG=en_US.UTF-8 \
      --env LC_ALL=en_US.UTF-8 \
      --bind "$(realpath "$output_dir"):/scRNA-seq/output" \
      --bind "$(realpath "$shared_mount_dir"):/scRNA-seq/shared_mount" \
      --bind "$(realpath "$CONFIG_FILE"):/scRNA-seq/src/config.json" \
      scrnaseq-env_latest.sif \
      /bin/bash -lc 'source /opt/conda/etc/profile.d/conda.sh && conda activate scrnaseq && python3 sc_pipeline/get_data_apptainer.py && Rscript sc_pipeline/script_apptainer.R'
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
    apptainer exec --cleanenv \
      --env LANG=en_US.UTF-8 \
      --env LC_ALL=en_US.UTF-8 \
      --bind "$(realpath "$output_dir"):/scRNA-seq/output" \
      --bind "$(realpath "$shared_mount_dir"):/scRNA-seq/shared_mount" \
      --bind "$(realpath "$CONFIG_FILE"):/scRNA-seq/src/config.json" \
      scrnaseq-env_latest.sif \
      /bin/bash -lc 'source /opt/conda/etc/profile.d/conda.sh && conda activate scrnaseq && python3 sc_pipeline/get_data_apptainer.py && Rscript sc_pipeline/script_apptainer.R'
fi