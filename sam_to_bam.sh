#!/bin/bash
# Convert STAR Aligned.out.sam files to coordinate-sorted, indexed BAMs.
#
# usage: bash sam_to_bam.sh [--delete] [--threads N] [file.sam ...]
#   no files given -> every Aligned.out.sam under ./shared_mount
#   --delete       -> remove each SAM after its BAM passes the checks below
#
# Run inside the star_env conda env (samtools 1.17):  conda activate star_env

THREADS=16
DELETE=false
SAMS=()

while [ $# -gt 0 ]; do
  case "$1" in
    --delete) DELETE=true ;;
    --threads) THREADS="$2"; shift ;;
    *) SAMS+=("$1") ;;
  esac
  shift
done

if ! command -v samtools >/dev/null; then
  echo "samtools not found -- run 'conda activate star_env' first" >&2
  exit 1
fi

if [ ${#SAMS[@]} -eq 0 ]; then
  mapfile -t SAMS < <(find ./shared_mount -name "Aligned.out.sam" | sort)
fi

for SAM in "${SAMS[@]}"; do
  DIR=$(dirname "$SAM")
  BAM="$DIR/Aligned.sortedByCoord.out.bam"
  echo "[$(date '+%F %T')] $SAM -> $BAM"

  if ! samtools sort -@ "$THREADS" -m 2G -T "$DIR/sort_tmp" -o "$BAM" "$SAM"; then
    echo "  FAILED: samtools sort, keeping SAM" >&2
    rm -f "$BAM" "$DIR"/sort_tmp*.bam
    continue
  fi
  samtools index -@ "$THREADS" "$BAM"

  # Verify: BAM is intact and holds the same number of records as the SAM
  SAM_N=$(grep -vc '^@' "$SAM")
  BAM_N=$(samtools view -c -@ "$THREADS" "$BAM")
  if samtools quickcheck "$BAM" && [ "$SAM_N" -eq "$BAM_N" ]; then
    echo "  OK: $BAM_N records, $(du -h "$SAM" | cut -f1) SAM -> $(du -h "$BAM" | cut -f1) BAM"
    if [ "$DELETE" = true ]; then
      rm "$SAM"
      echo "  deleted $SAM"
    fi
  else
    echo "  MISMATCH: SAM $SAM_N vs BAM $BAM_N records, keeping SAM" >&2
  fi
done

echo "[$(date '+%F %T')] done"
