#!/bin/bash
# usage: source run_seuratv5.sh [config.json]   (defaults to ./config.json)
#
# Species subsetting for human+pig chimeric data (README_chimeric_analysis.md):
# give a fastq_alignment sample SPECIES_BARCODES (barcodes_<species>.txt from
# 03_classify_species.py; a path or {"lane1": path, ...}) and SOUPX_RHO (a
# number or {"lane1": rho, ...}).  SoupX then runs on the whole library with
# that rho, and the corrected counts are subset to those barcodes before all
# downstream steps.  The barcode files usually live outside shared_mount, so
# their folders are mounted read-only into the container at the same path.
CONFIG_FILE="${1:-./config.json}"
if [ ! -f "$CONFIG_FILE" ]; then
  echo "Config file $CONFIG_FILE not found"
  return 1 2>/dev/null || exit 1
fi

echo "Step 3: Updating config.json"
if [ -f "./sc_pipeline/src/config.json" ]; then
  rm ./sc_pipeline/src/config.json
fi
cp "$CONFIG_FILE" ./sc_pipeline/src/config.json

echo "Step 3.1: Checking species-subset settings (SPECIES_BARCODES / SOUPX_RHO)"
# prints one directory per line to mount; exits non-zero on a bad config
BARCODE_DIRS=$(python3 - "$CONFIG_FILE" <<'EOF'
import json, os, sys
cfg = json.load(open(sys.argv[1]))
dirs, errors = set(), []
for key, s in cfg.get("fastq_alignment", {}).items():
    if not isinstance(s, dict):
        continue
    bc, rho = s.get("SPECIES_BARCODES"), s.get("SOUPX_RHO")
    if bc is None and rho is None:
        continue
    name = s.get("NAME", key)
    if str(s.get("COMBINE_LANES", False)).lower() == "true":
        for k, v in (("SPECIES_BARCODES", bc), ("SOUPX_RHO", rho)):
            if isinstance(v, dict):
                errors.append(f"{name}: COMBINE_LANES is true, so {k} must be a single value, "
                              f"not a per-lane map")
    if bc is not None and s.get("species") not in ("pig", "human"):
        errors.append(f"{name}: species must be 'pig' or 'human' (the species kept), "
                      f"not {s.get('species')!r}")
    bcs = bc if isinstance(bc, dict) else ({"all lanes": bc} if bc else {})
    for lane, p in bcs.items():
        if not os.path.isabs(p):
            errors.append(f"{name} {lane}: SPECIES_BARCODES must be an absolute path: {p}")
        elif not os.path.isfile(p):
            errors.append(f"{name} {lane}: SPECIES_BARCODES not found: {p}")
        else:
            n = sum(1 for l in open(p) if l.strip())
            print(f"  {name} {lane}: {n} barcodes from {p}", file=sys.stderr)
            dirs.add(os.path.dirname(p))
    rhos = rho if isinstance(rho, dict) else ({"all lanes": rho} if rho is not None else {})
    for lane, r in rhos.items():
        try:
            ok = 0 < float(r) < 1
        except (TypeError, ValueError):
            ok = False
        if not ok:
            errors.append(f"{name} {lane}: SOUPX_RHO must be a number between 0 and 1, got {r!r}")
        else:
            print(f"  {name} {lane}: SoupX rho = {float(r)}", file=sys.stderr)
    if isinstance(bc, dict) and isinstance(rho, dict) and set(bc) != set(rho):
        errors.append(f"{name}: SPECIES_BARCODES lanes {sorted(bc)} != SOUPX_RHO lanes {sorted(rho)}")
for e in errors:
    print("ERROR: " + e, file=sys.stderr)
if errors:
    sys.exit(1)
print("\n".join(sorted(dirs)))
EOF
) || { echo "Fix the species-subset settings in $CONFIG_FILE"; return 1 2>/dev/null || exit 1; }

BARCODE_MOUNTS=()
while IFS= read -r d; do
  [ -n "$d" ] || continue
  # mount at the same path the config uses, so the R code can read it as-is
  BARCODE_MOUNTS+=(--mount "type=bind,source=$(realpath "$d"),target=$d,readonly")
done <<< "$BARCODE_DIRS"
if [ ${#BARCODE_MOUNTS[@]} -gt 0 ]; then
  echo "  mounting ${#BARCODE_MOUNTS[@]} barcode folder(s) read-only"
else
  echo "  none set: standard run (STARsolo filtered cells, SoupX autoEstCont)"
fi

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
docker run -d \
  --memory="400g" \
  --memory-swap="400g" \
  --mount type=bind,source="$output_dir",target=/scRNA-seq/output \
  --mount type=bind,source="$shared_mount_dir",target=/scRNA-seq/shared_mount \
  --mount type=bind,source="$(realpath ./sc_pipeline/src/config.json)",target=/scRNA-seq/src/config.json \
  "${BARCODE_MOUNTS[@]}" \
  --entrypoint=/bin/sh \
  seuratv5 -c "python3 /scRNA-seq/get_data.py $DATA_FLAG && /opt/conda/envs/scrnaseq/bin/Rscript /scRNA-seq/script.R"
