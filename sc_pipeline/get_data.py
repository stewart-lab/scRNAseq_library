import os
import shutil
import json
import warnings
import gzip
import re

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Open and read the existing config.json
config_path = os.path.join(script_dir, "src/config.json")
with open(config_path, "r") as f:
    config = json.load(f)

config["lanes"] = []
gene_full_dirs = []

# Process directories
search_dir = "./shared_mount"
if config.get("parse", False) == True:
    # Parse Biosciences (split-pipe) layout: <lane>/<sample>/DGE_filtered|DGE_unfiltered/
    # e.g. shared_mount/<OUT_DIR>/S39/H1D1/DGE_filtered/count_matrix.mtx
    # OUT_DIR may point directly at one run's lanes, or at a parent folder holding
    # several timestamped run copies - walk the whole tree either way and rely on
    # the declared lane names (not directory depth) to pick out the real lanes.
    fastq_alignment_parse = config.get("fastq_alignment_parse", {})
    parse_out_dir = fastq_alignment_parse.get("OUT_DIR", "")
    parse_root = os.path.join(search_dir, parse_out_dir)

    # Only accept lanes explicitly declared as a RUN_* NAME in the config, so
    # unrelated/stray directories anywhere under OUT_DIR are never picked up.
    expected_lane_names = {
        run["NAME"] for run in fastq_alignment_parse.values()
        if isinstance(run, dict) and "NAME" in run
    }
    if not expected_lane_names:
        warnings.warn("No RUN_* NAME entries found in fastq_alignment_parse; check your config file.")

    # If a declared lane/sample shows up more than once (e.g. leftover copies in
    # older timestamped run folders), keep only the most recently modified one.
    latest_by_lane_sample = {}
    for root, dirs, files in os.walk(parse_root):
        # skip the aggregated cross-lane output, we only want per-lane data
        if "combined" in root.split(os.sep):
            continue
        if "DGE_filtered" in dirs and "DGE_unfiltered" in dirs:
            lane_name = os.path.basename(os.path.dirname(root))
            if lane_name not in expected_lane_names:
                print(f"Skipping {root}: lane '{lane_name}' not declared in fastq_alignment_parse")
                continue
            sample_name = os.path.basename(root)
            mtime = os.path.getmtime(root)
            key = (lane_name, sample_name)
            if key not in latest_by_lane_sample or mtime > latest_by_lane_sample[key][0]:
                latest_by_lane_sample[key] = (mtime, root)

    if not latest_by_lane_sample:
        warnings.warn(f"No lanes matching fastq_alignment_parse NAMEs found under {parse_root}")

    for (lane_name, sample_name), (mtime, root) in latest_by_lane_sample.items():
        full_path = os.path.normpath(os.path.join("/scRNA-seq", root))
        name = f"{lane_name}_{sample_name}"
        print(name)
        gene_full_dirs.append({"name": name, "base_directory": full_path})
elif config["process_matrix"] == True:
    for dir in os.listdir(search_dir):
        gene_full_path = os.path.join("/scRNA-seq", search_dir, dir) #"/scRNA-seq",
        #print(gene_full_path)
        clean_path = os.path.normpath(gene_full_path)
        match = re.match(r"(.+?_lane\d+)", clean_path.split("/")[-1])
        name = match.group(1) if match else clean_path.split("/")[-1]
        print(name)
        gene_full_dirs.append({"name": name, "base_directory": clean_path})
        
        # Process files in GeneFull directory
        for file in os.listdir(clean_path):
            file_path = os.path.join(clean_path, file)
            if file.startswith("UniqueAndMult") or file.endswith(".gz.gz"):
                os.remove(file_path)
                continue
                    
            if (file.endswith(".tsv") or file.startswith("matrix")) and not file.endswith(".gz"):
                gz_file_path = f"{file_path}.gz"
                if not os.path.exists(gz_file_path):
                    try:
                        with open(file_path, "rb") as f_in, gzip.open(gz_file_path, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        print(f"Successfully gzipped {file_path}")
                    except Exception as e:
                        print(f"Failed to gzip {file_path}: {e}")
else:        
    for root, dirs, files in os.walk(search_dir):
        for dir in dirs:
            if dir == "GeneFull":
                gene_full_path = os.path.join("/scRNA-seq", root, dir)
                clean_path = os.path.normpath(gene_full_path)
                # Extract lane name from the path - go up two levels from GeneFull
                lane_dir = root.split("/")[-2] if len(root.split("/")) > 1 else root.split("/")[-1]
                match = re.match(r"(.+?_lane\d+)", lane_dir)
                name = match.group(1) if match else lane_dir
                gene_full_dirs.append({"name": name, "base_directory": clean_path})
            
                # Process files in GeneFull directory
                for subroot, subdirs, subfiles in os.walk(clean_path):
                    for file in subfiles:
                        file_path = os.path.join(subroot, file)
                        if file.startswith("UniqueAndMult") or file.endswith(".gz.gz"):
                            os.remove(file_path)
                            continue
                    
                        if (file.endswith(".tsv") or file.startswith("matrix")) and not file.endswith(".gz"):
                            gz_file_path = f"{file_path}.gz"
                            if not os.path.exists(gz_file_path):
                                try:
                                    with open(file_path, "rb") as f_in, gzip.open(gz_file_path, "wb") as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                    print(f"Successfully gzipped {file_path}")
                                except Exception as e:
                                    print(f"Failed to gzip {file_path}: {e}")

if not gene_full_dirs:
    warnings.warn("No GeneFull directories found. The config.json will be empty for lanes.")

# Update config.json
config["lanes"] = gene_full_dirs
with open(os.path.join(script_dir, "src/config.json"), "w") as f:
    json.dump(config, f, indent=4)

