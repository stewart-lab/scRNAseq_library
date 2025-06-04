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
if config["process_matrix"] == True:
    for dir in os.listdir(search_dir):
        gene_full_path = os.path.join("/scRNA-seq", search_dir, dir) #"/scRNA-seq",
        #print(gene_full_path)
        match = re.match(r"(.+?_lane\d+)", gene_full_path.split("/")[-1])
        name = match.group(1) if match else gene_full_path.split("/")[-1]
        print(name)
        gene_full_dirs.append({"name": name, "base_directory": gene_full_path})
        
        # Process files in GeneFull directory
        for file in os.listdir(gene_full_path):
            file_path = os.path.join(gene_full_path, file)
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
                match = re.match(r"(.+?_lane\d+)", root.split("/")[-2])
                name = match.group(1) if match else root.split("/")[-2]
                gene_full_dirs.append({"name": name, "base_directory": gene_full_path})
            
                # Process files in GeneFull directory
                for subroot, subdirs, subfiles in os.walk(gene_full_path):
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

