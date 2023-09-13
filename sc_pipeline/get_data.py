import argparse
import os
import shutil
import subprocess
import json
import warnings
import gzip

# Parse arguments
parser = argparse.ArgumentParser(description="Download and extract data.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "--data",
    choices=["REH", "GAMM_S1", "GAMM_S2"],
    default=None,
    help="Specify data type to download and extract.",
)
group.add_argument(
    "--fastq", action="store_true", help="Look in the shared_volume directory."
)
args = parser.parse_args()

# Mapping of data types to their download links
download_links = {
    "REH": "https://www.morgridge.net/dmz/jfreeman/REH_DATA.tar.gz",
    "GAMM_S1": "https://www.morgridge.net/dmz/jfreeman/GAMM_S1_DATA.tar.gz",
    "GAMM_S2": "https://www.morgridge.net/dmz/jfreeman/GAMM_S2_DATA.tar.gz",
}

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Open and read the existing config.json
config_path = os.path.join(script_dir, "src/config.json")
with open(config_path, "r") as f:
    config = json.load(f)

config["lanes"] = []

gene_full_dirs = []

def search_gene_full_dirs(search_dir):
    for root, dirs, files in os.walk(search_dir):
        for dir in dirs:
            if dir.endswith("GeneFull"):
                gene_full_path = os.path.join("/scRNA-seq", root, dir)
                gene_full_dirs.append(
                    {"name": dir, "base_directory": gene_full_path}
                )


if args.data:
    assert (
        args.data in download_links
    ), "Please specify a correct --data argument. Use one of the following options: ['REH', 'GAMM_S1', 'GAMM_S2']"

    data_dir = os.path.join(script_dir, "DATA")

    # Check if the directory already exists and remove it if it does
    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)

    # Make the DATA directory
    os.makedirs(data_dir, exist_ok=True)

    # Change directory to the newly created DATA directory
    os.chdir(data_dir)

    # Use wget to download the file
    subprocess.run(
        [
            "wget",
            "-r",
            "--no-check-certificate",
            "--no-parent",
            download_links[args.data],
        ],
        check=True,
    )

    # Move the downloaded file to the DATA directory
    source_file = os.path.join(
        data_dir, "www.morgridge.net", "dmz", "jfreeman", f"{args.data}_DATA.tar.gz"
    )
    destination_file = os.path.join(data_dir, f"{args.data}_DATA.tar.gz")
    shutil.move(source_file, destination_file)

    # Remove the extra directories
    shutil.rmtree(os.path.join(data_dir, "www.morgridge.net"))

    # Extract the tar.gz file
    subprocess.run(["tar", "-zxvf", f"{args.data}_DATA.tar.gz"], check=True)

    # Get the path of the nested directory
    nested_dir = os.path.join(data_dir, "dmz", "jfreeman", f"{args.data}_DATA")

    # Check if the nested directory exists
    if os.path.isdir(nested_dir):
        # Move each file in the nested directory to the DATA directory
        for filename in os.listdir(nested_dir):
            dest_file_path = os.path.join(data_dir, filename)
            if os.path.isfile(dest_file_path):
                os.remove(dest_file_path)
            shutil.move(os.path.join(nested_dir, filename), data_dir)

        # Remove the extra directories
        shutil.rmtree(os.path.join(data_dir, "dmz"))

    # Remove the tar.gz file
    os.remove(f"{args.data}_DATA.tar.gz")

    # Change back to the script directory
    os.chdir(script_dir)
    search_gene_full_dirs(data_dir)
    # Update species based on data
    if args.data in ["GAMM_S1", "GAMM_S2"]:
        config["filter_cells"]["species"] = "pig"
    elif args.data == "REH":
        config["filter_cells"]["species"] = "human"

# Define the directory to search based on the --fastq flag
search_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "DATA")
if args.fastq:
    search_dir = "./shared_volume"
    for root, dirs, files in os.walk(search_dir):
        for dir in dirs:
            if dir == "GeneFull":
                gene_full_path = os.path.join("/scRNA-seq", root, dir)  # Changed path
                gene_full_dirs.append(
                    {"name": dir, "base_directory": gene_full_path}
                )
                for subroot, subdirs, subfiles in os.walk(gene_full_path):
                    for file in subfiles:
                        if file.endswith(".tsv") or file.endswith(".mtx"):
                            file_path = os.path.join(subroot, file)
                            if not os.path.exists(f"{file_path}.gz"):  # Check if already gzipped
                                try:
                                    with open(file_path, "rb") as f_in, gzip.open(
                                        f"{file_path}.gz", "wb"
                                    ) as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                    print(f"Successfully gzipped {file_path}")
                                except Exception as e:
                                    print(f"Failed to gzip {file_path}: {e}")
                            else:
                                print(f"{file_path} is already gzipped, skipping.")

# Adjust the 'lanes' list in the config to match the list of directories
config["lanes"] = gene_full_dirs

if not gene_full_dirs:
    warnings.warn(
        "No GeneFull directories found. The config.json will be empty for lanes."
    )

# Update the config.json
with open(os.path.join(script_dir, "src/config.json"), "r") as f:
    config = json.load(f)
config["lanes"] = gene_full_dirs
with open(os.path.join(script_dir, "src/config.json"), "w") as f:
    json.dump(config, f, indent=4)
