import argparse
import os
import shutil
import subprocess
import json

# Parse arguments
parser = argparse.ArgumentParser(description="Download and extract data.")
parser.add_argument(
    "--data",
    choices=["reh", "gamm"],
    default=None,
    help="Specify data type to download and extract.",
)
args = parser.parse_args()

# Only execute if the --data argument is 'reh'
if args.data == "reh":
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Make the DATA directory
    data_dir = os.path.join(script_dir, "DATA")
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
            "https://www.morgridge.net/dmz/jfreeman/REH_DATA.tgz",
        ],
        check=True,
    )

    # Move the downloaded file to the DATA directory
    source_file = os.path.join(
        data_dir, "www.morgridge.net", "dmz", "jfreeman", "REH_DATA.tgz"
    )
    shutil.move(source_file, "./REH_DATA.tgz")

    # Remove the extra directories
    shutil.rmtree(os.path.join(data_dir, "www.morgridge.net"))

    # Extract the tar.gz file
    subprocess.run(["tar", "-zxvf", "REH_DATA.tgz"], check=True)

    # Remove the tar.gz file
    os.remove("REH_DATA.tgz")

    # Change back to the script directory
    os.chdir(script_dir)

    # Load the config.json file
    with open(os.path.join(script_dir, "config.json"), "r") as f:
        config = json.load(f)

    # Recursively search through the directories
    gene_full_dirs = []
    for root, dirs, files in os.walk(data_dir):
        for dir in dirs:
            if dir.endswith("GeneFull"):
                gene_full_dirs.append(
                    {"name": dir, "base_directory": os.path.join(root, dir)}
                )

    # Adjust the 'lanes' list in the config to match the list of directories
    config["lanes"] = gene_full_dirs

    # Dump the updated config back to the file
    with open(os.path.join(script_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=4)
