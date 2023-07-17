import argparse
import os
import shutil
import subprocess
import json


# Parse arguments
parser = argparse.ArgumentParser(description="Download and extract data.")
parser.add_argument(
    "--data",
    choices=["REH", "GAMM_S1", "GAMM_S2"],
    default=None,
    help="Specify data type to download and extract.",
)
args = parser.parse_args()

# Mapping of data types to their download links
download_links = {
    "REH": "https://www.morgridge.net/dmz/jfreeman/REH_DATA.tar.gz",
    "GAMM_S1": "https://www.morgridge.net/dmz/jfreeman/GAMM_S1_DATA.tar.gz",
    "GAMM_S2": "https://www.morgridge.net/dmz/jfreeman/GAMM_S2_DATA.tar.gz",
}
if args.data and args.data not in download_links:
    print(
        "Data flag given but argument doesn't match any of our projects. Please see README.Md"
    )
    exit(1)

# If --data argument is not provided or doesn't match any of the choices, script execution will stop here
assert (
    args.data in download_links
), "Please specify a correct --data argument. Use one of the following options: ['REH', 'GAMM_S1', 'GAMM_S2']"

# Only execute if the --data argument is in the download_links
if args.data in download_links:
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
            download_links[args.data],
        ],
        check=True,
    )

    # Move the downloaded file to the DATA directory
    # Move the downloaded file to the DATA directory
    source_file = os.path.join(
        data_dir, "www.morgridge.net", "dmz", "jfreeman", f"{args.data}_DATA.tar.gz"
    )
    destination_file = os.path.join(data_dir, f"{args.data}_DATA.tar.gz")
    shutil.move(source_file, destination_file)

    # Remove the extra directories
    shutil.rmtree(os.path.join(data_dir, "www.morgridge.net"))

    # Extract the tar.gz file
    # ...

    # Extract the tar.gz file
    subprocess.run(["tar", "-zxvf", f"{args.data}_DATA.tar.gz"], check=True)

    # Get the path of the nested directory
    nested_dir = os.path.join(data_dir, "dmz", "jfreeman", f"{args.data}_DATA")

    # Check if the nested directory exists
    if os.path.isdir(nested_dir):
        # Move each file in the nested directory to the DATA directory
        for filename in os.listdir(nested_dir):
            shutil.move(os.path.join(nested_dir, filename), data_dir)

        # Remove the extra directories
        shutil.rmtree(os.path.join(data_dir, "dmz"))

    # ...

    # Remove the tar.gz file
    os.remove(f"{args.data}_DATA.tar.gz")

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
