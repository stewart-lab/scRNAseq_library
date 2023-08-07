import os
import argparse
import subprocess
from pysradb import SRAweb
import pandas as pd


def get_srp_from_gse(gse):
    db = SRAweb()
    try:
        srp_df = db.gse_to_srp(gse)
        if srp_df.empty:
            print(f"Error: No SRP found for GSE: {gse}")
            return None
        else:
            return srp_df["study_accession"].iloc[0]  # get the first SRP
    except Exception as e:
        print("Unexpected error occurred:", e)
        return None


def get_srx_from_gsm(gsm):
    db = SRAweb()
    try:
        srx_df = db.gsm_to_srx(gsm)
        if srx_df.empty:
            print(f"Error: No SRX found for GSM: {gsm}")
            return None
        else:
            return srx_df["experiment_accession"].tolist()  # get all SRXs
    except Exception as e:
        print("Unexpected error occurred:", e)
        return None


def download_fastq_from_srp(srp, out_dir):
    db = SRAweb()
    try:
        df = db.sra_metadata(srp, detailed=True)
        os.chdir(out_dir)
        for _, row in df.iterrows():
            print(f"Downloading {row['run_accession']}")
            download_df = pd.DataFrame(row).transpose()
            db.download(df=download_df, skip_confirmation=True)
            convert_sra_to_fastq(row["run_accession"])
    except Exception as e:
        print("Error occurred while downloading:", e)
    finally:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))


def convert_sra_to_fastq(sra_file):
    try:
        command = f"parallel-fastq-dump -s {sra_file} --split-files --gzip --outdir ."
        subprocess.check_call(command, shell=True)
        print(f"Converted {sra_file} to FASTQ format.")
    except Exception as e:
        print(f"Error occurred while converting {sra_file} to FASTQ: ", e)


def download_fastq_from_srx(srx, out_dir):
    db = SRAweb()
    try:
        df = db.sra_metadata(srx, detailed=True)
        os.chdir(out_dir)
        for _, row in df.iterrows():
            print(f"Downloading {row['run_accession']}")
            download_df = pd.DataFrame(row).transpose()
            db.download(df=download_df, skip_confirmation=True)
            convert_sra_to_fastq(row["run_accession"])
    except Exception as e:
        print("Error occurred while downloading:", e)
    finally:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))


def gse_or_gsm_to_fastq(out_dir, gse=None, gsm=None):
    if gse is not None:
        srp = get_srp_from_gse(gse)
        if srp is not None:
            download_fastq_from_srp(srp, out_dir)
    elif gsm is not None:
        srx_list = get_srx_from_gsm(gsm)
        if srx_list is not None:
            for srx in srx_list:
                download_fastq_from_srx(srx, out_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download and convert SRA files from a GSE accession."
    )
    parser.add_argument("--gse", help="The GSE accession.")
    parser.add_argument("--gsm", help="The GSM accession.")
    parser.add_argument(
        "out_dir",
        help="The directory to output the FASTQ files. This is a required argument.",
    )
    args = parser.parse_args()

    if args.gse is None and args.gsm is None:
        parser.error("At least one of --gse or --gsm must be provided.")

    gse_or_gsm_to_fastq(gse=args.gse, gsm=args.gsm, out_dir=args.out_dir)
