import pandas as pd
from scipy.io import mmread, mmwrite
import random
import numpy as np
from pathlib import Path
import gzip
from sys import argv

DATA = Path("Data")
S1_1 = DATA / "S1_1"
OUT = DATA / "default"
Path.mkdir(OUT, exist_ok=True)


def read_data(indir=S1_1):
    try:
        barcodes = pd.read_csv(indir / "barcodes.tsv.gz", header=None)
        features = pd.read_csv(indir / "features.tsv.gz", header=None)
        with gzip.open(indir / "matrix.mtx.gz", "rb") as f:
            mtx = mmread(f)
    except FileNotFoundError:
        barcodes = pd.read_csv(indir / "barcodes.tsv", header=None)
        features = pd.read_csv(indir / "features.tsv", header=None)
        mtx = mmread(indir / "matrix.mtx").tocsc()
    return (barcodes, features, mtx)


def downsample(barcodes, features, mtx):
    # Set the seed for reproducibility
    random.seed(123)
    np.random.seed(123)

    # Select 100 random barcodes and features
    selected_barcodes = barcodes.sample(100)

    sum_by_feature = np.array(mtx.sum(axis=1)).flatten()
    top_features_indices = sum_by_feature.argsort()[-100:]
    selected_features = features.iloc[top_features_indices]

    # Subset the matrix to the selected barcodes and features
    selected_mtx = mtx[top_features_indices][:, selected_barcodes.index]
    return (selected_barcodes, selected_features, selected_mtx)


def write_data(barcodes, features, mtx, outdir=OUT):
    # Write the downsampled barcodes and features to new .tsv files
    barcodes.to_csv(
        DATA / "barcodes.tsv.gz", index=False, header=False, compression=gzip
    )
    features.to_csv(
        DATA / "features.tsv.gz", index=False, header=False, compression=gzip
    )

    # Write the downsampled matrix to a new .mtx file
    mmwrite(DATA / "matrix.mtx.gz", mtx, compression=gzip)


def read_sample_and_write(indir=S1_1, outdir=Path("Data/default")):
    (barcodes, features, mtx) = read_data(dir)
    (barcodes_ds, features_ds, mtx_ds) = downsample(barcodes, features, mtx)
    write_data(barcodes_ds, features_ds, mtx_ds, dir)


if __name__ == "__main__":
    if len(argv) > 1:
        indata = argv[1]
    else:
        indata = S1_1

    if len(argv) > 2:
        outdata = Path(argv[2])
    else:
        outdata = DATA / f"outdata_{indata.name}"
    read_sample_and_write(indata, outdata)
