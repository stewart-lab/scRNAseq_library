import pandas as pd
from scipy.io import mmread, mmwrite
import random
import numpy as np
from pathlib import Path

DATA = Path("Data")
S1_1 = DATA / "S1_1"
SAMPLES = S1_1 / "barcodes.tsv"
GENES = S1_1 / "features.tsv"
MTX = S1_1 / "matrix.mtx"


def read_data():
    barcodes = pd.read_csv(SAMPLES, header=None)
    features = pd.read_csv(GENES, header=None)
    mtx = mmread(MTX).tocsc()
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


def write_data(barcodes, features, mtx):
    # Write the downsampled barcodes and features to new .tsv files
    barcodes.to_csv(DATA / "barcodes_ds.tsv", index=False, header=False)
    features.to_csv(DATA / "features_ds.tsv", index=False, header=False)

    # Write the downsampled matrix to a new .mtx file
    mmwrite(DATA / "matrix_ds.mtx", mtx)


if __name__ == "__main__":
    (barcodes, features, mtx) = read_data()
    (barcodes_ds, features_ds, mtx_ds) = downsample(barcodes, features, mtx)
    write_data(barcodes_ds, features_ds, mtx_ds)
