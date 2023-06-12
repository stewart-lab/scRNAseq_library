# scRNA-seq Conda Environment Setup Guide

This README provides instructions on how to create a Conda environment using the `scRNAseq_venv.yml` file. 

## Prerequisites

1. [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.

## Instructions

Follow the steps below to create your Conda environment:

1. **Clone the repository**

   Start by cloning the repository


2. **Navigate to the repository folder**

   Change your current directory to the <cloned_repo_dir>/install


3. **Creating the Conda environment**

Now run the following code
```bash
source install_venv.sh ENV_NAME [ENV_PATH]
```
where ENV_NAME is the name of your environment and ENV_PATH is the path to your (non-default?) conda environment directory (optional). If ENV_PATH is not specified, the environment will be created in the default conda environment directory.
```


# scRNAseq_library
Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
