# scRNA-seq Conda Environment Setup Guide

This README provides instructions on how to create a Conda environment using the `scRNAseq_venv.yml` file. 

## Prerequisites

1. [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.
   1. Check with `conda --version` in your terminal. https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html 
2. VSCode is higly recommended as well. You can download it [here](https://code.visualstudio.com/download).
   1. The R extension needs to be installed in VSCode. https://marketplace.visualstudio.com/items?itemName=Ikuyadeu.r

## Instructions

Follow the steps below to create your Conda environment:

1. **Clone the repository**

```bash
git clone git@github.com:stewart-lab/scRNAseq_library.git # Start by cloning the repository
```


2. **Navigate to the repository folder**

   Change your current directory to the <cloned_repo_dir>/install


3. **Creating the Conda environment**

Now run the following code
```bash
source install_venv.sh ENV_NAME [ENV_PATH]
```
where ENV_NAME is the name of your environment and ENV_PATH is the path to your (non-default?) conda environment directory (optional). If ENV_PATH is not specified, the environment will be created in the default conda environment directory.

4. Edit the script.rmd file to include your data and run the script. 

## References:

Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
