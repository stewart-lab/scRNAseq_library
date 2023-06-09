# scRNA-seq Conda Environment Setup Guide

This README provides instructions on how to create a Conda environment using the `scRNAseq_venv.yml` file. 

## Prerequisites

1. [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.

## Instructions

Follow the steps below to create your Conda environment:

1. **Clone the repository**

   Start by cloning the repository which contains the `scRNAseq_venv.yml` file. 


2. **Navigate to the repository folder**

   Change your current directory to the cloned repository


3. **Create the Conda environment**

   Use the following command to create the Conda environment:
   
```bash
conda env create -f scRNAseq_venv.yml
```


   This command creates a new Conda environment using the configuration specified in the `scRNA-seq.yml` file.

   **IMPORTANT:** Before running the command, open the `scRNA-seq.yml` file in a text editor and change the `prefix:` line (LOCATED AY BOTTOM) to match the location where your Conda is installed. This      location is usually the `envs` directory inside your Anaconda or Miniconda installation directory. For example, if Anaconda is installed at `/home/user/anaconda3`, the prefix          should be `/home/user/anaconda3/envs`.

4. **Activate R within the Conda environment**
```bash
conda activate scRNA-seq
R
```

That's it! You now have a Conda environment set up with the configuration specified in the `scRNAseq_venv.yml` file. You can start using it for your single-cell RNA sequencing analysis.


# scRNAseq_library
Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
