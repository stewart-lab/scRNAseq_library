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


3. **Creating the Conda environment**

   Click on the install_venv.sh file and set the "conda_envs_path" where your conda environment directory is located (eg. ("/home/jfreeman/bin/miniconda3/envs"))
   
   Additionally, set your name of your environment
   ![Screenshot 2023-06-12 at 1 09 29 PM](https://github.com/stewart-lab/scRNAseq_library/assets/95723801/c6e4f7d3-7ebb-4ceb-8481-beee3bc2b33b)

Now run the following code
```bash
source install_venv.sh 
```

That's it! You now have a Conda environment set up with the configuration specified in the `environment_template.yml` file. You can start using it for your single-cell RNA sequencing analysis.


# scRNAseq_library
Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
