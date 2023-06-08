# scRNA-seq Conda Environment Setup Guide

This README provides instructions on how to create a Conda environment using the `scRNA-seq.yml` file. 

## Prerequisites

1. [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.

## Instructions

Follow the steps below to create your Conda environment:

1. **Clone the repository**

   Start by cloning the repository which contains the `scRNA-seq.yml` file. This can be done with the following command:




Replace `<repository_url>` with the URL of the repository.

2. **Navigate to the repository folder**

Change your current directory to the cloned repository:




Replace `<repository_folder>` with the name of the directory that was created when you cloned the repository.

3. **Create the Conda environment**

Use the following command to create the Conda environment:


This command creates a new Conda environment using the configuration specified in the `scRNA-seq.yml` file.

**IMPORTANT:** Before running the command, open the `scRNA-seq.yml` file in a text editor and change the `prefix:` line to match the location where your Conda is installed. This location is usually the `envs` directory inside your Anaconda or Miniconda installation directory. For example, if Anaconda is installed at `/home/user/anaconda3`, the prefix should be `/home/user/anaconda3/envs`.

That's it! You now have a Conda environment set up with the configuration specified in the `scRNA-seq.yml` file. You can start using it for your single-cell RNA sequencing analysis.





# scRNAseq_library
Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
