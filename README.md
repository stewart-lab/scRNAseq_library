markdown
Copy code
# scRNA-seq Docker Environment Setup Guide

This README provides instructions on how to set up a Docker environment for single-cell RNA sequencing (scRNA-seq) analysis.

## Prerequisites

1. [Docker](https://www.docker.com/products/docker-desktop) installed on your system.
   - Verify your installation by typing `docker --version` in your terminal. If Docker is not installed, you can download it from [here](https://www.docker.com/products/docker-desktop).
   
2. Git installed on your system.
   - Verify your installation by typing `git --version` in your terminal. If Git is not installed, you can download it from [here](https://git-scm.com/downloads).

## Instructions

Follow the steps below to set up your Docker environment:

1. **Clone the repository**

```bash
git clone git@github.com:stewart-lab/scRNAseq_library.git # Start by cloning the repository
```

2. **Navigate to the repository folder**
Change your current directory to the top level of the cloned repository using

```bash
cd <cloned_repo_dir>.
```

3. **Build the Docker image**
Run the following command:

```bash
docker build -t scrna-seq .
```

This will create a Docker image named scrna-seq.

4. **Run the Docker image**
Start an interactive terminal session within the Docker container:

```bash
docker run -it scrna-seq
```

5. **Retrieve data (if necessary)**
If you need data for your analysis, you can run the get_data.py script with one of the following --data arguments: REH, GAMM_S1, or GAMM_S2.

```bash
python get_data.py --data GAMM_S1 # Example for getting GAMM_S1 data
```

6. **Install R packages**
Next, run the script to install the necessary R packages:

```bash
source install_R_packages.sh
```

7. **Run the Analysis Script**
Finally, execute the R command to render your analysis results as a PDF:

```bash
R -e "rmarkdown::render('src/script.rmd', output_format = 'pdf_document')"
```
## References:

Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
