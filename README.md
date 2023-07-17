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
Start an interactive terminal session within the Docker container and map the output directory of the docker filesystem to the repository:

```bash
docker run -v /output:scRNA-seq/output -it scrna-seq
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
# scRNA-seq Analysis Configuration Guide

This README provides a brief description of the configuration file used in the single-cell RNA sequencing (scRNA-seq) analysis.

Here is the explanation for each key in the configuration file:

- `title`: The title for the analysis.
  
- `annotation_reference`: Indicates whether annotation reference is used in the analysis.
  
- `DE_method`: The method used for differential expression analysis. In this case, "Seurat" is used.

- `lanes`: A list to specify the lanes used in the analysis.

- `prep_seurat_and_soupX`: A dictionary to specify parameters for Seurat and SoupX preparation. 
  - `dims`: The number of dimensions.
  - `umap.method`: The method used for UMAP.
  
- `filter_cells`: A dictionary for specifying cell filtering parameters.
  - `lower.nFeature`: The lower limit of features.
  - `upper.nFeature`: The upper limit of features.
  - `max.percent.mt`: The maximum percentage of mitochondrial content.
  - `species`: The species for the analysis.
  
- `normalize_data`: A dictionary for specifying normalization parameters.
  - `min_size`: The minimum size.
  - `min_mean`: The minimum mean.
  - `feature`: The feature to normalize.

- `feature_selection`: A dictionary for specifying feature selection parameters.
  - `n_features`: The number of features to select.
  - `analysis_type`: The type of analysis to use.

- `run_and_visualize_pca`: A dictionary for specifying PCA parameters.
  - `top_n_dims`: The top n dimensions for PCA.
  - `heatmap_dims`: The number of dimensions for the heatmap.
  - `num_cells`: The number of cells to use.

- `perform_batch_correction`: A dictionary for specifying batch correction parameters.
  - `dims.use`: The number of dimensions to use.
  - `max_iter`: The maximum number of iterations.

- `perform_clustering`: A dictionary for specifying clustering parameters.
  - `save_obj_before_clustering`: Whether to save object before clustering.
  - `num.replicate`: The number of replicates.
  - `dims`: The number of dimensions to use.
  - `dims_umap`: The number of dimensions for UMAP.
  - `resolution`: The resolution for clustering.
  - `algorithm`: The algorithm used for clustering.
  - `umap.method`: The method used for UMAP.

- `find_differentially_expressed_features`: A dictionary for specifying parameters to find differentially expressed features.
  - `min_pct`: The minimum percentage for filtering.
  - `logfc_threshold`: The threshold for log fold-change.
  - `top_n`: The top n features to select.

- `score_and_plot_markers`: A dictionary for specifying parameters for scoring and plotting markers.
  - `top_n_markers`: The top n markers to use.
  - `known_markers`: Whether to use known markers.
  - `known_markers_path`: The path to the known markers.

Please adjust the parameters as per your requirements. For additional details on each of these parameters, refer to the Seurat and SoupX documentation.

# References: 

Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
