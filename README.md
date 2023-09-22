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
git clone https://github.com/stewart-lab/scRNAseq_library.git # Start by cloning the repository
```

2. **Navigate to the repository folder**
Change your current directory to the top level of the cloned repository using

```bash
cd scRNAseq_library/
```

3. **Run the Docker image**
Start an interactive terminal session within the Docker container and map the output directory of the docker filesystem to the repository:

```bash
source run_docker.sh
```

4. **Retrieve data (if necessary)**
If you need data for your analysis, you can run the get_data.py script with one of the following --data arguments: REH, GAMM_S1, or GAMM_S2.

```bash
python get_data.py --data GAMM_S1 # Example for getting GAMM_S1 data
```

5. **Turn on conda environment**
Next, run this conda command to activate the correct environment:

```bash
conda activate scrnaseq
```

6. **Run the Analysis Script**
Finally, execute the R command to render your analysis results as a PDF:

```bash
Rscript script.R
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
  - `analysis_type`: The type of analysis to use: "Scry" or "Seurat"

 - `scale_data`: A dictionary for specifying scale data parameters
   - `vars.2.regress`: genes to regress out- currently either "cell.cycle" or "NA"

- `run_and_visualize_pca`: A dictionary for specifying PCA parameters.
  - `top_n_dims`: The top n dimensions for PCA.
  - `heatmap_dims`: The number of dimensions for the heatmap.
  - `num_cells`: The number of cells to use.
  - `dims`: The number of dimensions to use for jackstraw.
  - `num.replicate`: The number of replicates for jackstraw plot.

- `run_umap`: a dictionary for running umap.
   - `dims_umap`: The number of dimensions to use in umap reduction
   - `umap.method`: Method to run umap: "umap-learn" or "uwot"
   - `umap.red`: Reduction method to use: "pca" or "harmony" (or special cases)

- `perform_batch_correction`: A dictionary for specifying batch correction parameters.
  - `dims.use`: The number of dimensions to use.
  - `max_iter`: The maximum number of iterations.

- `perform_clustering`: A dictionary for specifying clustering parameters.
  - `reduction`: Type of reduction to use for clustering: "harmony", "umap", or "pca"
  - `resolution`: The resolution for clustering (higher for more clusters, lower for less clusters).
  - `algorithm`: The algorithm used for clustering: "leiden"
  - `dims_snn`: Number of dimensions to use for KNN graph

- `find_differentially_expressed_features`: A dictionary for specifying parameters to find differentially expressed features.
  - `min_pct`: The minimum percentage for filtering.
  - `logfc_threshold`: The threshold for log fold-change.
  - `top_n`: The top n features to select.

- `score_and_plot_markers`: A dictionary for specifying parameters for scoring and plotting markers.
  - `top_n_markers`: The top n markers to use.
  - `known_markers`: Whether to use known markers.
  - `known_markers_path`: The path to the known markers: "../known_marker_lists/Gamm_lab_Consolidated_markerList.txt"

Please adjust the parameters as per your requirements. For additional details on each of these parameters, refer to the Seurat and SoupX documentation.

# Automated annotation

## Clustifyr

Clustifyr can take either a marker list or a reference seurat object/ dataset to annotate clusters in query data.

The introduction vignette can be found here: https://rnabioco.github.io/clustifyr/articles/clustifyR.html

To run, modify script with your data and marker list

```
src/clustifyr.R
```
- R requirements
```
BiocManager::install("clustifyr")
library(clustifyr)
library(ggplot2)
library(cowplot)
library(Seurat)
```
## scPred

scPred uses a reference object/ dataset to predict annotations of clusters in query data based on similarity of cell expression to the reference. The default algorithm is SVM-radial, however many different models/algorithms can be applied from the caret package: https://topepo.github.io/caret/available-models.html. 

The introduction vignette for scPred can be found here: https://powellgenomicslab.github.io/scPred/articles/introduction.html

To run, your query and reference data must first be processed the same way. Cross-species predictions need an ortholog file:

```
src/preprocess_crossspecies.Rmd
```
- R requirements

```
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(BiocParallel)
library(DropletUtils)
library(cowplot)
library(harmony)
library(SoupX)
library(scDblFinder)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
```

Next run the scPred script with your preprocessed query and reference objects. scPred will divide the reference into training and testing objects, where the model is trained on the training set, and then applied to the testing set. Watch for cell types that don't predict well in the test set- this may mean the model for that cell type isn't good, and you can try a different one. After a final model is built (you can have different algorithms for each cell type if you want), then you can apply to the query data. To run, update with your reference and query objects, you can also specify different algorithms/models after first running your training data with the SVM-radial algorithm.

```
src/scPred_GAMM.Rmd
```

- R requirements
```
library(devtools)
devtools::install_github("powellgenomicslab/scPred")
library(scPred)
library(Seurat)
library(magrittr)
library(harmony)
library(rmarkdown)
library(jsonlite)
library(purrr)
library(scran)
library(patchwork)
library(dplyr)
library(reticulate)
```

## Seurat mapping

Seurat mapping also uses a reference object/ dataset to predict annotations of clusters based on similairty of cell expression to the reference. It does this by doing a canonical correlation analysis to find "anchor" cells between the reference and query, then annotated clusters based on these anchor cells. 

The Seurat mapping vignette can be found here: https://satijalab.org/seurat/articles/integration_mapping

Again, query and reference data need to be preprocessed the same way. First run preprocess_crossspecies.Rmd, and cross-species predictions need an ortholog file. 

```
src/preprocess_crossspecies.Rmd
```

Next run the seurat mapping script on your preprocessed data. Seurat will integrate the two datasets together to find anchors, then transfer the anontation data from reference to query. Seurat lets you view the query in either query space or reference space.

```
src/seurat_mapping_GAMM.Rmd
```

- R requirements
```
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(BiocParallel)
library(DropletUtils)
library(cowplot)
library(harmony)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
library(scPred)
```

# Cell type composition analysis

## sccomp

sccomp measures cell type composition across datasets to see if the numbers of certain cell type have changed due to developmetn stage, condition, species, etc. It uses a beta-binomial distribution for each cell type that then sums to one.

The tutorial for sccomp is here: https://github.com/stemangiola/sccomp

To run sccomp you need processed and labeled data sets. Since you are comparing cell type compositionality, you need to have annotated data with the same cell type marker list or same reference set. 

```
src/sccomp.R
```

- R requirements
```
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("sccomp")
library(sccomp)
library(Seurat)
library(tidyverse)
library(loo)
```

# References: 

Star solo paper: https://doi.org/10.1101/2021.05.05.442755

clustifyr paper: https://doi.org/10.12688/f1000research.22969.2

scPred paper: https://doi.org/10.1186/s13059-019-1862-5

Seurat paper: https://doi.org/10.1016/j.cell.2019.05.031

sccomp paper: https://doi.org/10.1073/pnas.2203828120

