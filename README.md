# scRNA-seq Analysis Guide

This README provides instructions on how to do a general single-cell RNA sequencing (scRNA-seq) analysis via Docker, as well as 
some specialized analyses such as automated annotation, cell type composition, and pseudotime.

## Prerequisites

1. [Docker](https://www.docker.com/products/docker-desktop) installed on your system.
   - Verify your installation by typing `docker --version` in your terminal. If Docker is not installed, you can download it from [here](https://www.docker.com/products/docker-desktop).
   
2. [Python](https://docs.python.org/3/using/index.html) installed on your system.
   - Verify your installation by typing `python --version` in your terminal.If Python is not installed, you can download it [here](https://docs.python.org/3/using/index.html)


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
3. **Configuring your pipeline run**
Before running the pipeline, make sure to configure your settings in the config.json file. For more details on how to set up the configuration, see the [scRNA-seq Analysis Configuration Guide](#scrna-seq-analysis-configuration-guide) below.


4. **Run the Pipeline**
This command will intialize the pipeline. You will be asked questions based upon your data and analysis requirements.:

```bash
source run_pipeline.sh 
```
5. **Analysis questions**
   - Have you loaded new data or would you like to realign? [y/N]:
     - If yes, the previous alignment will be deleted and pipeline will look to the specified files in the DATA_DIR to realign with STARsolo
     - If no, edits to the config file that are pipeline specfic (clustering, MT filtering, scaling) will be updated and a new time-stamped output will be generated
     - If no the following follow-up question will be asked: If you'd like to load a stored experiment select data. If you have aligned FASTQs loaded and changed pipeline parameters, select fastq [data/fastq]:
          - If data, that means you would like to load one of our pre-aligned datasets and you must select between the three: [REH,GAMM_S1,GAMM_S2]
          - If fastq, the alignment step will be skipped but your presumabley new config parameters will be applied to the latest time-stamped run 
    
    
   

# scRNA-seq Analysis Configuration Guide

This README provides a brief description of the configuration file used in the single-cell RNA sequencing (scRNA-seq) analysis.

## Configuration Keys

### General Settings

- `title`: The title for the analysis.
- `annotation_reference`: Indicates whether annotation reference is used in cluster labeling. (e.g., "FALSE")
- `DE_method`: The method used for differential expression analysis. (we recommend "Scran")
- `species`: The species for the analysis. (e.g., "human", "pig")
- `lanes`: THIS WILL BE AUTOPOPULATED VIA THE PIPELINE
- `DATA_DIR`: The directory where the FASTQs are stored. (e.g., "/isiseqruns/jfreeman_tmp_home/scRNA_FASTQS/")

### Fastq Alignment Settings

- `fastq_alignment`: A dictionary for specifying fastq alignment parameters.
  - `NUM_LANES`: Number of lanes. (e.g., 1)
  - `OUTPUT_PREFIX`: Output prefix. (e.g., "TEST")
  - `READ_FILE1_PREFIX`: Prefix for read file 1. (e.g., "/data/AtlasOfTheHumanRetina/SRR10130821_R2.fastq.gz")
  - `READ_FILE2_PREFIX`: Prefix for read file 2. (e.g., "/data/AtlasOfTheHumanRetina/SRR10130821_R1.fastq.gz")
  - `CHEMISTRY_VERSION`: Chemistry version. (e.g., "V3")
  - `SOLO_TYPE`: Solo type. (e.g., "CB_UMI_Simple")
  - `SOLO_FEATURES`: Solo features. (e.g., "Gene GeneFull SJ Velocyto")
  - `SOLO_CELL_FILTER`: Solo cell filter. (e.g., "EmptyDrops_CR")
  - `SOLO_MULTI_MAPPERS`: Solo multi-mappers. (e.g., "EM")
  - `READ_FILES_COMMAND`: Read files command. (e.g., "zcat")
  - `SOLO_UMI_DEDUP`: Solo UMI deduplication. (e.g., "1MM_CR")
  - `RUN_THREAD_N`: Number of threads to run. (e.g., 8)

### Data Preparation and Analysis

- `prep_seurat_and_soupX`: A dictionary to specify parameters for Seurat and SoupX preparation.
  - `dims`: The number of dimensions. (e.g., 30)
  - `umap.method`: The method used for UMAP. (e.g., "umap-learn")

- `filter_cells`: A dictionary for specifying cell filtering parameters.
  - `lower.nFeature`: The lower limit of features. (e.g., 200)
  - `upper.nFeature`: The upper limit of features. (e.g., 10000)
  - `max.percent.mt`: The maximum percentage of mitochondrial content. (e.g., 20)

- `normalize_data`: A dictionary for specifying normalization parameters.
  - `min_size`: The minimum size. (e.g., 100)
  - `min_mean`: The minimum mean. (e.g., 0.1)
  - `feature`: The feature to normalize. (e.g., "ECHS1")

- `feature_selection`: A dictionary for specifying feature selection parameters.
  - `n_features`: The number of features to select. (e.g., 2000)
  - `analysis_type`: The type of analysis to use. (e.g., "Scry")

- `scale_data`: A dictionary for specifying scale data parameters.
  - `vars.2.regress`: Genes to regress out. (e.g., "cell.cycle")

- `run_and_visualize_pca`: A dictionary for specifying PCA parameters.
  - `top_n_dims`: The top n dimensions for PCA. (e.g., 2)
  - `heatmap_dims`: The number of dimensions for the heatmap. (e.g., 15)
  - `num_cells`: The number of cells to use. (e.g., 500)
  - `dims`: The number of dimensions to use for jackstraw. (e.g., 20)
  - `num.replicate`: The number of replicates for jackstraw plot. (e.g., 100, or can be NA). If 'NA' jackstraw is not run.

- `run_umap`: A dictionary for running UMAP.
  - `dims_umap`: The number of dimensions to use in UMAP reduction. (e.g., 20)
  - `umap.method`: Method to run UMAP. (e.g., "umap-learn")
  - `umap.red`: Reduction method to use. (e.g., "pca")

- `perform_batch_correction`: A dictionary for specifying batch correction parameters.
  - `dims.use`: The number of dimensions to use. (e.g., 20)
  - `max_iter`: The maximum number of iterations. (e.g., 50)

- `perform_clustering`: A dictionary for specifying clustering parameters.
  - `reduction`: Type of reduction to use for clustering. (e.g., "harmony")
  - `resolution`: The resolution for clustering. (e.g., 0.5)
  - `algorithm`: The algorithm used for clustering. (e.g., "leiden")
  - `dims_snn`: Number of dimensions to use for KNN graph. (e.g., 10)

- `find_differentially_expressed_features`: A dictionary for specifying parameters to find differentially expressed features.
  - `min_pct`: The minimum percentage for filtering. (e.g., 0.25)
  - `logfc_threshold`: The threshold for log fold-change. (e.g., 0.25)
  - `top_n`: The top n features to select. (e.g., 11)

- `score_and_plot_markers`: A dictionary for specifying parameters for scoring and plotting markers.
  - `top_n_markers`: The top n markers to use. (e.g., 10)
  - `known_markers`: Whether to use known markers. (e.g., "TRUE")
  - `known_markers_path`: The path to the known markers. (e.g., "../known_marker_lists/Gamm_lab_Consolidated_markerList.txt")

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
# Pseudotime/ trajectory analysis

Cells are often in transition from one cell type to another, and pseudotime captures relationships between clusters, beginning with the least differentiated state to the most mature/terminal state(s).

## Palantir and CellRank

Palantir models trajectories of differentiated cells by treating cell fate as a probabilistic process and leverages entropy to measure cell plasticity along the trajectory. CellRank uses Palantir in
it's pseudotime kernel, but can also use RNA velocity, similarity, cytotrace, time-series, or matebolic labeling to calculate trajectories. Here we use it with Palantir. Together they identify initial
and terminal states, cell fate probabilities, and driver genes.

Tutorials:
- Palantir: https://github.com/dpeerlab/Palantir
- CellRank: https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/general/100_getting_started.html

To run Palantir and CellRank you first have to have a Seurat object that is clustered and annotated (see above). 

Next convert your Seurat object to an anndata object:

```
# before running, change the working directory and the input and output filenames in the script.
src/convert_seurat2anndata.R
```

- R requirements
```
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
library(SeuratDisk)
```

Once you have your anndata object, set up your python environment:

```
# set up with pseudotime_requirements.txt
source src/install_pseudotime_env.sh
# activate
source pst_venv/bin/activate
```

Now you are ready to run Palantir and CellRank

- open src/pseudotime_GAMM.ipynb is vscode and update the following variables:

```
# variables
DATA_DIR = "your_directory" # directory where anndata object is
ADATA_FILE = "gamms2_cca_pred.h5ad" # the name of your h5ad (anndata) file
ANNOT_TYPE = "manual" # the type of annotation used, options are: "seurat_map", "clustifyr", "manual"
CROSS_SPECIES = "TRUE" # is this a cross-species annotation? "TRUE" or "FALSE"
NC = 8 # number of components that are used to find terminal cells. In general, lower for few terminal cell types, higher for many terminal cell types
```

- now run src/pseudotime_GAMM.ipynb
- figures will be saved to the data directory

# References: 

Star solo paper: https://doi.org/10.1101/2021.05.05.442755

clustifyr paper: https://doi.org/10.12688/f1000research.22969.2

scPred paper: https://doi.org/10.1186/s13059-019-1862-5

Seurat paper: https://doi.org/10.1016/j.cell.2019.05.031

sccomp paper: https://doi.org/10.1073/pnas.2203828120

Palantir paper: https://doi.org/10.1038/s41587-019-0068-4

CellRank2 paper: https://doi.org/10.1101/2023.07.19.549685

