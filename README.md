# scRNA-seq Docker Environment Setup Guide

This README provides instructions on how to set up a Docker environment for single-cell RNA sequencing (scRNA-seq) analysis.

## Prerequisites

1. [Docker](https://www.docker.com/products/docker-desktop) installed on your system.
   - Verify your installation by typing `docker --version` in your terminal. If Docker is not installed, you can download it from [here](https://www.docker.com/products/docker-desktop).
   
2. Python installed on your system.
   - Verify your installation by typing `python --version` in your terminal.

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

# scRNA-seq Analysis Configuration Guide

This README provides a brief description of the configuration file used in the single-cell RNA sequencing (scRNA-seq) analysis.

## Configuration Keys

### General Settings

- `title`: The title for the analysis.
- `annotation_reference`: Indicates whether annotation reference is used in the analysis. (e.g., "FALSE")
- `DE_method`: The method used for differential expression analysis. (e.g., "Scran")
- `species`: The species for the analysis. (e.g., "human")
- `lanes`: A list to specify the lanes used in the analysis. (e.g., [])
- `DATA_DIR`: The directory where the data is stored. (e.g., "/isiseqruns/jfreeman_tmp_home/scRNA_FASTQS/")

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
  - `num.replicate`: The number of replicates for jackstraw plot. (e.g., 100)

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


# References: 

Star solo paper: doi: https://doi.org/10.1101/2021.05.05.442755
