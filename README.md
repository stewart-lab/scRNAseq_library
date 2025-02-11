# scRNA-seq Analysis Guide

This README provides instructions on how to do a general single-cell RNA sequencing (scRNA-seq) analysis via Docker, including cross-species analysis and automated annotation via marker lists.

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
- `DE_method`: The method used for differential expression analysis. "Seurat" or "Scran" (we recommend "Scran")
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
  - `isBarcodeFollowedbyReads`: true or false, Parameter to indicate cDNA reads on barcode file
  - `clip5pNbases`: If bases should be clipped, what is the range to be kept. Reads outside of this range are clipped. (e.g.,"0 91")

### Data Preparation and Analysis

- `prep_seurat_and_soupX`: A dictionary to specify parameters for Seurat and SoupX preparation.
  - `dims`: The number of dimensions. (e.g., 30)
  - `umap.method`: The method used for UMAP. (e.g., "umap-learn")
  - `tfidfMin`: Minimum value of tfidf to accept for a marker gene to estimate background contamination. (e.g., 1) A higher tf-idf value implies a more specific marker.
  - `min.cells`: Include features (genes) detected in at least this many cells. (e.g., 3)

- `filter_cells`: A dictionary for specifying cell filtering parameters.
  - `lower.nFeature`: The lower limit of features (genes). (e.g., 200)
  - `upper.nFeature`: The upper limit of features (genes). (e.g., 25000)
  - `max.percent.mt`: The maximum percentage of mitochondrial content. (e.g., 20)

- `ortholog_subset`:
  - `ortholog_file`: File containg list of orthologs for cross-species analysis. From Ensemble, should contain tab-delimited: ref.gene.stable.ID	ref.gene.name	query.gene.stable.ID	query.gene.name
  - `ref_species`: Name of reference species. Default is human.

- `normalize_data`: A dictionary for specifying normalization parameters. We use Scran normalization.
  - `min_size`: The minimum size. (e.g., 100)
  - `min_mean`: The minimum mean. (e.g., 0.1)
  - `feature`: The feature to normalize. (e.g., "ECHS1")

- `feature_selection`: A dictionary for specifying feature selection parameters.
  - `n_features`: The number of features to select. (e.g., 2000)
  - `analysis_type`: The type of analysis to use, we recommend Scry. ("Scry" or "Seurat")

- `scale_data`: A dictionary for specifying scale data parameters.
  - `vars.2.regress`: Genes to regress out. (e.g., "cell.cycle")
  - `marker.path.s`: Path to cell cycle S genes. Default in repo: "../cell_cycle_vignette/cell_cycle_orthologs_s.genes.txt""
  - `marker.path.g2m`: Path to cell cycle G2M genes. Default in repo: "../cell_cycle_vignette/cell_cycle_orthologs_g2m.genes.txt"

- `run_and_visualize_pca`: A dictionary for specifying PCA parameters.
  - `top_n_dims`: The top n dimensions of PCA to display. (e.g., 2)
  - `heatmap_dims`: The number of dimensions for the heatmap. (e.g., 15)
  - `num_cells`: The number of cells to use. (e.g., 500)
  - `dims`: The number of dimensions to use for jackstraw. (e.g., 20)
  - `num.replicate`: The number of replicates for jackstraw plot. (e.g., 100, or can be NA). If 'NA' jackstraw is not run.

- `run_umap`: A dictionary for running UMAP.
  - `dims_umap`: The number of dimensions to use in UMAP reduction. (e.g., 20)
  - `umap.method`: Method to run UMAP. (e.g., "umap-learn")
  - `umap.red`: Reduction method to use. (e.g., "pca" or "harmony")

- `perform_batch_correction`: A dictionary for specifying batch correction parameters.
  - `dims.use`: The number of dimensions to use. (e.g., 20)
  - `max_iter`: The maximum number of iterations. (e.g., 50)

- `perform_clustering`: A dictionary for specifying clustering parameters.
  - `reduction`: Type of reduction to use for KNN graph. (e.g., "pca" or "harmony")
  - `resolution`: The resolution for clustering. Lower means fewer clusters, higher means more clusters. (e.g., 0.5)
  - `algorithm`: The algorithm used for clustering. (e.g., "leiden")
  - `dims_snn`: Number of dimensions to use for KNN graph. (e.g., 10)

- `find_differentially_expressed_features`: A dictionary for specifying parameters to find differentially expressed features. Used if using Seurat to identify DE, if using Scran then score_and_plot_markers is used.
  - `min_pct`: The minimum percentage for filtering. (e.g., 0.25)
  - `logfc_threshold`: The threshold for log fold-change. (e.g., 0.25)
  - `top_n`: The top n features to select. (e.g., 11)

- `get_metadata`: Adding in known metadata i.e. if a reference is being used and this reference has already been annotated.
  - `metadata_file_ref`: Reference metadata file. Should be tab-delimited with cells as first column and subsequent cell information in the following columns. Example files in metadata folder. Put file in metadata folder, then "../metadata/filename.txt".
  - `metadata_file_query`: Query metadata file. Same format as above.
  - `metadata_subset1`: Subset metadata for reference. Subsets for a name in "source" column. If `NA`, full metadata is used.
  - `metadata_subset2`: Subset metadata for query. Same format as above.

- `score_and_plot_markers`: A dictionary for specifying parameters for scoring and plotting markers. DE genes scored and found by Scran.
  - `top_n_markers`: The top n markers to use (cut off for how many markers to find). (e.g., `100`)
  - `known_markers`: Whether to use known markers. If FALSE, manual annotation cannot be done and only returns DE gene list. (`TRUE` or `FALSE`)
  - `known_markers_path`: The path to the known markers. (e.g., `../known_marker_lists/Gamm_lab_Consolidated_markerList.txt`)
  - `cluster_type`: Cluster type to determine DE genes/ markers. (e.g. `seurat_clusters`,`orig.ident`)
  - `pairwise`: Do you want to calculate all pairwise comparisons between clusters? (`TRUE` or `FALSE`)
  - `logFC_thresh`: Cohen's D log fold change threshold, only DE genes above this threshold are kept. (e.g. `0.25`)
  - `auc_thresh`: Area-under-the-curve threshold (the probability that a randomly chosen observation from one group is greater than a random). (e.g. `0.49`)
 
- `process_known_markers`: A dictonary to determine how to annotate clusters with known markers
  - `annot_type`: Type of annotation, manual being using markers in the top n_rank to annotate. Other options are related to Gamm paper. (`manual`,`d40`,`d120`)
  - `n_rank`: Lowest rank based on the log fold change to consider when annotating cell type (e.g. `10`)

Please adjust the parameters as per your requirements. For additional details on each of these parameters, refer to the Seurat, Scran, SoupX documentation.
Note: **Clustifyr** used to do automated annotation via provided marker list in addition to manual annotation.

### Output files

- SoupX
   - `post_soupx_qc_combined.pdf`: nCount by nFeature after SoupX
- scDblFinder
   - `merged_doublet_table.txt`: Number of doublets and singlets called
   - `after_dbl_removal_and_merge.pdf`: nCount by nFeature after doublet removal
- Mitochondrial filtering
   - `percent_mt_unfiltered.pdf`: percent.mt by nCount before mitochondrial filtering
   - `percent_mt_filtered.pdf`: percent.mt by nCount after mitochondrial filtering
- Normalization
   - `violin_pre_norm.pdf`: Expression levels of a gene before normalization
   - `violin_post_norm.pdf`: Expression levels of a gene after normalization
- Scaling
   - `pca_before_cc_regression.pdf`: PCA before cell cycle regression
   - `pca_after_cc_regression.pdf`: PCA after cell cycle regression
- PCA
   - `pca_heat_map.pdf`: Heatmap of the top n PCs to look for gene variability
   - `jack_straw.pdf`: Top n PCs variability plot with p-values
   - `top_n_dims_with_genes.pdf`: PC1 and PC2 gene variablity
   - `elbow_pca.pdf`: standard deviation by PC plot
- Batch correction
   - `batch_uncorrected_pca.pdf`: PCA before batch correction.
   - `batch_corrected_pca.pdf`: PCA after batch correction.
- Umap plots
   - `umap_plot.pdf`: Umap with cell cycle genes
   - `umap_lanes.pdf`: Umap colored by sample
   - `umap_clusters.pdf`: Umap of clusters
   - `labeled-clusters.pdf`: Umap of labeled clusters (if marker list is used)
   - `clustifyr_marker_annotation_umap.pdf`: Umap of clusters annotated by Clustifyr (if marker list is used)
- DE gene/ Marker files:
   - See Scran's `scoreMarkers` for details on columns in output files: https://rdrr.io/github/MarioniLab/scran/man/scoreMarkers.html
   - `Top100genes_clust` files: contain the top 100 genes for a particular cluster against all other clusters, sorted by median.logFC.cohen and subsetted by the median Cohen's D logFC threshold and median AUC threshold.
   - If a known markers list is given, this list is merged with the `Top100genes_clust` file. These are the `KnownDE.markers_clust_` files.
   - The top ranked markers from the merged list are used to make feature umap plots highlighting the gene expression of that marker. Default is the ranked in the top 10.
   - Pairwise comparisons are also made if `pairwise=TRUE`. Each pairwise comparison between every cluster is made, and DE genes kept with the genes are higher than the median Cohen's D logFC threshold and median AUC threshold. These files start with `DEgenes_clust_X.vs_Y`.
- `sc_pipeline.pdf`: Code and output from the sc_pipeline
- Seurat objects:
   - `seurat_obj_labeled.rds` contains manual annotation (if manual annotation was selected, otherwise contains just clusters)
   - `seurat_obj_clustifyr.rds` contains clustifyr annotation if marker files were provided.


# References: 

Star solo paper: https://doi.org/10.1101/2021.05.05.442755

SoupX paper: https://doi.org/10.1093/gigascience/giaa151

scDblFinder paper: https://doi.org/10.12688/f1000research.73600.2

Seurat paper: https://doi.org/10.1016/j.cell.2019.05.031

Scran paper: https://doi.org/10.1186/s13059-016-0947-7

Clustifyr paper: https://doi.org/10.12688/f1000research.22969.2

