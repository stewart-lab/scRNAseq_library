# scRNA-seq Analysis Guide

This README provides instructions on how to do a general single-cell RNA sequencing (scRNA-seq) analysis via Docker, including cross-species analysis and automated annotation via marker lists.

## Prerequisites

1. [Docker](https://www.docker.com/products/docker-desktop) installed on your system.
   - Verify your installation by typing `docker --version` in your terminal. If Docker is not installed, you can download it from [here](https://www.docker.com/products/docker-desktop).
   
2. [Python](https://docs.python.org/3/using/index.html) installed on your system.
   - Verify your installation by typing `python --version` in your terminal.If Python is not installed, you can download it [here](https://docs.python.org/3/using/index.html)

**Instead of Docker, you can use Apptainer**

3. [Apptainer](https://apptainer.org/docs/user/main/quick_start.html) installed on your system.
   - verify by typing `apptainer --version` in your terminal.
     
4. If using apptainer, you may want [tmux](https://github.com/tmux/tmux/wiki) in order to run in the background. If you have tmux, when 

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


4. **Run the Pipeline on Docker**
This command will intialize the pipeline. You will be asked questions based upon your data and analysis requirements.:

```bash
# if you want to run the pipeline together with both alignment and single cell analysis. Note this way runs the aligner part interactively, so for large jobs it may be better to run in detached mode- see below
source run_pipeline.sh

# or to run separately run aligner first (this is in detached mode)
source run_aligner.sh

# when aligner finishes then run seurat in detached mode
source run_seuratv5.sh
```

5. **Analysis questions**
   - Would you like to run alignment? [y/N]:
     - If yes, the previous alignment will be deleted and pipeline will look to the specified fastq files in the DATA_DIR to realign with STARsolo
     - If no, the seurat program will look for the star solo output in the shared_mount directory. Edits to the config file will be applied and a new time-stamped output will be generated.

## Running Apptainer
1. **Pull docker image with apptainer** First after apptainer is installed, pull the docker image using apptainer. This will convert the docker image to an apptainer containter
```bash
apptainer pull docker://stewartlab/sc_aligner_v2_no_genomes
```

This creates a .sif file that is your apptainer

2. **Configuring your pipeline run**
Before running the pipeline, make sure to configure your settings in the config.json file. For more details on how to set up the configuration, see the [scRNA-seq Analysis Configuration Guide](#scrna-seq-analysis-configuration-guide) below.

3. **Run Apptainer**
```bash
# run aligner
source run_apptainer_aligner.sh 
```
4. **analysis questions**
   - Would you like to run alignment? [y/N]:
     - If yes, the previous alignment will be deleted and pipeline will look to the specified fastq files in the DATA_DIR to realign with STARsolo
     - If no, the seurat program will look for the star solo output in the shared_mount directory. Edits to the config file will be applied and a new time-stamped output will be generated.
       
   - Do you need to build a genome index? [y/N]:
      - If yes, you must provide the genome fasta file and genome gtf file as well as the genome directory where they are located in the config file.
      - If no, this assumes you have already built the genome index, and the program will look for the genome index in the "GENOME_INDEX_DIR". This should be a subdirectory from the genome directory ("GENOME_DIR") in the config file.

   - Run in detached tmux session? [y/N]:
      - If you have tmux installed, you can run in detached mode which means you do not have to keep the terminal open while alignment is running.

# scRNA-seq Analysis Configuration Guide

This README provides a brief description of the configuration file used in the single-cell RNA sequencing (scRNA-seq) analysis.

## Configuration Keys

### General Settings

- `title`: The title for the analysis.
- `annotation_reference`: Indicates whether annotation reference is used in cluster labeling. (e.g., "FALSE")
- `DE_method`: The method used for differential expression analysis. "Seurat" or "Scran" (we recommend "Scran")
- `species`: The species for the analysis. (e.g., "human", "pig")
- `Number_of_samples`: How many different samples you are running (not different lanes)
- `lanes`: THIS WILL BE AUTOPOPULATED VIA THE PIPELINE
- `DATA_DIR`: The directory where the FASTQs are stored. (e.g., "/isiseqruns/jfreeman_tmp_home/scRNA_FASTQS/")
- `GENOME_DIR`: The directory where your genome files are stored. (e.g., "/w5home/bmoore/genomes/human_genome_38/")
- `GENOME_FASTA`: The name of the primary genome fasta file (this should be in the GENOME_DIR, e.g., "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")
- `GENOME_GTF`: The name of the gtf file (this should be in the GENOME_DIR, e.g.,"Homo_sapiens.GRCh38.108.filtered.gtf")
- `GENOME_INDEX_DIR`: The name of the subdirectory of the GENOME_DIR where the genome index is. If building the index, leave as the default: "genome_build/"
- `RUN_THREAD_N`: Number of cpus to use (e.g. 16)

### Fastq Alignment Settings

- `fastq_alignment`: A dictionary for specifying fastq alignment parameters. Each item is a separate dictionary for each sample
   - "SAMPLE_1": dictionary for the first sample settings. For more than one sample, copy the sample_1 dictionary, change to SAMPLE_2, and edit settings. You can have as many samples as your docker/apptainer memory and cpus permit.
  
  Within each sample, edit parameters:
  - `NAME`: Name of sample
  - `NUM_LANES`: Number of lanes. (e.g., 1)
  - `cDNA_LANE1:`: Name of fastq read file with cDNA. (e.g., "/data/AtlasOfTheHumanRetina/SRR10130821_R2.fastq.gz") # to add more lanes, simply add another cDNA and barcode lane, i.e. cDNA_LANE2 and BARCODE_LANE2 with the corresponding fastq files
  - `BARCODE_LANE1`: Name of fastq read file with barcode. (e.g., "/data/AtlasOfTheHumanRetina/SRR10130821_R1.fastq.gz")
  - `CHEMISTRY_VERSION`: Chemistry version. (Versions 2,3, and 4 are supported: e.g., "V3")
  - `SOLO_TYPE`: Solo type. (e.g., "CB_UMI_Simple")
  - `SOLO_FEATURES`: Solo features. (e.g., "Gene GeneFull SJ Velocyto")
  - `SOLO_CELL_FILTER`: Solo cell filter. (e.g., "EmptyDrops_CR")
  - `SOLO_MULTI_MAPPERS`: Solo multi-mappers. (e.g., "EM")
  - `READ_FILES_COMMAND`: Read files command. (if reads are gzipped: "zcat")
  - `SOLO_UMI_DEDUP`: Solo UMI deduplication. (e.g., "1MM_CR")
  - `RUN_THREAD_N`: Number of threads to run. (e.g., 8)
  - `isBarcodeFollowedbyReads`: true or false, Parameter to indicate cDNA reads on barcode file
  - `clip5pNbases`: If bases should be clipped, what is the range to be kept. Reads outside of this range are clipped. (e.g.,"0 91")
  - `soloStrand`: Which direction, Forward or Reverse, default is: "Forward"

### Data Preparation and Analysis

- `prep_seurat_and_soupX`: A dictionary to specify parameters for Seurat and SoupX preparation.
  - `dims`: The number of dimensions. (e.g., 30)
  - `umap.method`: The method used for UMAP. (e.g., "umap-learn")
  - `tfidfMin`: Minimum value of tfidf to accept for a marker gene to estimate background contamination. (e.g., 1) A higher tf-idf value implies a more specific marker.
  - `min.cells`: Include features (genes) detected in at least this many cells. (e.g., 3)

- `filter_cells`: A dictionary for specifying cell filtering parameters.
  - `lower.nFeature`: The lower limit of features (genes). (e.g., 200)
  - `upper.nFeature`: The upper limit of features (genes). (e.g., 25000)
  - `max.percent.mt`: The maximum percentage of mitochondrial content allowed. (over this percentage, cells are filtered out- e.g., 20)

- `ortholog_subset`:
  - `ortholog_file`: File containg list of orthologs for cross-species analysis. From Ensemble, should contain tab-delimited: ref.gene.stable.ID	ref.gene.name	query.gene.stable.ID	query.gene.name
  - `ref_species`: Name of reference species. Default is human.

- `normalize_data`: A dictionary for specifying normalization parameters. We use Scran normalization.
  - `min_size`: The minimum size. (e.g., 100)
  - `min_mean`: The minimum mean. (e.g., 0.1)
  - `feature`: The name of the feature (gene name) to show normalization in a figure. (e.g., "ECHS1")

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




