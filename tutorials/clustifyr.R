install.packages("BiocManager")

BiocManager::install("clustifyr")

library(clustifyr)
library(ggplot2)
library(cowplot)

# Matrix of normalized single-cell RNA-seq counts
pbmc_matrix <- clustifyr::pbmc_matrix_small

# meta.data table containing cluster assignments for each cell
# The table that we are using also contains the known cell identities in the "classified" column
pbmc_meta <- clustifyr::pbmc_meta

# input: an SingleCellExperiment or Seurat object or a matrix of normalized single-cell RNA-seq counts
# metadata: a meta.data table containing the cluster assignments for each cell (not required if a Seurat object is given)
# ref_mat: a reference matrix containing RNA-seq expression data for each cell type of interest
# query_genes: a list of genes to use for comparison (optional but recommended)

# Calculate correlation coefficients for each cluster (spearman by default)
vargenes <- pbmc_vargenes[1:500]

res <- clustify(
  input = pbmc_matrix, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata = pbmc_meta, # meta.data table containing cell clusters
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  ref_mat = cbmc_ref, # matrix of RNA-seq expression data for each cell type
  query_genes = vargenes # list of highly varible genes identified with Seurat
)

# Peek at correlation matrix
res[1:5, 1:5]

# Call cell types
res2 <- cor_to_call(
  cor_mat = res,                  # matrix correlation coefficients
  cluster_col = "seurat_clusters" # name of column in meta.data containing cell clusters
)
res2[1:5, ]

# Insert into original metadata as "type" column
pbmc_meta2 <- call_to_metadata(
  res = res2,                     # data.frame of called cell type for each cluster
  metadata = pbmc_meta,           # original meta.data table containing cell clusters
  cluster_col = "seurat_clusters" # name of column in meta.data containing cell clusters
)

# plot_cor_heatmap() function to plot the correlation coefficients for each cluster and each cell type.
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# Create heatmap of correlation coefficients using clustifyr() output
library(circlize)
library("viridis")
viridis(6)
col_fun = colorRamp2(seq(0, 1, length = 6), rev(viridis(6)))
plot_cor_heatmap(cor_mat = res, col=col_fun)

# plot cluster identities and corelation coefficients
# Overlay correlation coefficients on UMAPs for the first two cell types
corr_umaps <- plot_cor(
  cor_mat = res,                     # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta,              # meta.data table containing UMAP or tSNE data
  data_to_plot = colnames(res)[1:2], # name of cell type(s) to plot correlation coefficients
  cluster_col = "seurat_clusters"    # name of column in meta.data containing cell clusters
)

plot_grid(
  plotlist = corr_umaps,
  rel_widths = c(0.47, 0.53)
)
# The plot_best_call() function can be used to label each cluster with the cell 
# type that gives the highest corelation coefficient. 
# Label clusters with clustifyr cell identities
clustifyr_types <- plot_best_call(
  cor_mat = res,          # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta,   # meta.data table containing UMAP or tSNE data
  do_label = TRUE,        # should the feature label be shown on each cluster?
  do_legend = TRUE,      # should the legend be shown?
  cluster_col = "seurat_clusters"
) +
  ggtitle("clustifyr cell types")

# Compare clustifyr results with known cell identities
known_types <- plot_dims(
  data = pbmc_meta,       # meta.data table containing UMAP or tSNE data
  feature = "classified", # name of column in meta.data to color clusters by
  do_label = TRUE,        # should the feature label be shown on each cluster?
  do_legend = TRUE       # should the legend be shown?
) +
  ggtitle("Known cell types")

plot_grid(known_types, clustifyr_types)

# Classify cells using known marker genes
# The clustify_lists() function allows cell types to be assigned based on known marker genes.
# Cell types can be assigned using several statistical tests including, hypergeometric, Jaccard, Spearman, and GSEA.

# Take a peek at marker gene table
cbmc_m

# Available metrics include: "hyper", "jaccard", "spearman", "gsea"
list_res <- clustify_lists(
  input = pbmc_matrix,             # matrix of normalized single-cell RNA-seq counts
  metadata = pbmc_meta,            # meta.data table containing cell clusters
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  marker = cbmc_m,                 # list of known marker genes
  metric = "pct"                   # test to use for assigning cell types
)

# View as heatmap, or plot_best_call
plot_cor_heatmap(
  cor_mat = list_res,              # matrix of correlation coefficients from clustify_lists()
  cluster_rows = FALSE,            # cluster by row?
  cluster_columns = FALSE,         # cluster by column?
  legend_title = "% expressed",     # title of heatmap legend
  col=col_fun
  )


# Downstream functions same as clustify()
# Call cell types
list_res2 <- cor_to_call(
  cor_mat = list_res,              # matrix correlation coefficients
  cluster_col = "seurat_clusters"  # name of column in meta.data containing cell clusters
)

# Insert into original metadata as "list_type" column
pbmc_meta3 <- call_to_metadata(
  res = list_res2,                 # data.frame of called cell type for each cluster
  metadata = pbmc_meta,            # original meta.data table containing cell clusters
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  rename_prefix = "list_"          # set a prefix for the new column
)

# Direct handling of SingleCellExperiment objects
res <- clustify(
  input = sce_small,          # an SCE object
  ref_mat = cbmc_ref,         # matrix of RNA-seq expression data for each cell type
  cluster_col = "cell_type1", # name of column in meta.data containing cell clusters
  obj_out = TRUE              # output SCE object with cell type inserted as "type" column
)

SingleCellExperiment::colData(res)[1:10,c("type", "r")]

# Direct handling of seurat v2 and v3 objects

res <- clustify(
  input = s_small3,       # a Seurat object
  ref_mat = cbmc_ref,    # matrix of RNA-seq expression data for each cell type
  cluster_col = "RNA_snn_res.1", # name of column in meta.data containing cell clusters
  obj_out = TRUE      # output Seurat object with cell type inserted as "type" column
)

res@meta.data[1:10, ]

# Building reference matrix from single cell expression matrix
# n its simplest form, a reference matrix is built by averaging expression (also
# includes an option to take the median) of a single cell RNA-seq expression matrix by cluster.

# from matrix
new_ref_matrix <- average_clusters(
  mat = pbmc_matrix,
  metadata = pbmc_meta$classified, # or use metadata = pbmc_meta, cluster_col = "classified"
  if_log = TRUE                    # whether the expression matrix is already log transformed
)

head(new_ref_matrix)

# generate reference matrix from `SingleCellExperiment` or `seurat` object
new_ref_matrix_sce <- object_ref(
  input = sce_small,               # SCE object
  cluster_col = "cell_type1"       # name of column in colData containing cell identities
)

new_ref_matrix_v3 <- seurat_ref(
  seurat_object = s_small3,        # SeuratV3 object
  cluster_col = "RNA_snn_res.1"    # name of column in meta.data containing cell identities
)

tail(new_ref_matrix_v3)

sessionInfo()
