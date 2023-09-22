# install packages
BiocManager::install("clustifyr")
# load libraries
library(clustifyr)
library(ggplot2)
library(cowplot)
library(Seurat)
# setwd
setwd("/Users/bmoore/Desktop/scRNAseq/GAMM/GAMM_S1/output_20230921_142919")
# load data set
gamms2 <- readRDS(file = "seurat_obj_labeled.rds")
table(gamms2$orig.ident)
table(gamms2$seurat_clusters)
table(gamms2@active.ident)
gamms2 <- AddMetaData(gamms2, metadata = gamms2@active.ident, col.name= "CellType")
table(gamms2$CellType)
# load marker list
markers <- read.csv2("../../known_markers/GammLab_Retinal-specific_genelist_mod.txt", sep="\t", header = TRUE)
# consolidated markers
markers <- read.csv2("../../known_markers/Gamm_lab_Consolidated_markerList.txt", sep="\t", header = TRUE)
# subset
markers_df <- data.frame(markers$gene,markers$Cell.type)
# rename columns and matrixize markers
colnames(markers_df) <- c("gene","cluster")
#marker_mtx <- matrixize_markers(markers_df)

# clustifyr annotation from marker list
# Available metrics include: "pct", "hyper", "jaccard", "spearman", "gsea"
list_res <- clustify_lists(
  input = gamms2,             # matrix of normalized single-cell RNA-seq counts
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  marker = markers_df,                 # list of known marker genes
  metric = "pct", # test to use for assigning cell types
  marker_inmatrix = FALSE,
  obj_out = FALSE
)

# View as heatmap, or plot_best_call
plot_cor_heatmap(
  cor_mat = list_res,              # matrix of correlation coefficients from clustify_lists()
  cluster_rows = TRUE,            # cluster by row?
  cluster_columns = TRUE,         # cluster by column?
  legend_title = "pct"     # title of heatmap legend
)

# Call cell types
list_res2 <- cor_to_call(
  cor_mat = list_res,              # matrix correlation coefficients
  cluster_col = "seurat_clusters"  # name of column in meta.data containing cell clusters
)

# Insert into metadata dataframe as "clustifyr_call" column
gamms2.1 <- call_to_metadata(
  res = list_res2,                 # data.frame of called cell type for each cluster
  metadata = gamms2@meta.data, 
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  rename_prefix = "clustifyr_call_full"          # set a prefix for the new column
)
# add to seurat object metadata
gamms2 <- AddMetaData(gamms2, metadata = gamms2.1)
# plot
pc <- DimPlot(gamms2, reduction = "umap", group.by = "clustifyr_call_full_type", label = TRUE,
              label.size = 3, repel = TRUE) + ggtitle("GammS1 clustifyr annotated full labels") +
              guides(fill = guide_legend(label.theme = element_text(size = 8)))
pdf("clustifyr_full_annotation_umap.pdf", width = 11, height = 6)
print(pc)
dev.off()

saveRDS(gamms2, file= "gamms1_clustifyr.rds")

## clustifyr annotation from reference

# load annotated human data
human_D205.seurat<- readRDS(file = "../../human_ref/output_preprocess20230912_144047_cc/human_D205_umap.rds")
table(human_D205.seurat$type)
# subset
human_D205.seurat<- subset(human_D205.seurat, subset = type != c('AC2'))
human_D205.seurat<- subset(human_D205.seurat, subset = type != c('miG'))
human_D205.seurat<- subset(human_D205.seurat, subset = type != c('T2'))
# get ref matrix
new_ref_matrix <- seurat_ref(
  seurat_object = human_D205.seurat,        # SeuratV3 object
  cluster_col = "type"    # name of column in meta.data containing cell identities
)
# run clustify (spearman corr is default)
res <- clustify(
  input = gamms2,       # a Seurat object
  ref_mat = new_ref_matrix,    # matrix of RNA-seq expression data for each cell type
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  obj_out = TRUE,      # output Seurat object with cell type inserted as "type" column
  rename_prefix= "clustify_pred",
  n_genes= 2000
)
table(res$clustify_pred_type)

pc <- DimPlot(res, reduction = "umap", group.by = "clustify_pred_type", label = TRUE,
              label.size = 3, repel = TRUE) + ggtitle("GammS2 clustifyr predicted labels- spearmans") +
  guides(fill = guide_legend(label.theme = element_text(size = 8)))
pdf("clustifyr_predicted_labels_umap.pdf", width = 10, height = 6)
print(pc)
dev.off()
table(res$clustify_pred_r)
saveRDS(res, file= "gamms2_clustifyr.rds")
