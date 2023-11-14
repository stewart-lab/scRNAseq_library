
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
use_condaenv(condaenv = 'scRNAseq_best', required = TRUE)
## to install SeuratDisk
# first install hdf5r:
# conda install -c conda-forge r-hdf5r
# Next install SeuratDisk:
# conda install -c conda-forge r-seuratdisk
# or
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# note for me conda install -c conda-forge r-seuratdisk did not work, but the github install did

library(SeuratDisk)

setwd("/w5home/bmoore/scRNAseq/LiFangChu/fluidigm_gup_expr_results/output_20231114_102348/")


# read in Seurat object
seurat.obj<- readRDS(file = "clustered_seurat_obj.rds")

# add in metadata from previous analysis
# metadata.gamm <- read.csv("../gamm_manual_annot_metadata_c0.5.txt", row.names = 1, 
#                           header = TRUE, sep = "\t")
# # check metadata with query data
# colnames(gamms2@assays$RNA@data)[1:10]
# rownames(metadata.gamm)[1:10]
# # add metadata to query data
# gamms2 <- AddMetaData(gamms2, metadata.gamm)


# first save seurat as h5 seurat file
SaveH5Seurat(seurat.obj, filename = "clustered_seurat_obj.h5Seurat")
# then convert to h5ad
Convert("clustered_seurat_obj.h5Seurat", dest = "h5ad")

# can now be read in by scanpy:
# import scanpy
# adata = scanpy.read_h5ad("pbmc3k.h5ad")