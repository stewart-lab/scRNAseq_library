library(Seurat)
library(SeuratData)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

# install PBMC data
# InstallData("pbmc3k")
# data("pbmc3k.final")
# pbmc3k.final

# read in Seurat object
gamms2<- readRDS(file = "gamms2_clustifyr.rds")

# first save seurat as h5 seurat file
SaveH5Seurat(gamms2, filename = "gamms2_clustifyr.h5Seurat")
# then convert to h5ad
Convert("gamms2_clustifyr.h5Seurat", dest = "h5ad")

# can now be read in by scanpy:
# import scanpy
# adata = scanpy.read_h5ad("pbmc3k.h5ad")