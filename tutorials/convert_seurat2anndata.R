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
gamms2<- readRDS(file = "gamms2_cca_pred.rds")

# metadata
# add in metadata from previous analysis
metadata.gamm <- read.csv("../gamm_manual_annot_metadata_c0.5.txt", row.names = 1, 
                          header = TRUE, sep = "\t")
# check metadata with query data
colnames(gamms2@assays$RNA@data)[1:10]
rownames(metadata.gamm)[1:10]
# add metadata to query data
gamms2 <- AddMetaData(gamms2, metadata.gamm)

# save rds object
saveRDS(gamms2, file= "gamms2_clustifyr.rds")
# first save seurat as h5 seurat file
SaveH5Seurat(gamms2, filename = "gamms2_cca_pred.h5Seurat")
# then convert to h5ad
Convert("gamms2_cca_pred.h5Seurat", dest = "h5ad")

# can now be read in by scanpy:
# import scanpy
# adata = scanpy.read_h5ad("pbmc3k.h5ad")