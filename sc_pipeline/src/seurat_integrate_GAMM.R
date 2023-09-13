---
title: "seurat_integration_for_cross_species"
author: "Beth Moore"
date: "2023-07-18"
output: html_document
---
# Preprocessing for cross species analysis
# to have cross species, we need ortholog genes
# we then need to subset both reference and query by these genes
# reference data will also be subseted by its metadata

# Load libraries, set wd, and source functions
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
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
BETHS_GAMM_DATA_DIR <- "/w5home/bmoore/scRNAseq/GAMM/"
use_condaenv("scRNAseq_best")
setwd(paste0(BETHS_GAMM_DATA_DIR, "human_data/reh_cellrep_2020/"))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output_", timestamp)
dir.create(output, showWarnings = FALSE)
file.copy("/w5home/bmoore/scRNAseq_library/config.json", file.path(output, "config.json"))
config <- fromJSON(file.path(output, "config.json"))
output <- paste0(output, "/")
source("/w5home/bmoore/scRNAseq_library/src/sc_pipeline_functions.R")
# do preprocessing on human and pig data (preprocess_crossspecies.Rmd)
# both reference and query should undergo normalization and find variable features,
# then only reference is scaled and undegoes dimentionality reduction

# read in seurat objects that were preprocessed
gamms2<- readRDS(file = "GAMM_S2_norm.rds")
genelist<- as.data.frame(rownames(x = gamms2))
colnames(genelist)[1] <- "gene"
metadata<- read.table("GammLab_Retinal-specific_genelist.txt", header=TRUE, sep="\t")
gamm.object.metadata <- merge(genelist,metadata,by="gene")
write.table(gamm.object.metadata, file="gamm_markers_in_orthologous_object.txt",col.names=TRUE, sep="\t",quote=F )
human_D205.seurat<- readRDS(file = "human_D205_umap.rds")
# prep for integration
ret.list<- list(gamm,human_D205.seurat)
features <- SelectIntegrationFeatures(object.list = ret.list, nfeatures = 3000)
ret.list <- PrepSCTIntegration(object.list = ret.list, anchor.features = features)

# # mapping and annotating query datasets

# reference umap plot
p2 <- DimPlot(human_D205.seurat, reduction = "umap", group.by = "type", label = TRUE, repel = TRUE) +
  NoLegend()
p2
# get query data and find anchors
gamm.anchors <- FindTransferAnchors(reference = human_D205.seurat, query = gamm,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = gamm.anchors, refdata = human_D205.seurat$type,
                            dims = 1:30)
gamm <- AddMetaData(gamm, metadata = predictions)

# evaluate how well our predicted cell type annotations match the full reference. 
# need manual annotation data
#pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
#table(pancreas.query$prediction.match)
# examine some canonical cell type markers for specific pancreatic islet cell populations.
table(gamm$predicted.id)
# need marker genes
#VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")
# projection of a query onto the reference UMAP structure.
human_D205.seurat <- RunUMAP(human_D205.seurat, dims = 1:30, reduction = "pca", return.model = TRUE)
gamm <- MapQuery(anchorset = gamm.anchors, reference = human_D205.seurat, query = gamm,
                           refdata = list(celltype = "type"), reference.reduction = "pca", reduction.model = "umap")

# visualize
p1 <- DimPlot(human_D205.seurat, reduction = "umap", group.by = "type", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(gamm, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
