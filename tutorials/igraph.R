install.packages("igraph")
BiocManager::install("graph")
library(igraph)
library(graph)
library(dplyr)
library(Seurat)
library(patchwork)
library(BiocParallel)
library(cowplot)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)

# setwd
setwd("/Users/bmoore/Desktop/scRNAseq/GAMM/GAMM_S1/output_20230921_142919/")
setwd("/Users/bmoore/Desktop/scRNAseq/GAMM/GAMM_S2/output_20230830_155530/")
# read in a seurat object
gamms2 <- readRDS(file = "gamms1_clustifyr.rds")
gamms2 <- readRDS(file = "gamms2_clustifyr.rds")
gamms2@graphs$RNA_nn
#gamms2@graphs$RNA_snn
summary(gamms2@graphs$RNA_nn)
# convert matrix to graph
# seurat FindNeighbors default is an undirected graph
g1<-graph_from_adjacency_matrix(gamms2@graphs$RNA_nn, mode="undirected")
#g2<-graph_from_adjacency_matrix(gamms2@graphs$RNA_nn, mode="directed")
# get a layout
g.layout<-layout_with_fr(g1)
#g2.layout<-layout_with_fr(g2)
# get metadata and color
# check idents
table(gamms2@active.ident)
# rename clusters if necessary
gamms2 <- RenameIdents(object = gamms2, `Muller Glia - Retinal Prog` = "Retinal Prog - Muller Glia"
                       )
# add metadata
gamms2 <- AddMetaData(gamms2, metadata = gamms2@active.ident, col.name= "CellType")
table(gamms2$CellType)
# read in metadata
metadat <- gamms2@meta.data
# color scale for Celltype
table(metadat$CellType)
# if have 10 different stages, define 10 colors for those.
sorted_celltypes <- sort(unique(metadat$CellType))
sorted_celltypes
coldef.celltype<-c("purple","red","orange","green","seagreen","cyan","dodgerblue","slateblue1","darkolivegreen1")
names(coldef.celltype) <- sorted_celltypes
col.celltype <- coldef.celltype[metadat$CellType]
# s1 annotation
coldef.celltype<-c("slateblue1","green","purple","seagreen")
sorted_celltypes <- sort(unique(metadat$CellType))
names(coldef.celltype) <- sorted_celltypes
col.celltype <- coldef.celltype[metadat$CellType]
table(col.celltype)
table(metadat$CellType)
coldef.celltype
#clustifyr annotation
table(metadat$clustifyr_call_consol_type)

# 5 cell types called
coldef.celltype<-c("purple","orange","dodgerblue","slateblue1","cyan")
sorted_celltypes <- sort(unique(metadat$clustifyr_call_consol_type))
names(coldef.celltype) <- sorted_celltypes
col.celltype <- coldef.celltype[metadat$clustifyr_call_consol_type]
# s1:
table(metadat$clustifyr_call_type)
# 4 cell types
sorted_celltypes <- sort(unique(metadat$clustifyr_call_type))
coldef.celltype<-c("yellow","orange","purple","orchid1")
names(coldef.celltype) <- sorted_celltypes
col.celltype <- coldef.celltype[metadat$clustifyr_call_type]
# clustifyr pred type
table(metadat$clustify_pred_type)
# 6 cell types called
coldef.celltype<-c("dodgerblue","orange","green","cyan","slateblue1","purple")
sorted_celltypes <- sort(unique(metadat$clustify_pred_type))
names(coldef.celltype) <- sorted_celltypes
col.celltype <- coldef.celltype[metadat$clustify_pred_type]

# plot
plot.igraph(g1, layout=g.layout,vertex.color=col.celltype,
            vertex.size=3,vertex.label=NA,main="scRNA NN graph")

# legend
legend("topleft",names(coldef.celltype),col=coldef.celltype,
       pch=16,cex=0.5,bty='n')

# plot g2
plot.igraph(g2, layout=g2.layout,vertex.color=col.celltype,
            vertex.size=3,vertex.label=NA,main="scRNA NN graph directed")

# legend
legend("topright",names(coldef.celltype),col=coldef.celltype,
       pch=16,cex=0.5,bty='n')


## make knn graph from expr data
# read in metadata
metadat <- gamms2@meta.data
# read in normalized expression
expdat <- gamms2@assays$data
