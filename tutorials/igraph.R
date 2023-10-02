install.packages("igraph")
BiocManager::install("graph")
library(igraph)
library(graph)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(BiocParallel)
library(cowplot)
library(harmony)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
#library(scPred)

# setwd
setwd("/Users/bmoore/Desktop/scRNAseq/GAMM/GAMM_S2/")
# read in a seurat object
gamms2 <- readRDS(file = "output_20230830_155530/gamms2_clustifyr.rds")
gamms2@graphs$RNA_nn
gamms2@graphs$RNA_snn
summary(gamms2@graphs$RNA_nn)
# convert matrix to graph
g1<-graph_from_adjacency_matrix(gamms2@graphs$RNA_nn, mode="undirected")
g2<-graph_from_adjacency_matrix(gamms2@graphs$RNA_snn, mode="undirected")
# get a layout
g.layout<-layout_with_fr(g1)
g2.layout<-layout_with_fr(g2)
# get metadata and color
# check idents
table(gamms2@active.ident)
# rename clusters if necessary
gamms2 <- RenameIdents(object = gamms2, `Muller Glia - Retinal Prog` = "Retinal Prog - Muller Glia",
                       `Cones - Pan PRs` = "Cones")
# add metadata
gamms2 <- AddMetaData(gamms2, metadata = gamms2@active.ident, col.name= "CellType")
table(gamms2$CellType)
# read in metadata
metadat <- gamms2@meta.data
# color scale for Celltype
table(metadat$CellType)
# if have 10 different stages, define 10 colors for those.
coldef.celltype<-c("black","red","green","blue","cyan","yellow","orange","pink","purple")
names(coldef.celltype) <- unique(metadat$CellType)
col.celltype <- coldef.celltype[metadat$CellType]

# plot
plot.igraph(g1, layout=g.layout,vertex.color=col.celltype,
            vertex.size=3,vertex.label=NA,main="scRNA NN graph")

# legend
legend("topright",names(coldef.celltype),col=coldef.celltype,
       pch=16,cex=0.5,bty='n')

# plot g2
plot.igraph(g2, layout=g2.layout,vertex.color=col.celltype,
            vertex.size=3,vertex.label=NA,main="scRNA sNN graph")

# legend
legend("topright",names(coldef.celltype),col=coldef.celltype,
       pch=16,cex=0.5,bty='n')


## make knn graph from expr data
# read in metadata
metadat <- gamms2@meta.data
# read in normalized expression
expdat <- gamms2@assays$data
