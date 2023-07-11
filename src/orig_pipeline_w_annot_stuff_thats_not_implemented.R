library(dplyr)
#BiocManager::install("Seurat")
library(Seurat)
#BiocManager::install("SeuratData")
#library(SeuratData)
library(patchwork)
#install.packages("remotes")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(scDblFinder)
library(scran)
library(BiocParallel)
#BiocManager::install("DropletUtils")
library(DropletUtils)
library(cowplot)
library(harmony)
library(ggplot2)

# Load filtered dataset from alignment

#gamm.data <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S1_mm_mt/GeneFull/filtered/")
# h5 file
#gamm.data<- Read10X_h5("/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/molecule_info.h5", use.names = TRUE, unique.features = TRUE)
# Initialize the Seurat object with the raw (non-normalized data).
#gamm <- CreateSeuratObject(counts = gamm.data, project = "gamm_s1", min.cells = 3, min.features = 200)

# load dataset with different batches
gamm.data1 <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-1_mm_mt/GeneFull/filtered/")
gamm.data2 <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-2_mm_mt/GeneFull/filtered/")

# Lets examine a few genes in the first thirty cells
#gamm.data[c("ALDH1A1", "C9orf40", "NMRK1"), 1:30]
# sparse matrix with reduced memory
#dense.size <- object.size(as.matrix(gamm.data))
#dense.size
#sparse.size <- object.size(gamm.data)
#sparse.size
#dense.size/sparse.size

# QC and filtering

# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes and ambient mRNA
# Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination

## ambient RNA count correction- SoupX- this has to happen first before other filtering
# because raw and filtered counts/ cells have to be the same
library(SoupX)
# get raw data for lanes 1 and 2
gamm.data1.raw <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-1_mm_mt/GeneFull/raw/")
gamm.data2.raw <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-2_mm_mt/GeneFull/raw/")

dim(gamm.data1)
dim(gamm.data1.raw)
dim(gamm.data2)
dim(gamm.data2.raw)
# GetAssayData(object = gamm[["RNA"]], slot = "data")
# make soup channel and profile the soup with raw and filtered data
sc1 = SoupChannel(gamm.data1.raw, gamm.data1)
sc2 = SoupChannel(gamm.data2.raw, gamm.data2)
# get basic clusters
# make seurat object without filtering
gamm1 <- CreateSeuratObject(counts = gamm.data1, project = "gamm_s2-1_clusters")
gamm2 <- CreateSeuratObject(counts = gamm.data2, project = "gamm_s2-2_clusters")
# quick cluster1
gamm1    <- SCTransform(gamm1, verbose = F)
gamm1    <- RunPCA(gamm1, verbose = F)
gamm1    <- RunUMAP(gamm1, dims = 1:30, verbose = F)
gamm1    <- FindNeighbors(gamm1, dims = 1:30, verbose = F)
gamm1    <- FindClusters(gamm1, verbose = T)
# quick cluster2
gamm2    <- SCTransform(gamm2, verbose = F)
gamm2    <- RunPCA(gamm2, verbose = F)
gamm2    <- RunUMAP(gamm2, dims = 1:30, verbose = F)
gamm2    <- FindNeighbors(gamm2, dims = 1:30, verbose = F)
gamm2    <- FindClusters(gamm2, verbose = T)
# extract clusters from seurat object and add to soup channel
meta1    <- gamm1@meta.data
meta2    <- gamm2@meta.data
umap1    <- gamm1@reductions$umap@cell.embeddings
umap2    <- gamm2@reductions$umap@cell.embeddings
sc1  <- setClusters(sc1, setNames(meta1$seurat_clusters, rownames(meta1)))
head(meta1)
sc2  <- setClusters(sc2, setNames(meta2$seurat_clusters, rownames(meta2)))
# Estimate contamination fraction
sc1  = autoEstCont(sc1)
sc2 = autoEstCont(sc2)
#Genes with highest expression in background. These are often enriched for ribosomal proteins.
head(sc1$soupProfile[order(sc1$soupProfile$est, decreasing = T), ], n = 20)
# Infer corrected table of counts and round to integer
out1 = adjustCounts(sc1, roundToInt = TRUE)
head(out1)
out2 = adjustCounts(sc2, roundToInt = TRUE)
# plot soup correction for a gene
plotChangeMap(sc1, out1, DR=umap1, "RPLP1")
plotChangeMap(sc2, out2, DR=umap2, "RPLP1")
# write
DropletUtils:::write10xCounts("gamm_soupX_filtS2-1.mtx", out1)
DropletUtils:::write10xCounts("gamm_soupX_filtS2-2.mtx", out2)
# replace filtered count values with soupX values
# gamm[["RNA"]] <- SetAssayData(gamm[["RNA"]],
#                               slot = "data", 
#                               new.data = out)
# done with soupx- remove raw data to save space -only need out files
rm(gamm.data1.raw,gamm.data2.raw,gamm1,gamm2,gamm.data1,gamm.data2,meta1,meta2,umap1,umap2,sc1,sc2)
# make separate seurat objects to get SingleCellExperiments for each lane
gamm1_seu <- CreateSeuratObject(counts = out1, project = "gamm_s1-1", min.cells = 3, min.features = 200)
gamm2_seu <- CreateSeuratObject(counts = out2, project = "gamm_s1-2", min.cells = 3, min.features = 200)
sce1 <- as.SingleCellExperiment(gamm1_seu)
sce2 <- as.SingleCellExperiment(gamm2_seu)
# remove outs and seurat objects here, only need sces for next step
rm(gamm1_seu,gamm2_seu,out1,out2)

## Doublet removalwith scDblFinder
## this should happen before filtering and normalization
# dim before doublet removal
dim(sce1)
dim(sce2)
set.seed(1234)
# run scDBLfinder on each SingleCellExperiment
sce1 <- scDblFinder(sce1)
table(call=sce1$scDblFinder.class)
sce2 <- scDblFinder(sce2)
table(call=sce2$scDblFinder.class)
# filter out doublets
sce1 <- sce1[,sce1$scDblFinder.class == "singlet"]
sce2 <- sce2[,sce2$scDblFinder.class == "singlet"]
# convert back to seurat
gamm1<- as.Seurat(sce1)
gamm2 <- as.Seurat(sce2)
gamm <- merge(gamm1, y = gamm2)
unique(gamm@meta.data$orig.ident)
# save
saveRDS(gamm, file = "GAMM_S2_doubletremoved.rds")

## Filter mitochondrail genes, high and low counts

# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function
# For human: all genes starting with MT- are mitochondrial genes
# gamm[["percent.mt"]] <- PercentageFeatureSet(gamm, pattern = "^MT-")
# no percent.mt in pig because they don't start with MT

# give specific features (genes) for pig data- 13 protein coding mt genes in pig:
mt.list <- c("ND1","ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L",
             "ND4", "ND5", "ND6", "CYTB")

# check if all features are in the matrix
all(mt.list %in% rownames(gamm))
# if false, need to determine which ones
for (mt in mt.list) {print(mt)
                      print(all(mt %in% rownames(gamm)))
                      }
percent_mt <- PercentageFeatureSet(gamm, features = mt.list, assay = 'RNA')
gamm[["percent.mt"]] <- percent_mt

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts
# Visualize QC metrics as a violin plot
VlnPlot(gamm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "percent.mt")+ theme(axis.text.x = element_text(angle = 45, hjust=1))
plot2 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ theme(axis.text.x = element_text(angle = 45, hjust=1))
plot1 + plot2

# subset data -  filter out cells that have unique feature counts over 2,500 or less than 200 and > 5% mitochondrial counts
gamm <- subset(gamm, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
# subset data by percent- maybe this is better?
#counts <- GetAssayData(seurat_obj, slot="counts", assay="RNA")   
#genes.percent.expression <- rowMeans(counts>0 )*100   
#genes.filter <- names(gene.percent.expressed[gene.percent.expressed>1])  #select genes expressed in at least 1% of cells
#counts.sub <- counts[genes.filter,]
#new_seurat_object <- CreateSeuratObject(counts=counts.sub)

# after you get lane 1 & 2 numbers, delete seurat object for them to save space
rm(gamm1,gamm2,sce1,sce2)
# replot
plot1 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "percent.mt")+ theme(axis.text.x = element_text(angle = 45, hjust=1))
plot2 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ theme(axis.text.x = element_text(angle = 45, hjust=1))
plot1 + plot2

# Normalize the data

# default is “LogNormalize” that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), 
# and log-transforms the result. 
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# scran normalization
library(scran)
sce <- as.SingleCellExperiment(gamm)
# get raw UMI counts
counts(sce)[40:50, 40:50]
# get log Normalized counts form Seurat
logcounts(sce)[40:50, 40:50]
#genes
rownames(sce)[1:5]
#cells
colnames(sce)[1:5]
# get cell metadata
colData(sce)[1:5, ]
# get gene metadata
rowData(sce)[1:5, ]
# cluster the cells for scran
clusters <- quickCluster(sce,
                         use.ranks = FALSE, # suggested by the authors
                         min.size = 100) # require at least 100 cells per cluster

table(clusters)
# Calculate size factors per cell
sce <- computeSumFactors(sce, 
                         clusters = clusters, 
                         min.mean = 0.1) # ignore low abundance genes 
summary(sizeFactors(sce))
#use scuttle to get lognormcounts function
library(scuttle)
# apply size factors to generate log normalized data # using log norm counts instead of this function -- sce <- normalize(sce)
sce <- logNormCounts(sce)

logcounts(sce)[40:50, 40:50]
# replace seurat normalized values with scran
gamm[["RNA"]] <- SetAssayData(gamm[["RNA"]],
                            slot = "data", 
                            new.data = logcounts(sce))

gamm$sizeFactors <- sizeFactors(sce)
UpdateSeuratObject(gamm)
# gene expression with reg counts
VlnPlot(gamm, "ECHS1", slot = "counts")
# gene expression with scran norm counts
VlnPlot(gamm, "ECHS1", slot = "data")
# save object
saveRDS(gamm, file = "GAMM_S2_filtered-norm.rds")


# ID highly variable features (feature selection)
# calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 
# seurat method
gamm <- FindVariableFeatures(gamm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gamm), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(gamm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot1
# feature selection with scry
# scry deviance for feature selection which works on raw counts [Germain et al., 2020]. Deviance can be computed in closed form and quantifies whether genes show a constant expression profile across cells as these are not informative. 
# BiocManager::install('scry')
library(scry)
m <- GetAssayData(gamm, slot = "counts", assay = "RNA")
devi <- scry::devianceFeatureSelection(m)
dev_ranked_genes <- rownames(gamm)[order(devi, decreasing = TRUE)]
topdev <- head(dev_ranked_genes, 2000)
# replace variable features with the deviance ranked genes
VariableFeatures(gamm) <- topdev
# VariableFeatures(gamm)
UpdateSeuratObject(gamm)
devtop10 <- head(VariableFeatures(gamm), 10)
devtop10
plot1<-VariableFeaturePlot(gamm)
plot2<-LabelPoints(plot = plot1, points = devtop10, repel = TRUE)
plot1 + plot2
# Scale the data
# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
#   
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(gamm)
gamm <- ScaleData(gamm, features = all.genes)

# Remove unwanted variation
# also use the ScaleData() function to remove unwanted sources of variation from a 
# single-cell dataset. For example, we could ‘regress out’ heterogeneity associated with 
# (for example) cell cycle stage, or mitochondrial contamination. # get list of cell cycle genes
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
saveRDS(gamm, file = "GAMM_S2_featSel-scaled.rds")
rm(m)
# Perform linear dimensional reduction on the variable features

gamm <- RunPCA(gamm, features = VariableFeatures(object = gamm))
# Examine and visualize PCA results a few different ways
print(gamm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gamm, dims = 1:2, reduction = "pca")
DimPlot(gamm, reduction = "pca")
#DimHeatmap(gamm, dims = 1, cells = 500, balanced = TRUE)
dh<- DimHeatmap(gamm, dims = 1:15, cells = 500, balanced = TRUE, fast=FALSE, combine=TRUE, nfeatures=5)
dh

# Harmony for batch control
# define batches
#gamm@meta.data$lanes <- c(rep("lane1", lane1), rep("lane2", lane2))
# check PCA plot
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gamm, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = gamm, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
# run harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
gamm <- gamm %>% 
  RunHarmony(group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50, plot_convergence = TRUE)
# access harmony embeddings
harmony_embeddings <- Embeddings(gamm, 'harmony')
harmony_embeddings[1:5, 1:5]
# check out new harmony plot
options(repr.plot.height = 5, repr.plot.width = 12)
p1a <- DimPlot(object = gamm, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2b <- VlnPlot(object = gamm, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1a,p2b)

# Determine the ‘dimensionality’ of the dataset
# to overcome technical noise, Seurat clusters cells based on their PCA scores, 
# JackStraw randomly permutes a subset (1%) of the data and reruns PCA, constructing a null
# distribution of feature scores. Significant PCs have strong enrichment of low p-val features
gamm <- JackStraw(gamm, num.replicate = 100)
gamm <- ScoreJackStraw(gamm, dims = 1:20)
# The JackStrawPlot() function provides a visualization tool for comparing the 
# distribution of p-values for each PC with a uniform distribution (dashed line). 
# ‘Significant’ PCs will show a strong enrichment of features with low p-values 
#(solid curve above the dashed line).
JackStrawPlot(gamm, dims = 1:20)
#  alternative heuristic method generates an ‘Elbow plot’: a ranking of principle 
# components based on the percentage of variance explained by each one (ElbowPlot() 
# function). In this example, we can observe an ‘elbow’ to suggest how many PCs capture
# the majority of the true signal. can use reduction='harmony' to use harmony embeddings
ElbowPlot(gamm, reduction='harmony')
ElbowPlot(gamm, reduction='pca')
# Cluster the Cells using only genes from chosen PCs
# to use harmony embeddings pass reduction='harmony'

# K-nearest neighbor (KNN) graph
gamm <- FindNeighbors(gamm, dims = 1:15, reduction='harmony')
# # Run non-linear dimensional reduction (UMAP/tSNE)
# The goal of these algorithms is to learn the underlying manifold of the data in order 
# to place similar cells together in low-dimensional space. Cells within the graph-based 
# clusters determined above should co-localize on these dimension reduction plots. 
# If you haven't installed UMAP: 
# reticulate::py_install(packages ='umap-learn')
gamm <- RunUMAP(gamm, dims = 1:15, reduction = "harmony")
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(gamm, reduction = "umap", group.by = "orig.ident", pt.size = .1)#, split.by = 'lanes')

# Clustering
# note need to pip install leiden and pandas for FindClusters to run
#reticulate::py_install(packages ='leidenalg') #pip install leidenalg
#reticulate::py_install(packages ='pandas')
gamm <- FindClusters(gamm, resolution = 0.5, algorithm = "leiden",reduction = "harmony")
# Look at cluster IDs of the first 5 cells
head(Idents(gamm), 5)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# plot individual clusters
DimPlot(gamm, reduction = "umap", label = T)
# save object
saveRDS(gamm, file = "GAMM_S2_dimRed-harmony-clusters.rds")
# read back in the object
gamm <- readRDS(file = "GAMM_dimRed-harmony-clusters.rds")

# get r environment
install.packages("renv")
renv::init()

# Finding differentially expressed features (cluster biomarkers)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
gamm.markers <- FindAllMarkers(gamm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gamm.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  print(n=36)
# ROC test for cluster markers, 0= random, 1= perfect)
#cluster1.markers <- FindMarkers(gamm, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#head(cluster1.markers, n = 5)
# visualize marker expression
VlnPlot(gamm, features = c("APOA1", "VIM"))
# you can plot raw counts as well
VlnPlot(gamm, features = c("APOA1", "VIM"), slot = "counts", log = TRUE)
# umap plot highlighting gene expression
plot<-FeaturePlot(gamm, features = c("APOA1", "VIM", "CENPE", "AGMO"),slot='data',reduction='umap')
#LabelClusters(plot = plot, id = 'ident',clusters= c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
plot
# or plot together with cluster map
plot1 <- UMAPPlot(gamm, group.by="orig.ident")
plot2 <- UMAPPlot(gamm, label = T)
plot3 <- FeaturePlot(gamm, c("APOA1", "VIM", "CENPE", "AGMO"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers
gamm.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# write out top10 markers
top10.df <- as.data.frame(top10)
write.table(top10.df, file="gamm.DE.markers.top10_S2.txt",quote = F, sep = "\t", row.names=F)
# run heatmap
DoHeatmap(gamm, features = top10$gene) + NoLegend()
# write out all markers
gamm.markers.df <- as.data.frame(gamm.markers)
write.table(gamm.markers.df, file="gamm.DE.markers_S2.txt",quote = F, sep = "\t", row.names=F)

# or plot together with cluster map
plot1 <- UMAPPlot(gamm, group.by="orig.ident")
plot2 <- UMAPPlot(gamm, label = T)
plot3 <- FeaturePlot(gamm, c("APOA1", "VIM", "CENPE", "AGMO"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers
gamm.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# write out top10 markers
top10.df <- as.data.frame(top10)
write.table(top10.df, file="gamm.DE.markers.top10_S2.txt",quote = F, sep = "\t", row.names=F)
# run heatmap
DoHeatmap(gamm, features = top10$gene) + NoLegend()
# write out all markers
gamm.markers.df <- as.data.frame(gamm.markers)
write.table(gamm.markers.df, file="gamm.DE.markers_S2.txt",quote = F, sep = "\t", row.names=F)

# Assigning cell type identity to clusters
# read in known GAMM retinoid markers
known.markers<- read.csv2("../GammLab_Retinal-specific_genelist.txt",sep = "\t", header = TRUE)
known.markers.df <- as.data.frame(known.markers)
# match with any DE markers from data by merging dataframes
marker_df <- merge(gamm.markers.df,known.markers.df,by="gene")
# write out marker df with known DE markers
write.table(marker_df, file="gamm.knownDE.markers_S2.txt",quote = F, sep = "\t", row.names=F)

# cell.type subsets
# first get unique cell type vector
cell.types <- unique(marker_df$Cell.type)
print(cell.types)
# check gene number between known markers and DE known markers
genesk<- unique(known.markers.df$gene)
k<-length(genesk)
genesDEk<- unique(marker_df$gene)
DEk<-length(genesDEk)
percent.markers.de<-(DEk/k)*100
percent.markers.de
# check a specific marker
df<-subset(marker_df,Cell.type == "Synaptic marker",select = c('gene'))
length(df$gene)
df2<-subset(known.markers.df,Cell.type == "Synaptic marker",select = c('gene'))
length(df2$gene)
# subset all and plot using for loop
library(ggplot2)
for (i in 1:length(cell.types)){
new_df<- subset(marker_df,Cell.type == cell.types[i],select = c('gene','Cell.type','cluster'))
new_vec<- unique(as.vector(new_df$gene))
pdf(paste0(cell.types[i], "_featureplot.pdf", collapse = ""),        # File name
    width = 8, height = 6, # Width and height in inches
    bg = "white")          # Background color
# umap plot highlighting gene expression
print(FeaturePlot(gamm, features = new_vec))

dev.off()

pdf(paste0(cell.types[i], "_dotplot.pdf", collapse = ""),         # File name
    width = 8, height = 6, # Width and height in inches
    bg = "white")          # Background color

# expression dot plot
dot.plot<-DotPlot(object = gamm, features = new_vec)
print(dot.plot +labs(title = cell.types[i]))

dev.off()

}
# score markers by pairwise comparisons across clusters 
library(scran)
sce.gamm <- as.SingleCellExperiment(gamm)
sce.gamm
marker.info <- scoreMarkers(sce.gamm, sce.gamm@colData@listData$seurat_clusters)
marker.info
colnames(marker.info[["1"]]) # statistics for cluster 1.
# rank marker by AUC
clust_1 <- marker.info[["1"]]
ordered <- clust_1[order(clust_1$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4]) # showing basic stats only, for brevity.
library(scater)
plotExpression(sce.gamm, features=head(rownames(ordered)), 
               x="seurat_clusters", colour_by="seurat_clusters")
# AUC only clster 1
auc.only <- clust_1[,grepl("AUC", colnames(clust_1))]
auc.only[order(auc.only$mean.AUC,decreasing=TRUE),]
# cohen's d (standardized log FC: the difference in the mean log-expression
# between groups is scaled by the average standard dev across groups)
cohen.only <- clust_1[,grepl("logFC.cohen", colnames(clust_1))]
cohen.only[order(cohen.only$mean.logFC.cohen,decreasing=TRUE),]
# using median cohen d:
ordered <- clust_1[order(clust_1$median.logFC.cohen,decreasing=TRUE),]
head(ordered[,1:4]) # showing basic stats only, for brevity.
# plot
plotExpression(sce.gamm, features=head(rownames(ordered)), 
               x="seurat_clusters", colour_by="seurat_clusters")
# using ranked cohen d
# top 5 genes (T=5)
ordered <- clust_1[order(clust_1$rank.logFC.cohen),]
top.ranked <- ordered[ordered$rank.logFC.cohen <= 5,]
rownames(top.ranked)
# heatmap
plotGroupedHeatmap(sce.gamm, features=rownames(top.ranked), group="seurat_clusters", 
                   center=TRUE, zlim=c(-3, 3))

# obtain full effects
marker.info <- scoreMarkers(sce.gamm, sce.gamm@colData@listData$seurat_clusters, full.stats=TRUE)
clust_1 <- marker.info[["1"]]
clust_1$full.AUC
# identify the genes that distinguish cluster 1 from other clusters with high VIM expression
vim.high <- c("13", "15", "16", "17") # based on inspection of the previous Figure.
subset <- clust_1$full.AUC[,colnames(clust_1$full.AUC) %in% vim.high]
to.show <- subset[computeMinRank(subset) <= 10,]
to.show

plotGroupedHeatmap(sce.gamm[,sce.gamm@colData@listData$seurat_clusters %in% vim.high],
                   features=rownames(to.show), group="seurat_clusters", center=TRUE, zlim=c(-3, 3))

colLabels(sce.gamm)
colLabels(sce.gamm)<-sce.gamm@colData@listData$seurat_clusters
plotGroupedHeatmap(sce.gamm[,colLabels(sce.gamm) %in% vim.high],
                   features=rownames(to.show), group="label", center=TRUE, zlim=c(-3, 3))
# assign cluster labels based on markers
marker.info <- scoreMarkers(sce.gamm, sce.gamm@colData@listData$seurat_clusters)
# write all marker info
marker.info.df <- as.data.frame(marker.info)
write.table(marker.info.df, file="marker.info.S1.txt",quote = F, sep = "\t", row.names=T)
# get clusters
clusters<- unique(gamm@meta.data$seurat_clusters)
clusters<- as.vector(clusters)
clusters
# read in known GAMM retinoid markers
known.markers<- read.csv2("../GammLab_Retinal-specific_genelist.txt",sep = "\t", header = TRUE, row.names = 1)
known.markers.df <- as.data.frame(known.markers)
# for loop to get info on each cluster
for (i in 1:length(clusters)){
  clust<- marker.info[[clusters[i]]]
  ordered <- clust[order(clust$median.logFC.cohen,decreasing=TRUE),]
  top100 <- ordered[1:100,]
  top100<-as.data.frame(top100)
  #top<- subset(clust, median.logFC.cohen >= 0.25)
  #top<- as.data.frame(top)
  # match with any DE markers from data by merging dataframes
  marker_df <- merge(top100,known.markers.df,by='row.names')
  if (nrow(marker_df) == 0){
    print(paste0("This data frame is empty: ", clusters[i]))
  }else{
  # write out marker df with known DE markers
  write.table(marker_df, file=paste0("gamm.knownDE.markers_S1_clust_",clusters[i],".txt", collapse = ""),quote = F, sep = "\t", row.names=F)
  # subset data
  new_df<- marker_df[,c('Row.names','rank.logFC.cohen','Cell.type')]
  new_vec<- unique(as.vector(new_df$Row.names))
  # # make plots
  pdf(paste0(clusters[i],"_featureplot.pdf", collapse = ""),        # File name
       width = 8, height = 11, # Width and height in inches
       bg = "white")          # Background color
   # umap plot highlighting gene expression
  print(FeaturePlot(gamm, features = new_vec), label=T)
  
  dev.off()
  # top 10 ranked
  new_df.ordered <- new_df[order(new_df$rank.logFC.cohen),]
  new_df.ordered<- subset(new_df.ordered, rank.logFC.cohen < 11)
  write.table(new_df.ordered, file=paste0("gamm.knownDE.markers_S1_clust_top10",clusters[i],".txt", collapse = ""),quote = F, sep = "\t", row.names=F)
  new_vec2<- unique(as.vector(new_df.ordered$Row.names))
  # make plots
  pdf(paste0(clusters[i],"_featureplot_top10ranks.pdf", collapse = ""),        # File name
      width = 8, height = 11, # Width and height in inches
      bg = "white")          # Background color
  # umap plot highlighting gene expression
  print(FeaturePlot(gamm, features = new_vec2), label=T)
  dev.off()
  # dot plots
  pdf(paste0(clusters[i], "_dotplot.pdf", collapse = ""),         # File name
       width = 8, height = 6, # Width and height in inches
       bg = "white")          # Background color
   dot.plot<-DotPlot(object = gamm, features = new_vec)
  print(dot.plot +labs(title = paste0("cluster_",clusters[i])))
   
  dev.off()
}}

# Annotate clusters based on markers
#S2
new.cluster.ids <- c("Retinal Prog-Muller glia-1","Pan PR-Synaptic-Neuronal proj", 
                    "Pan PR-Synaptic-Neuronal proj-Rod-1","Pan PR-Rod-Neuronal proj-1",
                    "Pan PR-Rod-Neuronal proj-2","Pan PR-Synaptic-Neuronal proj-Rod-2",
                    "mix1","Bipolar cells fetal-Synaptic-1","Pan PR-Rod-Neuronal proj-3",
                    "Phototransduction?","Phototransduction?-PanPR","Cone-PanPR-Neuronal proj-Ganglion cell(fetal)",
                    "mix2","Retinal Progenitor","Muller glia-Retinal prog-Rod-PanPR","Retinal Prog-Muller glia-2",
                    "Retinal Prog-Muller glia-3","Bipolar cells fetal-Synaptic-2")
#S1
new.cluster.ids <- c("Retinal Prog-Muller glia-1", "Neuronal proj-Ganglion cell-PanPRs", 
                     "Neuronal proj-PanPRs","Ganglion cell-1","Ganglion cell-2",
                     "Ganglion cell (fetal)-1","Retinal Progenitor","Neuronal proj-Ganglion cell",
                     "Retinal Prog-Muller glia-2","Ganglion cell (fetal)-2",
                     "Ganglion cell (fetal)-Synaptic markers","Retinal Prog-Muller glia-Microglia (fetal)",
                     "Neuronal projection","Retinal Progenitor (fetal)","Retinal Prog-Muller glia-3")
names(new.cluster.ids) <- levels(gamm)
gamm <- RenameIdents(gamm, new.cluster.ids)
DimPlot(gamm, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# save object
saveRDS(gamm, file = "GAMM_S1_labeled-clusters.rds")
sessionInfo()

