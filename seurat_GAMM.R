library(dplyr)
#BiocManager::install("Seurat")
library(Seurat)
#BiocManager::install("SeuratData")
#library(SeuratData)
library(patchwork)
#install.packages("remotes")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#library(DoubletFinder)
library(scran)
library(BiocParallel)
#BiocManager::install("DropletUtils")
library(DropletUtils)
library(cowplot)
library(harmony)

# starting directory
dir.name <- "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/"
# Load filtered dataset from alignment

# this is for one output of both lanes already combined during alignment
# gamm.data <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S1_mm_mt/GeneFull/filtered/")
# to read in a h5 file
#gamm.data<- Read10X_h5("/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/molecule_info.h5", use.names = TRUE, unique.features = TRUE)

# load dataset with different batches
gamm.data1 <- Read10X(data.dir = paste(dir.name,"output_S1-1_mm_mt/GeneFull/filtered/", sep = "", collapse = ""))
gamm.data2 <- Read10X(data.dir = paste(dir.name,"output_S1-2_mm_mt/GeneFull/filtered/", sep = "", collapse = ""))

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
gamm.data1.raw <- Read10X(data.dir = paste(dir.name,"output_S1-1_mm_mt/GeneFull/raw/", sep = "", collapse = ""))
gamm.data2.raw <- Read10X(data.dir = paste(dir.name,"output_S1-2_mm_mt/GeneFull/raw/", sep = "", collapse = ""))
# GetAssayData(object = gamm[["RNA"]], slot = "data")
# make soup channel and profile the soup with raw and filtered data
sc1 = SoupChannel(gamm.data1.raw, gamm.data1)
sc2 = SoupChannel(gamm.data2.raw, gamm.data2)
# get basic clusters
# make seurat object without filtering
gamm1 <- CreateSeuratObject(counts = gamm.data1, project = "gamm_s1-1_clusters")
gamm2 <- CreateSeuratObject(counts = gamm.data2, project = "gamm_s1-2_clusters")
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
DropletUtils:::write10xCounts("gamm_soupX_filtS1-1.mtx", out1)
DropletUtils:::write10xCounts("gamm_soupX_filtS1-2.mtx", out2)
# replace filtered count values with soupX values
# gamm[["RNA"]] <- SetAssayData(gamm[["RNA"]],
#                               slot = "data", 
#                               new.data = out)

## Doublet removal ## this should happen before filtering and normalization i think

# make a combined seurat object with adjusted data and filter
gamm <- CreateSeuratObject(counts = cbind(out1,out2), project = "gamm_s1", min.cells = 3, min.features = 200)
# make singular seurat objects to get accurate cell numbers for each lane
gamm1 <- CreateSeuratObject(counts = out1, project = "gamm_s1-1", min.cells = 3, min.features = 200)
gamm2 <- CreateSeuratObject(counts = out2, project = "gamm_s1-2", min.cells = 3, min.features = 200)
# done with soupx- remove raw data to save space
rm(gamm.data1.raw,gamm.data2.raw)

# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
# We use the set of all genes starting with MT- as a set of mitochondrial genes
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#gamm[["percent.mt"]] <- PercentageFeatureSet(gamm, pattern = "^MT-")
# no percent.mt in pig because they don't start with MT

# give specific features (genes) for pig data- 13 protein coding mt genes in pig:
mt.list <- c("ND1","ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L",
             "ND4", "ND5", "ND6", "CYTB")
# or ensemble ids
# mt.list2 <- c("ENSSSCG00000018065","ENSSSCG00000018069","ENSSSCG00000018075",
#              "ENSSSCG00000018078","ENSSSCG00000018080","ENSSSCG00000018081",
#              "ENSSSCG00000018082","ENSSSCG00000018084","ENSSSCG00000018086",
#              "ENSSSCG00000018087","ENSSSCG00000018091","ENSSSCG00000018092",
#              "ENSSSCG00000018094")
# check if all features are in the matrix
all(mt.list %in% rownames(gamm))
# if false, need to determine which ones
for (mt in mt.list) {print(mt)
                      print(all(mt %in% rownames(gamm)))
                      }
percent_mt <- PercentageFeatureSet(gamm, features = mt.list, assay = 'RNA')
gamm[["percent.mt"]] <- percent_mt
gamm[["percent.mt"]]
# for lanes 1 and 2
percent_mt1 <- PercentageFeatureSet(gamm1, features = mt.list, assay = 'RNA')
gamm1[["percent.mt"]] <- percent_mt1
percent_mt2 <- PercentageFeatureSet(gamm2, features = mt.list, assay = 'RNA')
gamm2[["percent.mt"]] <- percent_mt2
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts
# Visualize QC metrics as a violin plot
VlnPlot(gamm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# subset data -  filter out cells that have unique feature counts over 2,500 or less than 200 and > 5% mitochondrial counts
gamm <- subset(gamm, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
# subset data by percent- maybe this is better?
#counts <- GetAssayData(seurat_obj, slot="counts", assay="RNA")   
#genes.percent.expression <- rowMeans(counts>0 )*100   
#genes.filter <- names(gene.percent.expressed[gene.percent.expressed>1])  #select genes expressed in at least 1% of cells
#counts.sub <- counts[genes.filter,]
#new_seurat_object <- CreateSeuratObject(counts=counts.sub)
# getting numbers of filtered cells from lanes 1 & 2
gamm1 <- subset(gamm1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
sce1 <- as.SingleCellExperiment(gamm1)
lane1<- colData(sce1)@nrows
gamm2 <- subset(gamm2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
sce2 <- as.SingleCellExperiment(gamm2)
lane2<- colData(sce2)@nrows
# after you get lane 1 & 2 numbers, delete seurat object for them to save space
rm(gamm1,gamm2,sce1,sce2)
# replot
plot1 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gamm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
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


# Perform linear dimensional reduction on the variable features
gamm <- RunPCA(gamm, features = VariableFeatures(object = gamm))
# Examine and visualize PCA results a few different ways
print(gamm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gamm, dims = 1:2, reduction = "pca")
DimPlot(gamm, reduction = "pca")
#DimHeatmap(gamm, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(gamm, dims = 1:15, cells = 500, balanced = TRUE)

# Harmony
# define batches
gamm@meta.data$lanes <- c(rep("lane1", lane1), rep("lane2", lane2))
# check PCA plot
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gamm, reduction = "pca", pt.size = .1, group.by = "lanes")
p2 <- VlnPlot(object = gamm, features = "PC_1", group.by = "lanes", pt.size = .1)
plot_grid(p1,p2)
# run harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
gamm <- gamm %>% 
  RunHarmony("lanes", plot_convergence = TRUE)
# access harmony embeddings
harmony_embeddings <- Embeddings(gamm, 'harmony')
harmony_embeddings[1:5, 1:5]
# check out new pca plot
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gamm, reduction = "harmony", pt.size = .1, group.by = "lanes")
p2 <- VlnPlot(object = gamm, features = "harmony_1", group.by = "lanes", pt.size = .1)
plot_grid(p1,p2)

# Determine the ‘dimensionality’ of the dataset
# to overcome technical noise, Seurat clusters cells based on their PCA scores, 
# with each PC essentially representing a ‘metafeature’ that combines information 
# across a correlated feature set. The top principal components therefore represent 
# a robust compression of the dataset
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
# check harmony umap for differences across lanes
gamm <- RunUMAP(gamm, dims = 1:15, reduction = "harmony")
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(gamm, reduction = "umap", group.by = "lanes", pt.size = .1)#, split.by = 'lanes')
# note need to pip install leiden and pandas for FindClusters to run
#reticulate::py_install(packages ='leidenalg') #pip install leidenalg
#reticulate::py_install(packages ='pandas')
gamm <- FindClusters(gamm, resolution = 0.5, algorithm = "leiden")
# Look at cluster IDs of the first 5 cells
head(Idents(gamm), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# The goal of these algorithms is to learn the underlying manifold of the data in order 
# to place similar cells together in low-dimensional space. Cells within the graph-based 
# clusters determined above should co-localize on these dimension reduction plots. 
# As input to the UMAP and tSNE, we suggest using the same PCs as input to the 
# clustering analysis.
# If you haven't installed UMAP: 
# reticulate::py_install(packages ='umap-learn')
#gamm <- RunUMAP(gamm, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# plot individual clusters
DimPlot(gamm, reduction = "umap", label = T)
# save object
saveRDS(gamm, file = "GAMM_test1_lanes.rds")
# read back in the object
#gamm <- readRDS(file = "GAMM_test1.rds")


# Finding differentially expressed features (cluster biomarkers)

# # find all markers of cluster 1 compared to all other clusters
# cluster1.markers <- FindMarkers(gamm, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# #cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# #head(cluster5.markers, n = 5)
# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# gamm.markers <- FindAllMarkers(gamm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# gamm.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, order_by = avg_log2FC) %>%
#   print(n=28)
# # ROC test for cluster markers, 0= random, 1= perfect)
# cluster1.markers <- FindMarkers(gamm, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# head(cluster1.markers, n = 5)
# # visualize marker expression
# VlnPlot(gamm, features = c("SFRP2", "TRPM3"))
# # you can plot raw counts as well
# VlnPlot(gamm, features = c("SFRP2", "TRPM3"), slot = "counts", log = TRUE)
# # umap plot highlighting gene expression
# FeaturePlot(gamm, features = c("SFRP2", "TRPM3", "ENSSSCG00065027583", "SST"),label = TRUE)
# # DoHeatmap() generates an expression heatmap for given cells and features. 
# # In this case, we are plotting the top 20 markers
# gamm.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(gamm, features = top10$gene) + NoLegend()
# 
# # Assigning cell type identity to clusters
# 
# # Here canonical markers easily match the unbiased clustering to known cell types
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# # save object
# saveRDS(pbmc, file = "output/pbmc3k_final.rds")
# sessionInfo()
