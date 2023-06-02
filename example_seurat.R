library(dplyr)
library(Seurat)
BiocManager::install("SeuratData")
library(SeuratData)
library(patchwork)
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(scran)
library(BiocParallel)

# Load the PBMC dataset

pbmc.data <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/seurat/data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# sparse matrix with reduced memory
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size

# QC and filtering

# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination

# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
# We use the set of all genes starting with MT- as a set of mitochondrial genes
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Doublet removal ## this should happen before normalization i think
seu <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# ambient RNA removal- SoupX


# Normalize the data

# default is “LogNormalize” that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), 
# and log-transforms the result. 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# scran normalization
library(scran)
sce <- as.SingleCellExperiment(pbmc)
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
pbmc[["RNA"]] <- SetAssayData(pbmc[["RNA"]],
                            slot = "data", 
                            new.data = logcounts(sce))

pbmc$sizeFactors <- sizeFactors(sce)
UpdateSeuratObject(pbmc)
# gene expression with reg counts
VlnPlot(pbmc, "CD3E", slot = "counts")
# gene expression with scran norm counts
VlnPlot(pbmc, "CD3E", slot = "data")

# ID highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# feature selection with scry
BiocManager::install('scry')
library(scry)
sce2 <- as.SingleCellExperiment(pbmc)
sce2 <- devianceFeatureSelection(sce2, assay = "counts")
sce2
rowData(sce2)$binomial_deviance
# add biological deviance ## not sure how to do this
pbmc$deviance<- rowData(sce2)$binomial_deviance

# Scale the data
# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
#   
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Remove unwanted variation
# also use the ScaleData() function to remove unwanted sources of variation from a 
# single-cell dataset. For example, we could ‘regress out’ heterogeneity associated with 
# (for example) cell cycle stage, or mitochondrial contamination.
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# new normalization workflow, SCTransform()

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# to overcome technical noise, Seurat clusters cells based on their PCA scores, 
# with each PC essentially representing a ‘metafeature’ that combines information 
# across a correlated feature set. The top principal components therefore represent 
# a robust compression of the dataset
# JackStraw randomly permutes a subset (1%) of the data and reruns PCA, constructing a null
# distribution of feature scores. Significant PCs have strong enrichment of low p-val features

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# The JackStrawPlot() function provides a visualization tool for comparing the 
# distribution of p-values for each PC with a uniform distribution (dashed line). 
# ‘Significant’ PCs will show a strong enrichment of features with low p-values 
#(solid curve above the dashed line). In this case it appears that there is a sharp 
# drop-off in significance after the first 10-12 PCs.
JackStrawPlot(pbmc, dims = 1:15)
#  alternative heuristic method generates an ‘Elbow plot’: a ranking of principle 
# components based on the percentage of variance explained by each one (ElbowPlot() 
# function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that 
# the majority of true signal is captured in the first 10 PCs.
ElbowPlot(pbmc)

# Cluster the Cells using only genes from chosen PCs
# methods embed cells in a graph structure - for example a K-nearest neighbor 
# (KNN) graph, with edges drawn between cells with similar feature expression patterns,
# and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or 
# ‘communities’.
# To cluster the cells, we next apply modularity optimization techniques such as 
# the Louvain algorithm (default) or SLM-- need to use Leiden here, to iteratively group cells together, with 
# the goal of optimizing the standard modularity function. The FindClusters() 
# function implements this procedure, and contains a resolution parameter that sets 
# the ‘granularity’ of the downstream clustering, with increased values leading 
# to a greater number of clusters. We find that setting this parameter between 
# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
# Optimal resolution often increases for larger datasets. The clusters can be 
# found using the Idents() function
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# note need to pip install leiden and pandas for this to run
# reticulate::py_install(packages ='leiden')
# reticulate::py_install(packages ='pandas')
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = "leiden")
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# The goal of these algorithms is to learn the underlying manifold of the data in order 
# to place similar cells together in low-dimensional space. Cells within the graph-based 
# clusters determined above should co-localize on these dimension reduction plots. 
# As input to the UMAP and tSNE, we suggest using the same PCs as input to the 
# clustering analysis.
# If you haven't installed UMAP: 
# reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
# save object
saveRDS(pbmc, file = "output/pbmc_tutorial.rds")
# read back in the object
pbmc <- readRDS(file = "output/pbmc_tutorial.rds")


# Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
# ROC test for cluster markers, 0= random, 1= perfect)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)
# visualize marker expression
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# umap plot highlighting gene expression
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters

# Here canonical markers easily match the unbiased clustering to known cell types
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# save object
saveRDS(pbmc, file = "output/pbmc3k_final.rds")
sessionInfo()
