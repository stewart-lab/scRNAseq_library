library(Seurat)
library(cowplot)
library(harmony)

load('data/pbmc_stim.RData')

# create one seurat object with both datasets
pbmc <- CreateSeuratObject(counts = cbind(stim.sparse, ctrl.sparse), project = "PBMC", min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)

# define data sets with variable stim
pbmc@meta.data$stim <- c(rep("STIM", ncol(stim.sparse)), rep("CTRL", ncol(ctrl.sparse)))

# check uncorrected PCs
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "stim", pt.size = .1)
plot_grid(p1,p2)

# run harmony with seurat object and specify vairable(s) to integrate out
options(repr.plot.height = 2.5, repr.plot.width = 6)
pbmc <- pbmc %>% 
  RunHarmony("stim", plot_convergence = TRUE)
# access harmony embeddings
harmony_embeddings <- Embeddings(pbmc, 'harmony')
harmony_embeddings[1:5, 1:5]
# check harmony plots
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "stim", pt.size = .1)
plot_grid(p1,p2)
# downstream analyess on harmony embeddings
pbmc <- pbmc %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
# check umap plot
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(pbmc, reduction = "umap", group.by = "stim", pt.size = .1, split.by = 'stim')
# id clusters
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = .1)
