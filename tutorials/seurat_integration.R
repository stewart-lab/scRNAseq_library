# Seurat integration

# set up packages
devtools::install_github('satijalab/seurat-data')
library(Seurat)
library(SeuratData)
library(patchwork)

# get data
# show available data from SeuratData
AvailableData()
# set time out
getOption('timeout')
options(timeout=2000)
# install ifnb dataset
InstallData("ifnb")
# load dataset
LoadData("ifnb")

# preprocess
# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

# perform integration
# id anchors
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                         anchor.features = features)
# integrate
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# perform integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(immune.combined, reduction = "umap", split.by = "stim")

# ID conserved markers
# For performing differential expression after integration, we switch back to the original
# data
BiocManager::install('multtest')
install.packages('metap')
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)
# annotate with conserved markers
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                          "CCL2", "PPBP"), min.cutoff = "q9")
immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
DimPlot(immune.combined, label = TRUE)
# dotplot
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
                                                                      "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
                                                                      "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()

# Identify differential expressed genes across conditions
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

# what genes change under different conditions for cells of the same type?
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)

# visualize genes change in expression due to condition
FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3,
            cols = c("grey", "red"))
plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# integration on datasets normalized with SCTransform

# Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
# As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
# Run the PrepSCTIntegration() function prior to identifying anchors
# When running FindIntegrationAnchors(), and IntegrateData(), set the normalization.method parameter to the value SCT.
# When running sctransform-based workflows, including integration, do not run the ScaleData() function

LoadData("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
# find anchors and combine
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
# run PCA and UMAP
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
# visualize
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
              repel = TRUE)
p1 + p2

# mapping and annotating query datasets
library(Seurat)
library(SeuratData)
# dataset load
InstallData("panc8")

# First, we split the combined object into a list, with each dataset as an element 
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
# preprocess
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}
# identify ‘anchors’ between the individual datasets.
# here there are 3 pancreatic cell datasets
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
# We then pass these anchors to the IntegrateData() function, which returns a Seurat object.
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
# Seurat object will contain a new Assay with the integrated expression matrix. 
# Note that the original (uncorrected values) are still stored in the object in the “RNA” assay,

# scale, run PCA, and visulaize with umap the integrated matrix
library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2

# Cell type classification using an integrated reference
# project metadata onto query object
# In data transfer, Seurat does not correct or modify the query expression data.
# In data transfer, Seurat has an option (set by default) to project the PCA structure 
# of a reference onto the query, instead of learning a joint structure with CCA. 
# We generally suggest using this option when projecting data between scRNA-seq datasets.

# TransferData() function to classify the query cells based on reference data (a vector of reference cell type labels).
# TransferData() returns a matrix with predicted IDs and prediction scores, which we can add to the query metadata.
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype,
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

# evaluate how well our predicted cell type annotations match the full reference. 
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
# examine some canonical cell type markers for specific pancreatic islet cell populations.
table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")
# projection of a query onto the reference UMAP structure.
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

# visualize
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2



