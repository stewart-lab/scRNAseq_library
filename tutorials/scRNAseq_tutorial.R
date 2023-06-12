# install library
BiocManager::install("scRNAseq")
BiocManager::install("scater")
BiocManager::install("SingleR")
BiocManager::install("celldex")
# load data
library(scRNAseq)
sce.nest <- NestorowaHSCData()
# load annotation
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
anno <- select(ens.mm.v97, keys=rownames(sce.nest), 
               keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
rowData(sce.nest) <- anno[match(rownames(sce.nest), anno$GENEID),]
# check data
sce.nest

# QC
unfiltered <- sce.nest
# QC on spike-in only
library(scater)
stats <- perCellQCMetrics(sce.nest)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent")
sce.nest <- sce.nest[,!qc$discard]
# examine number of cells discarded:
colSums(as.matrix(qc))
# plot
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="altexps_ERCC_percent",
              colour_by="discard") + ggtitle("ERCC percent"),
  ncol=2
)

# Normalization
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce.nest)
sce.nest <- computeSumFactors(sce.nest, clusters=clusters)
sce.nest <- logNormCounts(sce.nest)
summary(sizeFactors(sce.nest))
# plot size factors to library size
plot(librarySizeFactors(sce.nest), sizeFactors(sce.nest), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")

# Variance modeling
# use spikeins to model technical noise
set.seed(00010101)
dec.nest <- modelGeneVarWithSpikes(sce.nest, "ERCC")
top.nest <- getTopHVGs(dec.nest, prop=0.1)
# plot
plot(dec.nest$mean, dec.nest$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec.nest)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
points(curfit$mean, curfit$var, col="red")

# Dimensionality reduction
set.seed(101010011)
sce.nest <- denoisePCA(sce.nest, technical=dec.nest, subset.row=top.nest)
sce.nest <- runTSNE(sce.nest, dimred="PCA")
# check PC number
ncol(reducedDim(sce.nest, "PCA"))

# Clustering
snn.gr <- buildSNNGraph(sce.nest, use.dimred="PCA")
colLabels(sce.nest) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(colLabels(sce.nest))
plotTSNE(sce.nest, colour_by="label")

# Marker detection
# get DE across clusters
markers <- findMarkers(sce.nest, colLabels(sce.nest), 
                       test.type="wilcox", direction="up", lfc=0.5,
                       row.data=rowData(sce.nest)[,"SYMBOL",drop=FALSE])
# marker genes for cluster 8
chosen <- markers[['8']]
best <- chosen[chosen$Top <= 10,]
aucs <- getMarkerEffects(best, prefix="AUC")
rownames(aucs) <- best$SYMBOL

library(pheatmap)
pheatmap(aucs, color=viridis::plasma(100))

# Cell type annotation with refernce
library(SingleR)
mm.ref <- celldex::MouseRNAseqData()

# Renaming to symbols to match with reference row names.
renamed <- sce.nest
rownames(renamed) <- uniquifyFeatureNames(rownames(renamed),
                                          rowData(sce.nest)$SYMBOL)
labels <- SingleR(renamed, mm.ref, labels=mm.ref$label.fine)
# most clusters not assigned to any single lineage
tab <- table(labels$labels, colLabels(sce.nest))
pheatmap(log10(tab+10), color=viridis::viridis(100))
