# trajectory anlaysis
library(scRNAseq)
# install packages
BiocManager::install("TSCAN")
BiocManager::install("slingshot")
BiocManager::install("tradeSeq")
BiocManager::install("velociraptor",force = TRUE)
BiocManager::install("zellkonverter")
BiocManager::install("anndata")
# see scRNAseq_tutorial.R for pre-processing, clustering, annotation
#sce.nest <- NestorowaHSCData()
sce.nest

# Cluster-based minimum spanning tree
# TSCAN uses the clustering to summarize the data into a smaller set of discrete 
# units, computes cluster centroids by averaging the coordinates of its member 
# cells, and then forms the minimum spanning tree (MST) across those centroids. 
library(scater)
by.cluster <- aggregateAcrossCells(sce.nest, ids=colLabels(sce.nest))
centroids <- reducedDim(by.cluster, "PCA")

# Set clusters=NULL as we have already aggregated above.
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst
# plot in tSNE space
line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")

plotTSNE(sce.nest, colour_by="label") + 
  geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge))
# obtain a pseudotime ordering by projecting the cells onto the MST with mapCellsToEdges(). 
# We move each cell onto the closest edge of the MST; the pseudotime is then calculated 
# as the distance along the MST to this new position from a “root node” with orderCells(). 
# For our purposes, we will arbitrarily pick one of the endpoint nodes as the root
map.tscan <- mapCellsToEdges(sce.nest, mst=mst, use.dimred="PCA")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)
# multiple sets of pseudotimes are reported for a branched trajectory
# visualize
common.pseudo <- averagePseudotime(tscan.pseudo) 
plotTSNE(sce.nest, colour_by=I(common.pseudo), 
         text_by="label", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge))
# quickPseudotime wrapper
pseudo.all <- quickPseudotime(sce.nest, use.dimred="PCA")
head(pseudo.all$ordering)
# MST with outgroup
pseudo.og <- quickPseudotime(sce.nest, use.dimred="PCA", outgroup=TRUE)
set.seed(10101)
plot(pseudo.og$mst)
# distance based on mutual nearest neighbor- different from shortest distance between centroids
pseudo.mnn <- quickPseudotime(sce.nest, use.dimred="PCA", dist.method='mnn')
mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)
plotTSNE(sce.nest, colour_by=I(mnn.pseudo), text_by="label", text_colour="red") +
  geom_line(data=pseudo.mnn$connected$TSNE, mapping=aes(x=TSNE1, y=TSNE2, group=edge))

# Principle curves
# a non-linear generalization of PCA where the axes of most variation are allowed to bend
library(slingshot)
sce.sling <- slingshot(sce.nest, reducedDim='PCA')
head(sce.sling$slingPseudotime_1)
embedded <- embedCurves(sce.sling, "TSNE")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])

plotTSNE(sce.sling, colour_by="slingPseudotime_1") +
  geom_path(data=embedded, aes(x=TSNE1, y=TSNE2), size=1.2)
# fit to cluster to accomadate bifurcating data
sce.sling2 <- slingshot(sce.nest, cluster=colLabels(sce.nest), reducedDim='PCA')
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)
sce.nest <- runUMAP(sce.nest, dimred="PCA")
reducedDim(sce.sling2, "UMAP") <- reducedDim(sce.nest, "UMAP")

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
}

gg
# We can use slingshotBranchID() to determine whether a particular cell is shared 
# across multiple curves or is unique to a subset of curves (i.e., is located “after” branching).
curve.assignments <- slingBranchID(sce.sling2)
table(curve.assignments)
# we can speed up the algorithm by approximating each principal curve with a fixed number of points.
sce.sling3 <- slingshot(sce.nest, cluster=colLabels(sce.nest), 
                        reducedDim='PCA', approx_points=100)
pseudo.paths3 <- slingPseudotime(sce.sling3)
head(pseudo.paths3)
# MST constructed with an OMEGA cluster to avoid connecting unrelated trajectoires
sce.sling4 <- slingshot(sce.nest, cluster=colLabels(sce.nest), 
                        reducedDim='PCA', approx_points=100, omega=TRUE)
pseudo.paths4 <- slingPseudotime(sce.sling4)
head(pseudo.paths4)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
gg <- plotUMAP(sce.sling4, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling4, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
}
gg
# Changes along a trajectory
library(TSCAN)
pseudo <- testPseudotime(sce.nest, pseudotime=tscan.pseudo[,1])[[1]]
pseudo$SYMBOL <- rowData(sce.nest)$SYMBOL
pseudo[order(pseudo$p.value),]
# filter out cluster 7
# Making a copy of our SCE and including the pseudotimes in the colData.
sce.nest2 <- sce.nest
sce.nest2$TSCAN.first <- pathStat(tscan.pseudo)[,1]
sce.nest2$TSCAN.second <- pathStat(tscan.pseudo)[,2]

# Discarding the offending cluster.
discard <- "7"
keep <- colLabels(sce.nest)!=discard
sce.nest2 <- sce.nest2[,keep]

# Testing against the first path again.
pseudo <- testPseudotime(sce.nest2, pseudotime=sce.nest2$TSCAN.first)
pseudo$SYMBOL <- rowData(sce.nest2)$SYMBOL
sorted <- pseudo[order(pseudo$p.value),]
# examine top down regulated genes
up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)
best <- head(up.left$SYMBOL, 10)
plotExpression(sce.nest2, features=best, swap_rownames="SYMBOL",
               x="TSCAN.first", colour_by="label")
# examine top upregulated genes
up.right <- sorted[sorted$logFC > 0,]
head(up.right, 10)
best <- head(up.right$SYMBOL, 10)
plotExpression(sce.nest2, features=best, swap_rownames="SYMBOL",
               x="TSCAN.first", colour_by="label")
# overal heatmap
on.first.path <- !is.na(sce.nest2$TSCAN.first)
plotHeatmap(sce.nest2[,on.first.path], order_columns_by="TSCAN.first", 
            colour_columns_by="label", features=head(up.right$SYMBOL, 50),
            center=TRUE, swap_rownames="SYMBOL")

# detect branching points with DE genes
starter <- "3"
tscan.pseudo2 <- orderCells(map.tscan, mst, start=starter)
# Making a copy and giving the paths more friendly names.
sub.nest <- sce.nest
sub.nest$TSCAN.first <- pathStat(tscan.pseudo2)[,1]
sub.nest$TSCAN.second <- pathStat(tscan.pseudo2)[,2]
sub.nest$TSCAN.third <- pathStat(tscan.pseudo2)[,3]

# Subsetting to the desired cluster containing the branch point.
keep <- colLabels(sce.nest) == starter
sub.nest <- sub.nest[,keep]

# Showing only the lines to/from our cluster of interest.
line.data.sub <- line.data[grepl("^3--", line.data$edge) | grepl("--3$", line.data$edge),]
ggline <- geom_line(data=line.data.sub, mapping=aes(x=TSNE1, y=TSNE2, group=edge))

gridExtra::grid.arrange(
  plotTSNE(sub.nest, colour_by="TSCAN.first") + ggline,
  plotTSNE(sub.nest, colour_by="TSCAN.second") + ggline,
  plotTSNE(sub.nest, colour_by="TSCAN.third") + ggline,
  ncol=3
)
#  apply testPseudotime() to each path involving cluster 3
pseudo1 <- testPseudotime(sub.nest, df=1, pseudotime=sub.nest$TSCAN.first)
pseudo1$SYMBOL <- rowData(sce.nest)$SYMBOL
pseudo1[order(pseudo1$p.value),]
pseudo2 <- testPseudotime(sub.nest, df=1, pseudotime=sub.nest$TSCAN.second)
pseudo2$SYMBOL <- rowData(sce.nest)$SYMBOL
pseudo2[order(pseudo2$p.value),]
pseudo3 <- testPseudotime(sub.nest, df=1, pseudotime=sub.nest$TSCAN.third)
pseudo3$SYMBOL <- rowData(sce.nest)$SYMBOL
pseudo3[order(pseudo3$p.value),]
# find genes significnat in one path but not in others
only3 <- pseudo3[which(pseudo3$FDR <= 0.05 & 
                         (pseudo2$p.value >= 0.05 | sign(pseudo1$logFC)!=sign(pseudo3$logFC)) &
                         (pseudo2$p.value >= 0.05 | sign(pseudo2$logFC)!=sign(pseudo3$logFC))),]
only3[order(only3$p.value),]
# upregulation of interesting genes such as Gata2, Cd9 and Apoe in this path, along with downregulation of Flt3
gridExtra::grid.arrange(
  plotTSNE(sub.nest, colour_by="Flt3", swap_rownames="SYMBOL") + ggline,
  plotTSNE(sub.nest, colour_by="Apoe", swap_rownames="SYMBOL") + ggline,
  plotTSNE(sub.nest, colour_by="Gata2", swap_rownames="SYMBOL") + ggline,
  plotTSNE(sub.nest, colour_by="Cd9", swap_rownames="SYMBOL") + ggline
)
# Generalized additive models (GAMs) are quite popular for pseudotime-based DE analyses as they are able to handle non-normal noise distributions and a greater diversity of non-linear trends.
# if we assume that our pseudotime values are comparable across paths of the MST, we can use the patternTest() function to test for significant differences in expression between paths
# Getting rid of the NA's; using the cell weights
# to indicate which cell belongs on which path.
nonna.pseudo <- pathStat(tscan.pseudo)
not.on.path <- is.na(nonna.pseudo)
nonna.pseudo[not.on.path] <- 0
cell.weights <- !not.on.path
storage.mode(cell.weights) <- "numeric"

# Fitting a GAM on the subset of genes for speed.
library(tradeSeq)
fit <- fitGAM(counts(sce.nest)[1:100,], 
              pseudotime=nonna.pseudo,
              cellWeights=cell.weights)

res <- patternTest(fit)
res$Symbol <- rowData(sce.nest)[1:100,"SYMBOL"]
res <- res[order(res$pvalue),]
head(res, 10)

# Finding the root

# the root of the trajectory is best set to the “start” of the differentiation process,
# i.e., the most undifferentiated state that is observed in the dataset. It is usually
# possible to identify this state based on the genes that are expressed at each point
# of the trajectory. However, when such prior biological knowledge is not available,
# we can fall back to the more general concept that undifferentiated cells have more
# diverse expression profiles.

## We quantify the diversity of expression by computing the entropy of each cell’s 
## expression profile, with higher entropies representing greater diversity. 

library(TSCAN)
#sce <- scuttle::mockSCE()
#ent <- perCellEntropy(sce)
# getting Error in .read_block_OLD(x, viewport, as.sparse = as.sparse) : 
# could not find function ".read_block_OLD"
entropy <- perCellEntropy(sce.nest)
ent.data <- data.frame(cluster=colLabels(sce.nest), entropy=entropy)
ggplot(ent.data, aes(x=cluster, y=entropy)) + 
  geom_violin() +
  coord_cartesian(ylim=c(7, NA)) +
  stat_summary(fun=median, geom="point")

# RNA velocity to ID root
## given gene, a high ratio of unspliced to spliced transcripts indicates that 
# that gene is being actively upregulated, under the assumption that the increase 
# in transcription exceeds the capability of the splicing machinery to process the 
# pre-mRNA. Conversely, a low ratio indicates that the gene is being downregulated 
# as the rate of production and processing of pre-mRNAs cannot compensate for the 
# degradation of mature transcripts. Thus, we can infer that cells with high and 
# low ratios are moving towards a high- and low-expression state, respectively, 
# allowing us to assign directionality to any trajectory or even individual cells.

library(scRNAseq)
sce.sperm <- HermannSpermatogenesisData(strip=TRUE, location=TRUE)
assayNames(sce.sperm)

# Quality control:
library(scuttle)
is.mito <- which(seqnames(sce.sperm)=="MT")
sce.sperm <- addPerCellQC(sce.sperm, subsets=list(Mt=is.mito), assay.type="spliced")
qc <- quickPerCellQC(colData(sce.sperm), sub.fields=TRUE)
sce.sperm <- sce.sperm[,!qc$discard]

# Normalization:
set.seed(10000)
library(scran)
sce.sperm <- logNormCounts(sce.sperm, assay.type="spliced")
dec <- modelGeneVarByPoisson(sce.sperm, assay.type="spliced")
hvgs <- getTopHVGs(dec, n=2500)

# Dimensionality reduction:
set.seed(1000101)
library(scater)
sce.sperm <- runPCA(sce.sperm, ncomponents=25, subset_row=hvgs)
sce.sperm <- runTSNE(sce.sperm, dimred="PCA")

# write as anndata object
library(zellkonverter)
library(anndata)
adata<-SCE2AnnData(sce.sperm)
write_h5ad(adata, "sce.sperm.h5ad")
# velociraptor performs calculations cia scvelo
library(velociraptor) # cannot install python packages
velo.out <- scvelo(sce.sperm, assay.X="spliced", 
                   subset.row=hvgs, use.dimred="PCA")
velo.out
