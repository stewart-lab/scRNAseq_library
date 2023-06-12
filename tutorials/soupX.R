# check out https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html


library(SoupX)
library(Seurat)
library(Matrix)

# get raw and filtered data (PBMC dataset)


tmpDir = tempdir(check = TRUE)
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "tod.tar.gz"))
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "toc.tar.gz"))
untar(file.path(tmpDir, "tod.tar.gz"), exdir = tmpDir)
untar(file.path(tmpDir, "toc.tar.gz"), exdir = tmpDir)

# load data into 10x Soup channel object
sc = load10X(tmpDir)
# or manually create soup channel
toc = Seurat::Read10X(file.path(tmpDir, "filtered_gene_bc_matrices", "GRCh38"))
tod = Seurat::Read10X(file.path(tmpDir, "raw_gene_bc_matrices", "GRCh38"))
sc = SoupChannel(tod, toc)
# pre-loaded and processed object
data(PBMC_sc)
sc = PBMC_sc
sc
# profile the soup
sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc = estimateSoup(sc)
# to change estimateSoup parameters:
library(Matrix)
toc = sc$toc
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
# Calculate soup profile
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops, soupProf)
# Add extra meta data to the SoupChannel object (clusters)- not necessary but helps estimate soup
data(PBMC_metaData)
sc = setClusters(sc, setNames(PBMC_metaData$Cluster, rownames(PBMC_metaData)))
# visualize with dimension reduction
sc = setDR(sc, PBMC_metaData[colnames(sc$toc), c("RD1", "RD2")])
library(ggplot2)
dd = PBMC_metaData[colnames(sc$toc), ]
mids = aggregate(cbind(RD1, RD2) ~ Annotation, data = dd, FUN = mean)
gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = Annotation), size = 0.2) + 
  geom_label(data = mids, aes(label = Annotation)) + ggtitle("PBMC 4k Annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

# Estimate contamination fraction
sc  = autoEstCont(sc)
# Infer corrected table of counts and rount to integer
out = adjustCounts(sc, roundToInt = TRUE)
# write
DropletUtils:::write10xCounts("./strainedCounts", out)
# load into seurat
library(Seurat)
srat = CreateSeuratObject(out)