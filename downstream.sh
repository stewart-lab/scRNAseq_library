## QC covariates

# filter out genes only expressed in < 3 cells
# filter high mitochondrial expressed cells
# filter ambient RNA (SoupX)
# Doublet detection (DoubletFinder, scDblFinder: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html)
# Normalization (scran)
# Feature selection (remove genes with constant expressioin- seurat)
# dimensionality reduction (PCA, UMAP)

## cluster on preprocessed data

# cluster (seurat-Leiden algorithm)

## cluster annotation
