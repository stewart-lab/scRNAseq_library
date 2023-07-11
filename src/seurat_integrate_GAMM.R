# load human data
# read in 10X
exp.matrix <- Read10X("ORG_D205_filtered_feature_bc_matrix/")
exp.matrix[,1:100]
# load in orthologs
orthologs <- read.csv2("pig-human_genesynbol_jack.txt_nodups.txt", sep="\t", header=TRUE)
# subset genes so they match only orthologous ones
human.orthos <- as.vector(orthologs$human.gene.name)
exp.matrix1<- exp.matrix[rownames(exp.matrix) %in% human.orthos, ]
# subset orthologs to they only match expr matrix
orthologs1<- orthologs[orthologs$human.gene.name %in% rownames(exp.matrix1), ]
# check length to make sure they are the same
length(rownames(exp.matrix1))
length(orthologs1$human.gene.name)
# switch to pig genes
# reorder based on exp matrix genes
ind_reorder <- match(rownames(exp.matrix1),orthologs1$human.gene.name)
orthologs.reorder<- orthologs1[ind_reorder,]
# now replace rownames of exp matrix with pig genes
pig.orthos <- as.vector(orthologs.reorder$Pig.gene.name)
row.names(exp.matrix1) <- pig.orthos
# check exp matrix row names
head(row.names(exp.matrix1))

# read in meta.data
metadata <- read.csv2("GSE142526_metadata/F6/cca_fetalvsorg_125CP_205_metadata.csv", sep=",", header=TRUE, row.names= 1)
# subset by D205 data
length(exp.matrix1@Dimnames[[2]])
metadata <- subset(metadata, metadata$orig.ident=="D205")
metadata<- as.data.frame(metadata)
# replace underscore with dash
row.names(metadata) <- sub("_", "-", row.names(metadata))
metadata.rows <-as.vector(row.names(metadata))

# subset counts data to only include cells in metadata
exp.matrix1<- exp.matrix1[, colnames(exp.matrix1) %in% metadata.rows]
dim(metadata)
dim(exp.matrix1)
# create seurat object
human_D205.seurat<- CreateSeuratObject(exp.matrix1, project = "human_D205", assay = "RNA",
                                       min.cells = 3, min.features = 200)
human_D205.seurat<-AddMetaData(object=human_D205.seurat, metadata=metadata)

# load gamm pig data sample 2
gamm.data1 <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-1_mm_mt/GeneFull/filtered/")
gamm.data2 <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-2_mm_mt/GeneFull/filtered/")

# subset pig data so that there is same genes
gamm.data1<- gamm.data1[rownames(gamm.data1) %in% pig.orthos, ]
gamm.data2<- gamm.data2[rownames(gamm.data2) %in% pig.orthos, ]
# both data sets must be pre-processed/ normalized the same way
# soupx
# scdblfinder

# make combined seurat objects for both lanes
gamm <- CreateSeuratObject(counts = cbind(gamm.data1,gamm.data2), project = "gamm_s2", min.cells = 3, min.features = 200)

# remove high mitochondral gene percentage
# give specific features (genes) for pig data- 13 protein coding mt genes in pig:
mt.list <- c("ND1","ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L",
             "ND4", "ND5", "ND6", "CYTB")
for (mt in mt.list) {print(mt)
  print(all(mt %in% rownames(gamm)))
}
for (mt in mt.list) {print(mt)
  print(all(mt %in% rownames(human_D205.seurat)))
}
# no mt genes here- need them in orthologs list
# get percent mt genes
percent_mt <- PercentageFeatureSet(gamm, features = mt.list, assay = 'RNA')
gamm[["percent.mt"]] <- percent_mt
percent_mt2 <- PercentageFeatureSet(human_D205.seurat, features = mt.list, assay = 'RNA')
human_D205.seurat[["percent.mt"]] <- percent_mt2

# plot features
plot1 <- FeatureScatter(human_D205.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_D205.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
# subset data -  filter out cells that have unique feature counts over 2,500 or less than 200 and > 5% mitochondrial counts
gamm <- subset(gamm, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
human_D205.seurat<- subset(human_D205.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
# use scTransform
gamm <- SCTransform(gamm)
human_D205.seurat <- SCTransform(human_D205.seurat)
# scran normalization didn't work- error with clustering
# prep for integration
ret.list<- list(gamm,human_D205.seurat)
features <- SelectIntegrationFeatures(object.list = ret.list, nfeatures = 3000)
ret.list <- PrepSCTIntegration(object.list = ret.list, anchor.features = features)

# # mapping and annotating query datasets
# reduce dim of reference
human_D205.seurat <- ScaleData(human_D205.seurat, verbose = FALSE)
human_D205.seurat <- RunPCA(human_D205.seurat, npcs = 30, verbose = FALSE)
human_D205.seurat <- RunUMAP(human_D205.seurat, reduction = "pca", dims = 1:30, verbose = FALSE)

p2 <- DimPlot(human_D205.seurat, reduction = "umap", group.by = "type", label = TRUE, repel = TRUE) +
  NoLegend()
p2
# get query data and find anchors
gamm.anchors <- FindTransferAnchors(reference = human_D205.seurat, query = gamm,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = gamm.anchors, refdata = human_D205.seurat$type,
                            dims = 1:30)
gamm <- AddMetaData(gamm, metadata = predictions)

# evaluate how well our predicted cell type annotations match the full reference. 
# need manual annotation data
#pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
#table(pancreas.query$prediction.match)
# examine some canonical cell type markers for specific pancreatic islet cell populations.
table(gamm$predicted.id)
# need marker genes
#VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")
# projection of a query onto the reference UMAP structure.
human_D205.seurat <- RunUMAP(human_D205.seurat, dims = 1:30, reduction = "pca", return.model = TRUE)
gamm <- MapQuery(anchorset = gamm.anchors, reference = human_D205.seurat, query = gamm,
                           refdata = list(celltype = "type"), reference.reduction = "pca", reduction.model = "umap")

# visualize
p1 <- DimPlot(human_D205.seurat, reduction = "umap", group.by = "type", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(gamm, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
