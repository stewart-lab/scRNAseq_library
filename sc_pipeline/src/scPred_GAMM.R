# load packages
library("scPred")
library("Seurat")
library("magrittr")
library("harmony")

# get reference and query data
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
head(row.names(exp.matrix1))
row.names(exp.matrix1) <- pig.orthos
# check exp matrix row names
head(row.names(exp.matrix1))

# read in meta.data
metadata <- read.csv2("GSE142526_metadata/F6/cca_fetalvsorg_125CP_205_metadata.csv", sep=",", header=TRUE, row.names= 1)
# subset by D205 data
length(exp.matrix1@Dimnames[[2]])
metadata <- subset(metadata, metadata$orig.ident=="D205")
metadata<- as.data.frame(metadata)
complete.cases(metadata)
metadata<-na.omit(metadata)

# get numbers for each celltype
table(metadata$type)
# drop AC1 because only 2
metadata <- subset(x = metadata, type != 'AC1')
# replace underscore with dash
row.names(metadata) <- sub("_", "-", row.names(metadata))
metadata.rows <-as.vector(row.names(metadata))
dim(metadata)

# subset counts data to only include cells in metadata
exp.matrix1<- exp.matrix1[, colnames(exp.matrix1) %in% metadata.rows]
dim(exp.matrix1)
exp.matrix1 = na.omit(exp.matrix1)
dim(exp.matrix1)

# load gamm pig data sample 2
gamm.data1 <- Read10X(data.dir = "/Users/bmoore/Desktop/GitHub/scRNAseq/GAMM/output_S2-1_mm_mt/GeneFull/filtered/")

# subset pig data so that there is same genes
dim(gamm.data1)
gamm.data1<- gamm.data1[rownames(gamm.data1) %in% pig.orthos, ]
dim(gamm.data1)

# subset human refernce based on gamm data
exp.matrix1<- exp.matrix1[rownames(exp.matrix1) %in% rownames(gamm.data1), ]
dim(exp.matrix1)
# create seurat object
human_D205.seurat<- CreateSeuratObject(exp.matrix1, project = "human_D205", assay = "RNA",
                                       min.cells = 3, min.features = 200)
human_D205.seurat<-AddMetaData(object=human_D205.seurat, metadata=metadata)


# both data sets must be pre-processed/ normalized the same way
# soupx
# scdblfinder

# make combined seurat object
gamm <- CreateSeuratObject(counts = gamm.data1, project = "gamm_s2", min.cells = 3, min.features = 200)
dim(gamm)
# run preprocessing on reference
#Single SCTransform command replaces NormalizeData, ScaleData, and FindVariableFeatures. 
# We will also correct for % MT genes and cell cycle scores using vars.to.regress
human_D205.seurat <- human_D205.seurat %>% 
  #SCTransform() %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
# plot the umap for cell type
# we don't do clustering because clusters already annotated
DimPlot(human_D205.seurat, group.by = "type", label = TRUE, repel = TRUE)

# training the classifyer

# getFeatureSpace will create a scPred object stored in the @misc slot. This 
# object will contained all required information to classify cells.
human_D205.seurat <- getFeatureSpace(human_D205.seurat, "type")

# Downsample the number of cells per identity class, since min is 17, set at 17
# human_D205.seurat.ds<- subset(x = human_D205.seurat, downsample = 17)

# train models for each cell type
# cv has to be set BELOW min class
human_D205.seurat <- trainModel(human_D205.seurat,resampleMethod = "cv",
                                number = 10, seed = 704)

# parameters
  # method: resampling method
  # p: For leave-group out cross-validation: the training percentage
  # search = "grid", "random"
  # sampling="down", "up" addtl sampling to resolve class imbalance
# if you get this error:
# There were missing values in resampled performance measures.
# try: 1. make sure no NAs in data 2. increasing the cvs, 3. down or up sample to make classes even

# get training prob for each cell-cell type
get_probabilities(human_D205.seurat) %>% head()
# use get_scpred method to retrieve the scPred object from the Seurat object.
# this gives stats on each cell type prediction
get_scpred(human_D205.seurat)
# visualize cell type probabilties
plot_probabilities(human_D205.seurat)
# can try other models from caret:
# https://topepo.github.io/caret/available-models.html
# pass in the model parameter
# can reclassify a subset of cells (like those that didn't work well)
human_D205.seurat <- trainModel(human_D205.seurat, model = "rf", 
                                   resampleMethod = "cv",
                                   number = 10, seed = 704,
                                   reclassify = c("AC2", "T1/T3","T2"))
get_scpred(human_D205.seurat)
plot_probabilities(human_D205.seurat)
human_D205.seurat <- trainModel(human_D205.seurat, model = "glm", #'logreg'
                                   resampleMethod = "cv",
                                   number = 10, seed = 704,
                                   reclassify = c("AC2", "T1/T3","T2"))
get_scpred(human_D205.seurat)
plot_probabilities(human_D205.seurat)

# classify cells

# An important requirement for classifying cells is using the same normalization 
# method for both the reference and the query datasets.
# normalize query
gamm <- NormalizeData(gamm)

# scPred now uses Harmony to align the query data onto the training low-dimensional 
# space used as reference. Once the data is aligned, cells are classified using the 
# pre-trained models.
# predict query
gamm <- scPredict(gamm, human_D205.seurat)

#plot the classifications over the aligned data.
DimPlot(gamm, group.by = "scpred_prediction", reduction = "scpred")
# run UMAP using the aligned data as an input
gamm <- RunUMAP(gamm, reduction = "scpred", dims = 1:30)
# plot the predicted labels for each cell type over the UMAP:
DimPlot(gamm, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
#  probabilities of each cell in the @meta.data slot of the query Seurat object.
# get labels
labels<- as.vector(unique(colnames(gamm@meta.data)))
labels<-labels[! labels %in% c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'scpred_max', 'scpred_prediction', 'scpred_no_rejection')]
# visualize the probabilities over the UMAP plot:
FeaturePlot(gamm, labels)
# verify performance with manual annotation
# verify model performance in query data
##crossTab(query, "cell_type", "scpred_prediction")
# check out proportion of cells
##crossTab(query, "cell_type", "scpred_prediction", output = "prop")

# access the classifiers
get_classifiers(human_D205.seurat)
# Each model can be normally treated using the caret enviroment.
caret::plot.train(get_classifiers(reference)[["NK cell"]])