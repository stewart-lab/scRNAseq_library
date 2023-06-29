# scPred
# install scPred
devtools::install_github("powellgenomicslab/scPred")

# load packages
library("scPred")
library("Seurat")
library("magrittr")
library("harmony")

# use Peripheral Mononuclear Blood Cells (PBMC) data 
# these are seurat objects already
reference <- scPred::pbmc_1
query <- scPred::pbmc_2

# run preprocessing on reference
reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
# plot the umap for cell type
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

# training the classifyer

# getFeatureSpace will create a scPred object stored in the @misc slot. This 
# object will contained all required information to classify cells.
reference <- getFeatureSpace(reference, "cell_type")
# train models for each cell type
reference <- trainModel(reference)
# get training prob for each cell-cell type
get_probabilities(reference) %>% head()
# use get_scpred method to retrieve the scPred object from the Seurat object.
# this gives stats on each cell type prediction
get_scpred(reference)
# visualize cell type probabilties
plot_probabilities(reference)
# can try other models from caret:
# https://topepo.github.io/caret/available-models.html
# pass in the model parameter
# can reclassify a subset of cells (like those that didn't work well)
reference <- trainModel(reference, model = "mda", reclassify = c("cMono", "ncMono"))
get_scpred(reference)
plot_probabilities(reference)
# classify cells

# An important requirement for classifying cells is using the same normalization 
# method for both the reference and the query datasets.
# normalize query
query <- NormalizeData(query)

# scPred now uses Harmony to align the query data onto the training low-dimensional 
# space used as reference. Once the data is aligned, cells are classified using the 
# pre-trained models.
# predict query
query <- scPredict(query, reference)

# scPred will store the final classifications in the scpred_prediction column of the Seurat meta data. 
# It will store the aligned data as a scpred reduction.
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
# run umap using aligned data as input
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
# plot the predicted labels for each cell type over the UMAP:
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
# compare with original labels
DimPlot(query, group.by = "cell_type", label = TRUE, repel = TRUE)
# scPred stores the probabilities of each cell in the @meta.data slot of the query Seurat object.
# visualize cell type probabilities on umap
FeaturePlot(query, c("scpred_B.cell", "scpred_CD4.T.cell", "scpred_CD8.T.cell", 
                     "scpred_cMono", "scpred_ncMono", "scpred_Plasma.cell", 
                     "scpred_cDC", "scpred_pDC"))
# verify model performance in query data
crossTab(query, "cell_type", "scpred_prediction")
# check out proportion of cells
crossTab(query, "cell_type", "scpred_prediction", output = "prop")

# access the classifiers
get_classifiers(reference)
# Each model can be normally treated using the caret enviroment.
caret::plot.train(get_classifiers(reference)[["NK cell"]])

# use a different prediction model
# logistic regression via glm()
reference <- trainModel(reference, model = "glm")
# don't realign
# Training and alignning the data are separate processes. Therefore, if a query 
# dataset has already being aligned to a reference data via scPred/harmony and the 
# prediction models have changed, then we can use recompute_alignment = FALSE to 
# avoid aligning step
query <- scPredict(query, reference, recompute_alignment = FALSE)
# use a different probability threshold
# default is 0.55
query <- scPredict(query, reference, recompute_alignment = FALSE, threshold = 0.9)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)

# parallel training
#  resampling performed for each model can be parallelized via doParallel
library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
reference <- trainModel(reference, model = "mda", allowParallel = TRUE)

stopCluster(cl)

devtools::session_info()
