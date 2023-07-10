# singleCellNet

install.packages("devtools")
devtools::install_github("pcahan1/singleCellNet")
library(singleCellNet)

# object conversion
# integrate single cell exp
#exp_type options can be: counts, data, and scale.data if they are available in your sce object
scefile <- extractSCE(sce_object, exp_slot_name = "counts") 
sampTab = scefile$sampTab
expDat = scefile$expDat
# integrate seurat object
#exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
seuratfile <- extractSeurat(seurat_object, exp_slot_name = "counts")
sampTab = seuratfile$sampTab
expDat = seuratfile$expDat

# load query data
stPark = utils_loadObject("sampTab_Park_MouseKidney_062118.rda")
expPark = utils_loadObject("expDat_Park_MouseKidney_062218.rda")
dim(expPark)

genesPark = rownames(expPark)

rm(expPark)
gc()
# Load the training data
expTMraw = utils_loadObject("expMatrix_TM_Raw_Oct_12_2018.rda")
dim(expTMraw)

# meta data
stTM = utils_loadObject("sampTab_TM_053018.rda")
dim(stTM)

stTM<-droplevels(stTM)

# Find genes in common to the data sets and limit analysis to these
commonGenes = intersect(rownames(expTMraw), genesPark)
length(commonGenes)

expTMraw = expTMraw[commonGenes,]

# Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=stTM, ncells=100, dLevel="newAnn")
stTrain = stList[[1]]
expTrain = expTMraw[,rownames(stTrain)]

# Train the classifier
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, 
                                  nTopGenes = 10, nRand = 70, nTrees = 1000, 
                                  nTopGenePairs = 25, dLevel = "newAnn", 
                                  colName_samp = "cell"))

# Apply to held out data
#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="newAnn") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = expTMraw[commonGenes,rownames(stTest)]

#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

# Assess classifier
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "newAnn", classQuery = "newAnn", nRand = 50)

plot_PRs(tm_heldoutassessment)
plot_metrics(tm_heldoutassessment)

# Classification result heatmap

#Create a name vector label used later in classification heatmap where the values 
#are cell types/ clusters and names are the sample names

nrand = 50
sla = as.vector(stTest$newAnn)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created

sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)

# Attribution plot
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="newAnn", sid="cell")

#Visualize average top pairs genes expression for training data
gpTab = compareGenePairs(query_exp = expTest, training_exp = expTrain, training_st = stTrain, classCol = "newAnn", sampleCol = "cell", RF_classifier = class_info$cnProc$classifier, numPairs = 20, trainingOnly= TRUE)

train = findAvgLabel(gpTab = gpTab, stTrain = stTrain, dLevel = "newAnn")

hm_gpa_sel(gpTab, genes = class_info$cnProc$xpairs, grps = train, maxPerGrp = 50)

# Apply to query data
expPark = utils_loadObject("expDat_Park_MouseKidney_062218.rda") 

nqRand = 50
system.time(crParkall<-scn_predict(class_info[['cnProc']], expPark, nrand=nqRand))

# visualize
sgrp = as.vector(stPark$description1)
names(sgrp) = as.vector(stPark$sample_name)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

# heatmap classification result
sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)

# Classification annotation assignment
# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.

stPark <- get_cate(classRes = crParkall, sampTab = stPark, dLevel = "description1", sid = "sample_name", nrand = nqRand)

# violin plot
sc_violinClass(sampTab = stPark, classRes = crParkall, sid = "sample_name", dLevel = "description1", addRand = nqRand)

# Skyline plot of classification results
library(viridis)
stKid2 = addRandToSampTab(crParkall, stPark, "description1", "sample_name")
skylineClass(crParkall, "T cell", stKid2, "description1",.25, "sample_name")
# convert to seurat object
library(Seurat)
stPark.obj <- CreateSeuratObject(counts= expPark, meta.data=stPark)

# Cross-species classification
# Load the mouse training and human query data
stQuery = utils_loadObject("stDat_beads_mar22.rda")
expQuery = utils_loadObject("6k_beadpurfied_raw.rda") # use Matrix if RAM low
dim(expQuery)

stTM = utils_loadObject("sampTab_TM_053018.rda")
expTMraw = utils_loadObject("expMatrix_TM_Raw_Oct_12_2018.rda") # reload training

# Load the ortholog table and convert human gene names to mouse ortholog names, 
# and limit analysis to genes in common between the training and query data.
oTab = utils_loadObject("human_mouse_genes_Jul_24_2018.rda")
dim(oTab)

aa = csRenameOrth(expQuery, expTMraw, oTab)
expQueryOrth = aa[['expQuery']]
expTrainOrth = aa[['expTrain']]

# Limit anlaysis to a subset of the TM cell types
cts = c("B cell",  "cardiac muscle cell", "endothelial cell", "erythroblast", "granulocyte", "hematopoietic precursor cell", "late pro-B cell", "limb_mesenchymal", "macrophage", "mammary_basal_cell", "monocyte", "natural killer cell", "T cell", "trachea_epithelial", "trachea_mesenchymal")

stTM2 = filter(stTM, newAnn %in% cts)
stTM2 = droplevels(stTM2)
rownames(stTM2) = as.vector(stTM2$cell) # filter strips rownames

expTMraw2 = expTrainOrth[,rownames(stTM2)]
dim(expTMraw2)

# Train Classifier
stList = splitCommon(stTM2, ncells=100, dLevel="newAnn")
stTrain = stList[[1]]
expTrain = expTMraw2[,rownames(stTrain)]

system.time(class_info2<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "newAnn", colName_samp = "cell"))

# Apply to held out data
#validate data
stTestList = splitCommon(stList[[2]], ncells=100, dLevel="newAnn") 
stTest = stTestList[[1]]
expTest = expTMraw2[,rownames(stTest)]

#predict
system.time(classRes_val_all2 <- scn_predict(class_info2[['cnProc']], expTest, nrand = 50))

# Assess classifier
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all2, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "newAnn", classQuery = "newAnn", nRand = 50)

plot_PRs(tm_heldoutassessment)
plot_metrics(tm_heldoutassessment)

# Classification result heatmap
nrand=50
sla = as.vector(stTest$newAnn)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand)
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand)
unique(sla)
unique(stTest$newAnn)
# heatmap classification result
sc_hmClass(classRes_val_all2, sla, max=300, font=7, isBig=TRUE)
# attribute plot
plot_attr(classRes_val_all2, stTest, nrand=nrand, dLevel="newAnn", sid="cell")

# apply to human query data
stQuery$description = as.character(stQuery$description)
stQuery[which(stQuery$description == "NK cell"), "description"] = "natural killer cell"

nqRand = 50
system.time(crHS <- scn_predict(class_info2[['cnProc']], expQueryOrth, nrand=nqRand))

# Assess classifier with external dataset
tm_pbmc_assessment = assess_comm(ct_scores = crHS, stTrain = stTrain, stQuery = stQuery, classTrain = "newAnn",classQuery="description",dLevelSID="sample_name")
plot_PRs(tm_pbmc_assessment)
plot_metrics(tm_pbmc_assessment)

# Classification result heatmap
sgrp = as.vector(stQuery$prefix)
names(sgrp) = as.vector(stQuery$sample_name)
grpRand = rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

sc_hmClass(crHS, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
# Classification violin plot
sc_violinClass(sampTab = stQuery, classRes = crHS, sid = "sample_name", dLevel = "description")
# Classification violin plot with selected cluster
sc_violinClass(stQuery, crHS, sid = "sample_name", dLevel = "description", ncol = 12, sub_cluster = "B cell")
# Attribution plot
plot_attr(crHS, stQuery, nrand=nqRand, sid="sample_name", dLevel="description")
# sub cluster
plot_attr(sampTab = stQuery, classRes = crHS, sid = "sample_name", dLevel = "description", nrand = 50, sub_cluster = c("B cell", "T cell"))
# umap by category
system.time(umPrep_HS<-prep_umap_class(crHS, stQuery, nrand=nqRand, dLevel="description", sid="sample_name", topPC=5))
plot_umap(umPrep_HS)
# Heatmap top pairs genes for training sample average
system.time(gpTab2 <- compareGenePairs(query_exp = expQueryOrth, training_exp = 
                                         expTrainOrth, training_st = stTrain, 
                                       classCol = "newAnn", sampleCol = "cell", 
                                       RF_classifier = class_info2$cnProc$classifier, 
                                       numPairs = 20, trainingOnly = FALSE))

sgrp = as.vector(stQuery$prefix)
names(sgrp) = rownames(stQuery)
train2 = findAvgLabel(gpTab2, stTrain = stTrain, dLevel = "newAnn")
sgrp = append(sgrp, train2)

hm_gpa_sel(gpTab2, genes = class_info2$cnProc$xpairs, grps = sgrp, maxPerGrp = 5)

# How to calibrate/make sense of a given SCN score

# this function aims to give you a sense of how precise/sensitive SCN is with the 
# assigned score of a given cell type for a cell
#tm_assess_matrix = tm_heldoutassessment$nonNA_PR

#tm_assess_matrix is a held_out assessment metric extracted from tm_heldoutassessment, 
# which is already stored in SCN.
#e_assess_matrix is also provided for a gastrulation SCN classifier 
score = 0.6
celltype = "B cell"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#[1] "SCN score of 0.6 for cell type B cell has precision of 0.979 ~ 0.979 and sensitivity of 0.93 ~ 0.93"
calibration
