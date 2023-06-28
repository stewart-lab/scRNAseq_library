# singleCellNet

install.packages("devtools")
devtools::install_github("pcahan1/singleCellNet")
library(singleCellNet)

# object conversion
#exp_type options can be: counts, data, and scale.data if they are available in your sce object
scefile <- extractSCE(sce_object, exp_slot_name = "counts") 
sampTab = scefile$sampTab
expDat = scefile$expDat

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
