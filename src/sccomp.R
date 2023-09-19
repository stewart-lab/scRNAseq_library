if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("sccomp")
library(sccomp)
library(Seurat)
library(tidyverse)
setwd("/Users/bmoore/Desktop/scRNAseq/GAMM/human_ref")
output <- "cell_composition/"

# read in data
# from seurat mapping
gamms2<- readRDS(file = "output_seurat_mapping_20230913_100651_cc/gamms2_cca_pred.rds")
# from scPred
gamms2 <- readRDS(file = "output_scPred_20230913_131047_filt/GAMM_S2_scpred_c0.5.rds")
# human reference data
human_D205.seurat<- readRDS(file = "output_20230711_155556/human_D205_umap.rds")
human_D205.seurat<- subset(human_D205.seurat, subset = type != 'AC2')
human_D205.seurat<- subset(human_D205.seurat, subset = type != 'T2')
human_D205.seurat<- subset(human_D205.seurat, subset = type != 'Midbrain')
human_D205.seurat<- subset(human_D205.seurat, subset = type != 'miG')
# example data
# example.data <- load("seurat_obj.rda")
# table(seurat_obj$sample)
# table(seurat_obj$type)
# table(seurat_obj$cell_group)
# gamm s1 and s2
gamms1 <- readRDS(file = "GAMM_S1/output_20230829_135606/seurat_obj_labeled.rds")
gamms2 <- readRDS(file = "GAMM_S2/output_20230830_155530/seurat_obj_labeled.rds")
# put labeled clusters in metadata
gamms1@active.ident
Idents(gamms1)
gamms1 <- AddMetaData(gamms1, metadata = gamms1@active.ident, col.name= "CellType")
table(gamms1$CellType)
gamms2@active.ident
gamms2 <- AddMetaData(gamms2, metadata = gamms2@active.ident, col.name= "CellType")
table(gamms2$CellType)
# check location of annotation
# seurat mapping
gamms2$predicted.id
#scPred
gamms2$scpred_prediction
#human
human_D205.seurat$type

# merge data sets
# pig and human
seurat.combined <- merge(gamms2, y = human_D205.seurat, add.cell.ids = c("pig", "human"), project = "combined")
head(colnames(seurat.combined))
table(seurat.combined$orig.ident)
unique(seurat.combined$orig.ident)
unique(sapply(X = strsplit(colnames(seurat.combined), split = "_"), FUN = "[", 1))
table(Idents(object = seurat.combined))
# gamm s1 and s2
seurat.combined <- merge(gamms1, y = gamms2, add.cell.ids = c("s1", "s2"), project = "combined")
head(colnames(seurat.combined))
table(seurat.combined$orig.ident)
unique(seurat.combined$orig.ident)
unique(sapply(X = strsplit(colnames(seurat.combined), split = "_"), FUN = "[", 1))
table(Idents(object = seurat.combined))
table(seurat.combined$CellType)
# set identity if needed
#Idents(object = seurat.combined) <- "orig.ident"
# table(Idents(object = seurat.combined))

# Rename identity classes
seurat.combined <- RenameIdents(object = seurat.combined, `human_D205` = "human",
                                `gamm_S2` = "pig")
table(Idents(object = seurat.combined))
# stash identities
seurat.combined[["idents"]] <- Idents(object = seurat.combined)
unique(seurat.combined$idents)
# gamm s1 and s2
Idents(object = seurat.combined) <- "orig.ident"
table(Idents(object = seurat.combined))
seurat.combined <- RenameIdents(object = seurat.combined, `output_S1-1_mm_mt_GeneFull` = "s1",
                                `output_S1-2_mm_mt_GeneFull` = "s1", `output_S2-1_mm_mt_GeneFull` = "s2",
                                `output_S2-2_mm_mt_GeneFull` = "s2")
table(Idents(object = seurat.combined))
# stash identities
seurat.combined[["idents"]] <- Idents(object = seurat.combined)
unique(seurat.combined$idents)

# check cell type annotations
# seurat mapping
unique(seurat.combined$predicted.id)
# scpred
unique(seurat.combined$scpred_prediction)
# human
unique(seurat.combined$type)

# subset to remove NAs
seurat.combined <- subset(x= seurat.combined, idents != "NA")
# seurat mapping
seurat.pig <- subset(x= seurat.combined, predicted.id != "NA")
unique(seurat.pig$predicted.id)
# scpred
seurat.pig <- subset(x= seurat.combined, scpred_prediction != "NA")
unique(seurat.pig$scpred_prediction)
# human
seurat.human <- subset(x= seurat.combined, type != "NA")
unique(seurat.human$type)

# plot cell type proportions for human and pig
library(devtools)
devtools::install_github("Oshlack/speckle")
library(speckle)
library(ggplot2)
library(limma)
# seurat mapping
pig.plot <- plotCellTypeProps(x = seurat.pig, clusters = seurat.pig$predicted.id, sample = seurat.pig$ident)
# scpred
pig.plot <- plotCellTypeProps(x = seurat.pig, clusters = seurat.pig$scpred_prediction, sample = seurat.pig$ident)
# human
hum.plot <- plotCellTypeProps(x = seurat.human, clusters = seurat.human$type, sample = seurat.human$orig.ident)
print(pig.plot + hum.plot)
# save plot
pdf(paste0(output, "celltype_proportions.pdf"), width = 8, height = 6)
print(pig.plot + hum.plot)
dev.off()
# gamm s1 and s2
prop.plot <- plotCellTypeProps(x = seurat.combined, clusters = seurat.combined$CellType, sample = seurat.combined$idents)
prop.plot
pdf(paste0(output, "celltype_proportions.pdf"), width = 8, height = 6)
print(prop.plot)
dev.off()
# put predicted id and type in new category together
colnames(seurat.combined@meta.data)
# seurat mapping
New_idents <- c(seurat.pig$predicted.id,seurat.human$type)
# scpred
New_idents <- c(seurat.pig$scpred_prediction,seurat.human$type)
table(New_idents)
seurat.combined@meta.data$"Newidents" <- as.factor(New_idents)
seurat.combined <- subset(x= seurat.combined, orig.ident != "NA")
# save seurat object
saveRDS(seurat.combined, file = paste0(output,"pig-human.combined.rds"))
saveRDS(seurat.combined, file = paste0(output,"gamm_s1-s2.combined.rds"))

# now we are ready to do sccomp
# model composition
comp.tibble <- seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ idents, 
    .sample =  orig.ident, 
    .cell_group = Newidents, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
# gamm s1 and s2
comp.tibble <- seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ idents, 
    .sample =  orig.ident, 
    .cell_group = CellType, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
# write results
comp.table <- as.data.frame(comp.tibble)
comp.table.1 <- comp.table[,1:8]
write.table(comp.table.1, file= paste0(output, "composition_result_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")
# get counts from significant result:
pr.table <- as.data.frame(comp.table["8","count_data"])
write.table(pr.table, file= paste0(output, "PR_composition_result_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")
# model contrast
contrast.tibble <- seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ 0 + idents, 
    contrasts =  c("pig - human", "human - pig"),
    .sample = orig.ident,
    .cell_group = Newidents, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
contrast.table <- as.data.frame(contrast.tibble)
contrast.table.1 <- contrast.table[,1:10]
write.table(contrast.table.1, file= paste0(output, "contrast_result_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")
# categorical factor (a Bayesian Anova)
# is the model with variables significantly different than the model without?
library(loo)

# Fit first model
model_with_factor_association = 
  seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ idents, 
    .sample =  orig.ident, 
    .cell_group = Newidents, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )
# gamm s1 and s2
model_with_factor_association = 
  seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ idents, 
    .sample =  orig.ident, 
    .cell_group = CellType, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )
# Fit second model
model_without_association = 
  seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ 1, 
    .sample =  orig.ident, 
    .cell_group = Newidents, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 , 
    enable_loo = TRUE
  )
# gamm s1 and s2
model_without_association = 
  seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ 1, 
    .sample =  orig.ident, 
    .cell_group = CellType, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 , 
    #enable_loo = TRUE
  )
# Compare models
comparison <- loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
summary(comparison)
comp_df <- as.data.frame(comparison)
write.table(comp_df, file= paste0(output, "model_comparison_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")
# is elpd_diff/se_diff > abs(5)?
c.res <- comp_df$elpd_diff[2]/comp_df$se_diff[2]
print(paste0("if ",c.res," is greater than abs(5), significant"))
#-7.4/2.5
# [1] -2.96, so not significant

# Differential variablity- model cell-group variability dependent on idents (pig or human)
res = 
  seurat.combined |>
  sccomp_glm( 
    formula_composition = ~ idents, 
    formula_variability = ~ idents,
    .sample = orig.ident,
    .cell_group = Newidents,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
res.table <- as.data.frame(res)
res.table.1 <- res.table[,1:17]
write.table(res.table.1, file= paste0(output, "celltypexspecies_composition_result_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")

# plot results
plots = plot_summary(res)
# A plot of group proportion, faceted by groups. The blue boxplots represent the 
# posterior predictive check. If the model is likely to be descriptively adequate 
# to the data, the blue box plot should roughly overlay with the black box plot, 
# which represents the observed data. The outliers are coloured in red. A box plot 
# will be returned for every (discrete) covariate present in formula_composition. 
#The colour coding represents the significant associations for composition and/or 
# variability.
pdf(paste0(output, "signif_associations_boxplot.pdf"), width = 8, height = 6)
print(plots$boxplot)
dev.off()
# A plot of estimates of differential composition (c_) on the x-axis and 
# differential variability (v_) on the y-axis. The error bars represent 95% 
# credible intervals. The dashed lines represent the minimal effect that the 
#vhypothesis test is based on. An effect is labelled as significant if bigger 
# than the minimal effect according to the 95% credible interval. Facets represent 
# the covariates in the model.
pdf(paste0(output, "credible_intervals_1d.pdf"), width = 8, height = 6)
print(plots$credible_intervals_1D)
dev.off()
# Plot 2D significance plot. Data points are cell groups. Error bars are the 95% 
# credible interval. The dashed lines represent the default threshold fold change 
# for which the probabilities (c_pH0, v_pH0) are calculated. pH0 of 0 represent 
# the rejection of the null hypothesis that no effect is observed. The differential 
# variability estimates are reliable only if the linear association between mean 
# and variability for (intercept) (left-hand side facet) is satisfied.
pdf(paste0(output, "credible_intervals_2d.pdf"), width = 8, height = 6)
print(plots$credible_intervals_2D)
dev.off()
# markov chain plot
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")

