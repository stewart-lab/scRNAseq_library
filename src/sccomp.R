if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("sccomp")
library(sccomp)
library(Seurat)
library(tidyverse)
setwd("/Users/bmoore/Desktop/scRNAseq/GAMM/human_ref")
# read in data
gamms2<- readRDS(file = "output_seurat_mapping_20230824_094352/gamms2_cca_pred.rds")
human_D205.seurat<- readRDS(file = "human_D205_umap.rds")
example.data <- load("seurat_obj.rda")
table(seurat_obj$sample)
table(seurat_obj$type)
table(seurat_obj$cell_group)
# check location of annotation
gamms2$predicted.id
human_D205.seurat$type

# merge data sets
seurat.combined <- merge(gamms2, y = human_D205.seurat, add.cell.ids = c("pig", "human"), project = "combined")
head(colnames(seurat.combined))
table(seurat.combined$orig.ident)
unique(seurat.combined$orig.ident)
unique(sapply(X = strsplit(colnames(seurat.combined), split = "_"), FUN = "[", 1))
table(Idents(object = seurat.combined))
# set identity
Idents(object = seurat.combined) <- "orig.ident"
table(Idents(object = seurat.combined))
# Rename identity classes
seurat.combined <- RenameIdents(object = seurat.combined, `D205` = "human",
                                `gamm_s1-1` = "pig")
table(Idents(object = seurat.combined))
# stash identities
seurat.combined[["idents"]] <- Idents(object = seurat.combined)
unique(seurat.combined$idents)
# check cell type annotations
unique(seurat.combined$predicted.id)
unique(seurat.combined$type)
# subset to remove NAs
seurat.combined <- subset(x= seurat.combined, idents != "NA")
seurat.pig <- subset(x= seurat.combined, predicted.id != "NA")
unique(seurat.pig$predicted.id)
seurat.human <- subset(x= seurat.combined, type != "NA")
unique(seurat.human$type)
# plot cell type proportions for human and pig
library(devtools)
devtools::install_github("Oshlack/speckle")
library(speckle)
library(ggplot2)
library(limma)
pig.plot <- plotCellTypeProps(x = seurat.pig, clusters = seurat.pig$predicted.id, sample = seurat.pig$ident)
hum.plot <- plotCellTypeProps(x = seurat.human, clusters = seurat.human$type, sample = seurat.human$orig.ident)
print(pig.plot + hum.plot)
# put predicted id and type in new category together
seurat.combined@meta.data
New_idents <- c(seurat.pig$predicted.id,seurat.human$type)
table(New_idents)
seurat.combined@meta.data$"Newidents" <- as.factor(New_idents)
# save seurat object
saveRDS(seurat.combined, file = "cell_composition/pig-human.combined.rds")
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
# write results
comp.table <- as.data.frame(comp.tibble)
comp.table.1 <- comp.table[,1:17]
write.table(comp.table.1, file= "composition_result_table.txt", quote = F, col.names = TRUE, row.names= F, sep= "\t")
# get counts from significant result:
pr.table <- as.data.frame(comp.table[12,18])
write.table(pr.table, file= "PR_composition_result_table.txt", quote = F, col.names = TRUE, row.names= F, sep= "\t")
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

# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
# is elpd_diff/se_diff > abs(5)?
-7.4/2.5
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
write.table(res.table.1, file= "celltypexspecies_composition_result_table.txt", quote = F, col.names = TRUE, row.names= F, sep= "\t")

# plot results
plots = plot_summary(res)
# A plot of group proportion, faceted by groups. The blue boxplots represent the 
# posterior predictive check. If the model is likely to be descriptively adequate 
# to the data, the blue box plot should roughly overlay with the black box plot, 
# which represents the observed data. The outliers are coloured in red. A box plot 
# will be returned for every (discrete) covariate present in formula_composition. 
#The colour coding represents the significant associations for composition and/or 
# variability.
plots$boxplot
# A plot of estimates of differential composition (c_) on the x-axis and 
# differential variability (v_) on the y-axis. The error bars represent 95% 
# credible intervals. The dashed lines represent the minimal effect that the 
#vhypothesis test is based on. An effect is labelled as significant if bigger 
# than the minimal effect according to the 95% credible interval. Facets represent 
# the covariates in the model.
plots$credible_intervals_1D
# Plot 2D significance plot. Data points are cell groups. Error bars are the 95% 
# credible interval. The dashed lines represent the default threshold fold change 
# for which the probabilities (c_pH0, v_pH0) are calculated. pH0 of 0 represent 
# the rejection of the null hypothesis that no effect is observed. The differential 
# variability estimates are reliable only if the linear association between mean 
# and variability for (intercept) (left-hand side facet) is satisfied.
plots$credible_intervals_2D
# marcov chain plot
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
