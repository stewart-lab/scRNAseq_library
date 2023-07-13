library(Seurat)
# get gamm metadata
setwd("/w5home/bmoore/scRNAseq/GAMM/human_data/reh_cellrep_2020/")
gamm <- readRDS(file = "GAMM_S2_clabeled-clusters_0.14.rds")
gamm@meta.data
gamm@active.ident
Idents(gamm)
gamm$CellType <- Idents(gamm)
saveRDS(gamm, file = "GAMM_S2_clabeled-clusters_0.14.rds")
gamm.metadata<- as.data.frame(gamm@meta.data)
write.table(gamm.metadata, file = "gamm_manual_annot_metadata_c0.14.txt", sep="\t", quote=F, row.names = T)
colnames(gamm@assays$RNA@data)
