

read_aligned_data <- function(data_directory){
    data <- Read10X(data_directory)
}

prep_seurat_and_soupX <- function(data.raw, data, project) {
  # Create SoupChannel object
  sc <- SoupChannel(data.raw, data)

  
  # Create a Seurat object without filtering
  gamm <- CreateSeuratObject(counts = data, project = project)
  
  # Perform transformations and find clusters
  gamm <- SCTransform(gamm, verbose = F)
  gamm <- RunPCA(gamm, verbose = F)
  gamm <- RunUMAP(gamm, dims = 1:30, verbose = F)
  gamm <- FindNeighbors(gamm, dims = 1:30, verbose = F)
  gamm <- FindClusters(gamm, verbose = T)
  
  # Extract clusters from Seurat object and add to SoupChannel
  meta <- gamm@meta.data
  umap <- gamm@reductions$umap@cell.embeddings
  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  
  # Estimate contamination fraction
  sc <- autoEstCont(sc)
  out = adjustCounts(sc, roundToInt = TRUE)
  list(gamm = gamm, meta = meta, umap = umap, out = out)
}


create_seurat_and_sce <- function(out1, out2) {

  # Create singular Seurat objects for each data set
  gamm1_seu <- CreateSeuratObject(counts = out1, project = "gamm_s1-1")
  gamm2_seu <- CreateSeuratObject(counts = out2, project = "gamm_s1-2")

  # Convert Seurat objects to SingleCellExperiment objects
  sce1 <- as.SingleCellExperiment(gamm1_seu)
  sce2 <- as.SingleCellExperiment(gamm2_seu)

  list(gamm1_seu = gamm1_seu, gamm2_seu = gamm2_seu, sce1 = sce1, sce2 = sce2)
}

run_scDblFinder_and_merge <- function(sce1, sce2) {
  set.seed(1234)

  # Run scDblFinder
  sce1 <- scDblFinder(sce1)
  sce2 <- scDblFinder(sce2)

  table1 <- table(call=sce1$scDblFinder.class)
  table2 <- table(call=sce2$scDblFinder.class)

  # Filter out doublets
  sce1 <- sce1[, sce1$scDblFinder.class == "singlet"]
  sce2 <- sce2[, sce2$scDblFinder.class == "singlet"]

  # Convert back to Seurat
  gamm1 <- as.Seurat(sce1)
  gamm2 <- as.Seurat(sce2)

  # Merge Seurat objects
  gamm <- merge(gamm1, y = gamm2)
  
  list(gamm = gamm, gamm1 = gamm1, gamm2 = gamm2, table1=table1, table2=table2)
}

filter_cells <- function(gamm, mt.list, lower.nFeature = 200, upper.nFeature = 8000, max.percent.mt = 5) {
  # Check if all features are in the matrix
  if (!all(mt.list %in% rownames(gamm))) {
    for (mt in mt.list) {
      print(paste(mt, " in rownames: ", all(mt %in% rownames(gamm))))
    }
  }

  # Calculate percent.mt and add to the metadata
  percent_mt <- PercentageFeatureSet(gamm, features = mt.list, assay = 'RNA')
  gamm[["percent.mt"]] <- percent_mt
  pre_filter_plot <- create_feature_scatter_plot(gamm, 'nCount_RNA', 'percent.mt')
  # Filter cells based on nFeature_RNA and percent.mt
  gamm <- subset(gamm, subset = nFeature_RNA > lower.nFeature & nFeature_RNA < upper.nFeature & percent.mt < max.percent.mt)

  return(list(gamm=gamm, pre_filter_plot=pre_filter_plot))
}

normalize_data <- function(gamm) {
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(gamm)

  # Cluster the cells for scran
  clusters <- quickCluster(sce, use.ranks = FALSE, min.size = 100)

  # Calculate size factors per cell
  sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1)

  # Apply size factors to generate log normalized data
  sce <- logNormCounts(sce)

  # Replace Seurat normalized values with scran
  gamm[["RNA"]] <- SetAssayData(gamm[["RNA"]],
                                slot = "data", 
                                new.data = logcounts(sce))

  gamm$sizeFactors <- sizeFactors(sce)
  gamm <- UpdateSeuratObject(gamm)

  gamm
}

feature_selection <- function(gamm, analysis_type) {
  if (analysis_type == "Seurat") {
    # Seurat method
    gamm <- FindVariableFeatures(gamm, selection.method = "vst", nfeatures = 2000)
  } else if (analysis_type == "Scry") {
    # Scry method
    m <- GetAssayData(gamm, slot = "counts", assay = "RNA")
    devi <- scry::devianceFeatureSelection(m)
    dev_ranked_genes <- rownames(gamm)[order(devi, decreasing = TRUE)]
    topdev <- head(dev_ranked_genes, 2000)
    # replace variable features with the deviance ranked genes
    VariableFeatures(gamm) <- topdev
    gamm <- UpdateSeuratObject(gamm)
  } else {
    stop("Invalid analysis_type. Please choose 'Seurat' or 'Scry'.")
  }
  
  gamm
}


scale_data <- function(gamm) {
    # Get all gene names
    all.genes <- rownames(gamm)
  
    # Scale the data
    gamm <- ScaleData(gamm, features = all.genes)
  
    return(gamm)
}

run_and_visualize_pca <- function(seurat_obj, top_n_dims=5, top_n_features=5, heatmap_dims=1:15, num_cells=500) {

  # Perform PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  # Print top dimensions and features
 print(seurat_obj[["pca"]], dims = 1:top_n_dims, nfeatures = top_n_features)
  
  # Visualize feature loadings
  feature_loadings <- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
  
  # Generate scatter plot of PCA results
  pca_scat_plot <- DimPlot(seurat_obj, reduction = "pca")
  
  # save PCA heat map
  
  pca_heat_map <- DimHeatmap(seurat_obj, nfeatures = 5,dims = heatmap_dims, cells = num_cells, balanced = TRUE, fast = FALSE, combine = TRUE)
  
  # return all objects created in this function 
  return(list( seurat_obj=seurat_obj,feature_loadings=feature_loadings, pca_scat_plot=pca_scat_plot, pca_heat_map=pca_heat_map))
  
}


perform_batch_correction <- function(seurat_obj,dims.use=1:20) {

 
  
  # Check PCA plot
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1_pre <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = "orig.ident")
  p2_pre <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = "orig.ident", pt.size = .1)
  
  # Run Harmony
  options(repr.plot.height = 2.5, repr.plot.width = 6)
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", dims.use=dims.use, max.iter.harmony = 50)
  
  # Access Harmony embeddings
  harmony_embeddings <- Embeddings(seurat_obj, 'harmony')
  
  # Check out new PCA plot
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1_post <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2_post <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
 
  
  return(list(seurat_obj=seurat_obj, harmony_embeddings=harmony_embeddings, p1_pre=p1_pre, p2_pre=p2_pre, p1_post=p1_post, p2_post=p2_post))
}


perform_clustering <- function(seurat_obj, num.replicate = 100, dims = 1:20, dims_umap = 1:15, resolution = 0.5) {

  # Perform JackStraw
  seurat_obj <- JackStraw(seurat_obj, num.replicate = num.replicate)
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = dims)
  
  # Visualize JackStraw results
  jack_straw <- JackStrawPlot(seurat_obj, dims = dims)
  
  # Generate Elbow plots
  elbow_harm <-ElbowPlot(seurat_obj, reduction='harmony')
  elbow_pca <- ElbowPlot(seurat_obj, reduction='pca')
  
  # Perform K-nearest neighbor (KNN) graph using harmony embeddings
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_umap, reduction='harmony')
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_umap, reduction = "harmony")
  
  # Plot UMAP results
  options(repr.plot.height = 4, repr.plot.width = 10)
  umap_lanes <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = .1)
  
  # Cluster cells
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = "leiden")

  # Visualize clusters
  options(repr.plot.height = 4, repr.plot.width = 10)
  umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .1)
  
  return(list(seurat_obj=seurat_obj, jack_straw=jack_straw, elbow_harm=elbow_harm, elbow_pca=elbow_pca, umap_lanes=umap_lanes, umap_clusters=umap_clusters))
}


find_differential_expression <- function(seurat_obj, logfc.threshold = 0.25, min.pct = 0.25, top_n = 10) {
  # Find markers for every cluster compared to all remaining cells, only positive ones
  seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)

  # Print the top 'n' markers for each cluster by avg_log2FC
  seurat_obj.markers %>%
    group_by(cluster) %>%
    slice_max(n = top_n, order_by = avg_log2FC) %>%
    print(n = top_n * length(unique(seurat_obj.markers$cluster)))

  # Find ROC markers for cluster 1
  cluster1.markers <- FindMarkers(seurat_obj, ident.1 = 1, logfc.threshold = logfc.threshold, test.use = "roc", only.pos = TRUE)

  # Visualize marker expression
  VlnPlot(seurat_obj, features = c("SFRP2", "TRPM3"))

  # Generate an expression heatmap for the top 'n' markers for each cluster
  seurat_obj.markers %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC) -> top_n_markers

  DoHeatmap(seurat_obj, features = top_n_markers$gene) + NoLegend()

  return(list("AllMarkers" = seurat_obj.markers, "Cluster1Markers" = cluster1.markers, "TopNMarkers" = top_n_markers))
}

create_feature_scatter_plot <- function(obj, feature1, feature2) {

 plot <- FeatureScatter(object = obj, feature1 = feature1, feature2 = feature2)
  return(plot)
}