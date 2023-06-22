
read_aligned_data <- function(base_directory, project_name) {
  filtered_data_directory <- paste0(base_directory, "/filtered/")
  raw_data_directory <- paste0(base_directory, "/raw/")

  list(
    filtered = Read10X(filtered_data_directory),
    raw = Read10X(raw_data_directory),
    project = project_name
  )
}

prep_seurat_and_soupX <- function(data.raw, data, project, dims = 1:30) {
  # Create SoupChannel object
  sc <- SoupChannel(data.raw, data)

  # Create a Seurat object without filtering
  gamm <- CreateSeuratObject(counts = data, project = project)
  
  # Perform transformations and find clusters
  gamm <- SCTransform(gamm, verbose = F)
  gamm <- RunPCA(gamm, verbose = F)
  
  # Use the provided 'dims' parameter in the following functions
  gamm <- RunUMAP(gamm, dims = dims, verbose = F)
  gamm <- FindNeighbors(gamm, dims = dims, verbose = F)
  
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

create_seurat_and_sce <- function(out, project, feature_set) {
  # Create singular Seurat object
  seu <- CreateSeuratObject(counts = out, project = project)

  # Convert Seurat objects to SingleCellExperiment object
  sce <- as.SingleCellExperiment(seu)
  
  # Post SoupX QC
  combine_feature_plots(filtered_list_of_objs, feature_set, file_name_base = "post_soupx_qc")
  
  list(seu = seu, sce = sce)
}

run_scDblFinder_and_filter <- function(sce) {
  set.seed(1234)
  sce <- scDblFinder(sce)
  table <- table(call=sce$scDblFinder.class)
  sce <- sce[, sce$scDblFinder.class == "singlet"]
  list(sce = sce, table = as.data.frame(table))
}

run_scDblFinder_and_merge <- function(sce_list, save_plot = TRUE, file_name = "after_dbl_removal_and_merge", path = "plot_outputs/") {

  result_list <- list()
  for (i in seq_along(sce_list)) {
    sce_list[[i]] <- scDblFinder(sce_list[[i]])
    result_list[[paste0("table", i)]] <- as.data.frame(table(call=sce_list[[i]]$scDblFinder.class))
    sce_list[[i]] <- sce_list[[i]][, sce_list[[i]]$scDblFinder.class == "singlet"]
    result_list[[paste0("gamm", i)]] <- as.Seurat(sce_list[[i]])
  }
  result_list$gamm <- merge(result_list$gamm1, y = result_list$gamm2)
  if (save_plot) {
    create_feature_scatter_plot(obj = result_list$gamm, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', 
                                file_name = file_name, path = path)
  }
  
  return(result_list)
}

filter_cells <- function(gamm, mt.list, lower.nFeature = 200, upper.nFeature = 8000, max.percent.mt = 5,
                         save_plots = TRUE, path = "plot_outputs/") {
  # Check if all features are in the matrix
  if (!all(mt.list %in% rownames(gamm))) {
    for (mt in mt.list) {
      print(paste(mt, " in rownames: ", all(mt %in% rownames(gamm))))
    }
  }

  # Calculate percent.mt and add to the metadata
  percent_mt <- PercentageFeatureSet(gamm, features = mt.list, assay = 'RNA')
  gamm[["percent.mt"]] <- percent_mt
  
  # Create pre-filter plot
  pre_filter_plot <- create_feature_scatter_plot(gamm, 'nCount_RNA', 'percent.mt', file_name = "mt_unfiltered", save = save_plots, path = path)
  
  # Filter cells based on nFeature_RNA and percent.mt
  gamm <- subset(gamm, subset = nFeature_RNA > lower.nFeature & nFeature_RNA < upper.nFeature & percent.mt < max.percent.mt)
  
  # Create post-filter plots
  post_filter_plot1 <- create_feature_scatter_plot(gamm, 'nCount_RNA', 'percent.mt', file_name = "mt_filtered1", save = save_plots, path = path)
  post_filter_plot2 <- create_feature_scatter_plot(gamm, 'nCount_RNA', 'nFeature_RNA', file_name = "mt_filtered2", save = save_plots, path = path)
  
  # Combine plots and delete individual ones
  seurat_objs_list <- list(seu1 = list(seu = gamm), seu2 = list(seu = gamm))
  feature_set1 <- list(feature1 = 'nCount_RNA', feature2 = 'percent.mt')
  feature_set2 <- list(feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
  combine_feature_plots(seurat_objs_list, feature_set1, feature_set2, same_feature_set = FALSE, file_name_base = "mt_filtered", path = path)
  
  return(gamm)
}

normalize_data <- function(gamm, feature = "ECHS1", path = "plot_outputs/", min_size = 100, min_mean = 0.1) {
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(gamm)

  # Cluster the cells for scran
  clusters <- quickCluster(sce, use.ranks = FALSE, min.size = min_size)

  # Calculate size factors per cell
  sce <- computeSumFactors(sce, clusters = clusters, min.mean = min_mean)

  # Apply size factors to generate log normalized data
  sce <- logNormCounts(sce)

  # Replace Seurat normalized values with scran
  gamm[["RNA"]] <- SetAssayData(gamm[["RNA"]],
                                slot = "data", 
                                new.data = logcounts(sce))

  gamm$sizeFactors <- sizeFactors(sce)
  gamm <- UpdateSeuratObject(gamm)
  
  # Create violin plots pre and post normalization
  vin_pre <- VlnPlot(gamm, feature, slot = "counts")
  vin_post <- VlnPlot(gamm, feature, slot = "data")
  
  # Save violin plots
  pdf(file = paste0(path, "violin_pre_norm.pdf"), width = 8, height = 6)
  print(vin_pre)
  dev.off()
  
  pdf(file = paste0(path, "violin_post_norm.pdf"), width = 8, height = 6)
  print(vin_post)
  dev.off()

  return(gamm)
}

feature_selection <- function(gamm, analysis_type, n_features = 2000) {
  if (analysis_type == "Seurat") {
    # Seurat method
    gamm <- FindVariableFeatures(gamm, selection.method = "vst", nfeatures = n_features)
  } else if (analysis_type == "Scry") {
    # Scry method
    m <- GetAssayData(gamm, slot = "counts", assay = "RNA")
    devi <- scry::devianceFeatureSelection(m)
    dev_ranked_genes <- rownames(gamm)[order(devi, decreasing = TRUE)]
    topdev <- head(dev_ranked_genes, n_features)
    # replace variable features with the deviance ranked genes
    VariableFeatures(gamm) <- topdev
    gamm <- UpdateSeuratObject(gamm)
  } else {
    stop("Invalid analysis_type. Please choose 'Seurat' or 'Scry'.")
  }
  
  return(gamm)
}

scale_data <- function(gamm) {
    # Get all gene names
    all.genes <- rownames(gamm)
  
    # Scale the data
    gamm <- ScaleData(gamm, features = all.genes)
  
    return(gamm)
}

run_and_visualize_pca <- function(seurat_obj, path = "plot_outputs/", top_n_dims=5, top_n_features=5, heatmap_dims=1:15, num_cells=500) {

  # Perform PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  # Print top dimensions and features
  print(seurat_obj[["pca"]], dims = 1:top_n_dims, nfeatures = top_n_features)
  
  # Visualize feature loadings
  pdf(paste0(path, "feature_loadings.pdf"), width = 8, height = 6)
  print(VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca"))
  dev.off()

  # Generate scatter plot of PCA results
  pdf(paste0(path, "pca_scatter_plot.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "pca"))
  dev.off()
  
  # Save PCA heat map
  pdf(paste0(path, "pca_heat_map.pdf"), width = 8, height = 6)
  print(DimHeatmap(seurat_obj, nfeatures = 5,dims = heatmap_dims, cells = num_cells, balanced = TRUE, fast = FALSE, combine = TRUE))
  dev.off()
  
  # Return the updated Seurat object
  return(seurat_obj)
}

perform_batch_correction <- function(seurat_obj, dims.use=1:20, max_iter=50, path = "plot_outputs/") {
  
  # Generate pre-batch correction PCA plot
  pdf(paste0(path, "batch_uncorrected_pca.pdf"), width = 8, height = 6)
  p1_pre <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = "orig.ident")
  p2_pre <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = "orig.ident", pt.size = .1)
  print(p1_pre + p2_pre)
  dev.off()
  
  # Run Harmony
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", dims.use=dims.use, max.iter.harmony = max_iter)
  
  # Access Harmony embeddings
  harmony_embeddings <- Embeddings(seurat_obj, 'harmony')
  
  # Generate post-batch correction PCA plot
  pdf(paste0(path, "batch_corrected_pca.pdf"), width = 8, height = 6)
  p1_post <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2_post <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  print(p1_post + p2_post)
  dev.off()
  
  # Return the updated Seurat object and harmony embeddings
  return(list(seurat_obj=seurat_obj, harmony_embeddings=harmony_embeddings))
}

perform_clustering <- function(seurat_obj, path = "plot_outputs/") {
  
  num_replicate = config$perform_clustering$num.replicate
  dims = 1:config$perform_clustering$dims
  dims_umap = 1:config$perform_clustering$dims_umap
  resolution = config$perform_clustering$resolution
  
  # Perform JackStraw
  seurat_obj <- JackStraw(seurat_obj, num.replicate = num_replicate)
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = dims)
  
  # Save JackStraw plot
  pdf(paste0(path, "jack_straw.pdf"), width = 8, height = 6)
  jack_straw <- JackStrawPlot(seurat_obj, dims = dims)
  print(jack_straw)
  dev.off()
  
  # Save Elbow plots
  pdf(paste0(path, "elbow_harmony.pdf"), width = 8, height = 6)
  elbow_harm <- ElbowPlot(seurat_obj, reduction='harmony')
  print(elbow_harm)
  dev.off()

  pdf(paste0(path, "elbow_pca.pdf"), width = 8, height = 6)
  elbow_pca <- ElbowPlot(seurat_obj, reduction='pca')
  print(elbow_pca)
  dev.off()

  # Perform K-nearest neighbor (KNN) graph using harmony embeddings
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_umap, reduction='harmony')
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_umap, reduction = "harmony")

  # Save UMAP lanes plot
  pdf(paste0(path, "umap_lanes.pdf"), width = 8, height = 6)
  umap_lanes <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = .1)
  print(umap_lanes)
  dev.off()
  
  # Cluster cells
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = "leiden")

  # Save UMAP clusters plot
  pdf(paste0(path, "umap_clusters.pdf"), width = 8, height = 6)
  umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .1)
  print(umap_clusters)
  dev.off()
  
  # Return the updated Seurat object
  return(seurat_obj)
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

create_feature_scatter_plot <- function(obj, feature1, feature2, save = TRUE, file_name = NULL, path = "plot_outputs/") {

  plot <- FeatureScatter(object = obj, feature1 = feature1, feature2 = feature2)
  
  # If save is TRUE, save the plot
  if (save) {
    # If file_name is not provided, use default naming scheme
    if (is.null(file_name)) {
      file_name <- paste0(feature1, "_vs_", feature2, "_scatter")
    }
    
    pdf(paste0(path, file_name, ".pdf"), width = 8, height = 6)
    print(plot)
    dev.off()
  }
  
  return(plot)
}

combine_feature_plots <- function(seurat_objs_list, feature_set1, feature_set2 = NULL, same_feature_set = TRUE, file_name_base = "post_SoupX_plot", path = "plot_outputs/") {
  
  # If same_feature_set is TRUE, use feature_set1 for both plots
  if (same_feature_set) {
    feature_set2 <- feature_set1
  }
  
  # Create individual plots and save them as PDFs
  plot1 <- create_feature_scatter_plot(seurat_objs_list[[1]]$seu, feature_set1$feature1, feature_set1$feature2, save = TRUE, file_name = paste0(file_name_base, "1"), path = path)
  plot2 <- create_feature_scatter_plot(seurat_objs_list[[2]]$seu, feature_set2$feature1, feature_set2$feature2, save = TRUE, file_name = paste0(file_name_base, "2"), path = path)
  
  # Add plots together and save as a new PDF
  pdf(paste0(path, file_name_base, "_combined.pdf"), width = 8, height = 6)
  combined_plot <- plot1 + plot2
  print(combined_plot)
  dev.off()
  
  # Delete the individual plot files
  file.remove(paste0(path, file_name_base, "1.pdf"))
  file.remove(paste0(path, file_name_base, "2.pdf"))
  
  return(combined_plot)
}