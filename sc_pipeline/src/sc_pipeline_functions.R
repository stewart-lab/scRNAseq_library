split_sce_list_by_sample <- function(filtered_sce_list) {
  # Removed parallelization, using simple map
  identifiers <- map_chr(filtered_sce_list, function(x) {
    ident_value <- colData(x)$ident[1]
    sub("_lane.*", "", ident_value)
  })
  unique_identifiers <- unique(identifiers)
  # Direct map instead of future_map
  split_lists <- map(unique_identifiers, function(condition) {
    filtered_sce_list[identifiers == condition]
  }) %>% setNames(unique_identifiers)
  return(split_lists)
}

run_scDblFinder_and_merge <- function(sce_list, path, save_plot = TRUE, file_name = "after_dbl_removal_and_merge") {
  set.seed(1234)
  # Removed parallelization, use map instead of future_map
  processed_results <- map(seq_along(sce_list), function(i) {
    sce <- scDblFinder(sce_list[[i]])
    dbl_table <- as.data.frame(table(call = sce$scDblFinder.class))
    sce_filtered <- sce[, sce$scDblFinder.class == "singlet"]
    seurat_obj <- as.Seurat(sce_filtered, data = NULL)
    list(dbl_table = dbl_table, seurat_obj = seurat_obj)
  })
  
  dbl_tables <- map(processed_results, "dbl_table")
  seurat_objects <- map(processed_results, "seurat_obj")
  
  merged_dbl_table <- do.call(rbind, dbl_tables)
  write.table(merged_dbl_table, file = file.path(path, "merged_doublet_table.txt"),
              sep = "\t", row.names = TRUE, quote = FALSE)
  
  merged_seurat_obj <- Reduce(function(x, y) merge(x, y), seurat_objects)
  
  if (save_plot) {
    create_feature_scatter_plot(
      obj = merged_seurat_obj,
      feature1 = "nCount_RNA",
      feature2 = "nFeature_RNA",
      file_name = file_name,
      path = path
    )
  }
  
  return(merged_seurat_obj)
}

filter_cells <- function(seurat_obj, sample_name, path, save_plots = TRUE) {
  params <- config$filter_cells
  species <- NULL
  for (sample in names(config$fastq_alignment)) {
    sample_details <- config$fastq_alignment[[sample]]
    if (sample_details$NAME == sample_name) {
      species <- sample_details$species
      break
    }
  }
  if (is.null(species)) stop("Species is NULL. Please check your config file.")
  
  # Removed parallelization in mitochondrial feature checks
  if (species == "pig") {
    mt.list <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
    valid_features <- map_lgl(mt.list, ~ .x %in% rownames(seurat_obj))
    mt.list <- mt.list[valid_features]
    percent_mt <- PercentageFeatureSet(seurat_obj, features = mt.list, assay = "RNA")
  } else if (species == "human") {
    percent_mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "RNA")
  } else {
    stop("Species not recognized: please input 'pig' or 'human'")
  }
  
  seurat_obj[["percent.mt"]] <- percent_mt
  
  if (save_plots) {
    # Direct call instead of future
    create_feature_scatter_plot(seurat_obj, "nCount_RNA", "percent.mt",
                                file_name = "percent_mt_unfiltered",
                                save = TRUE, path = path)
  }
  
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > params$lower.nFeature &
             nFeature_RNA < params$upper.nFeature &
             percent.mt < params$max.percent.mt
  )
  
  if (save_plots) {
    create_feature_scatter_plot(seurat_obj, "nCount_RNA", "percent.mt",
                                file_name = "percent_mt_filtered",
                                save = TRUE, path = path)
  }
  
  return(seurat_obj)
}

create_feature_scatter_plot <- function(obj, feature1, feature2, file_name, save = TRUE, path) {
  p <- FeatureScatter(obj, feature1 = feature1, feature2 = feature2)
  if (save) {
    pdf(file.path(path, paste0(file_name, ".pdf")))
    print(p)
    dev.off()
  }
  return(p)
}

normalize_data <- function(seurat_obj, path) {
  min_size <- config$normalize_data$min_size
  min_mean <- config$normalize_data$min_mean
  feature <- config$normalize_data$feature
  
  sce <- as.SingleCellExperiment(seurat_obj)
  clusters <- quickCluster(sce, use.ranks = FALSE, min.size = min_size)
  sce <- computeSumFactors(sce, clusters = clusters, min.mean = min_mean)
  sce <- logNormCounts(sce)
  
  seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = logcounts(sce))
  seurat_obj$sizeFactors <- sizeFactors(sce)
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  
  vin_pre <- VlnPlot(seurat_obj, feature, slot = "counts")
  vin_post <- VlnPlot(seurat_obj, feature, slot = "data")
  
  pdf(file.path(path, "violin_pre_norm.pdf"), width = 8, height = 6)
  print(vin_pre)
  dev.off()
  
  pdf(file.path(path, "violin_post_norm.pdf"), width = 8, height = 6)
  print(vin_post)
  dev.off()
  
  return(seurat_obj)
}

feature_selection <- function(seurat_obj) {
  n_features <- config$feature_selection$n_features
  analysis_type <- config$feature_selection$analysis_type
  
  if (analysis_type == "Seurat") {
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = n_features)
  } else if (analysis_type == "Scry") {
    m <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
    devi <- scry::devianceFeatureSelection(m)
    dev_ranked_genes <- rownames(seurat_obj)[order(devi, decreasing = TRUE)]
    topdev <- head(dev_ranked_genes, n_features)
    VariableFeatures(seurat_obj) <- topdev
    seurat_obj <- UpdateSeuratObject(seurat_obj)
  } else {
    stop("Invalid analysis_type. Please choose 'Seurat' or 'Scry'.")
  }
  
  return(seurat_obj)
}


scale_data <- function(seurat_obj, path, sample_name = NULL) {
  vars.2.regress <- config$scale_data$vars.2.regress
  marker.path.s <- config$scale_data$marker.path.s
  marker.path.g2m <- config$scale_data$marker.path.g2m
  
  if (is.null(sample_name)) sample_name <- unique(seurat_obj$orig.ident)[1]
  
  species <- NULL
  for (sample in names(config$fastq_alignment)) {
    sample_details <- config$fastq_alignment[[sample]]
    if (sample_details$NAME == sample_name) {
      species <- sample_details$species
      break
    }
  }
  
  if (is.null(species)) {
    warning("Species not found in config, defaulting to 'human'")
    species <- "human"
  }
  
  all.genes <- rownames(seurat_obj)
  
  if (vars.2.regress == "cell.cycle") {
    cell.cycle.markers.s <- read.csv2(marker.path.s, sep = "\t", header = TRUE, row.names = 1)
    cell.cycle.markers.g2m <- read.csv2(marker.path.g2m, sep = "\t", header = TRUE, row.names = 1)
    varslist <- c(cell.cycle.markers.s, cell.cycle.markers.g2m)
    
    if (species == "human") {
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
    } else if (species == "pig") {
      s.genes <- varslist[[4]]$pig.gene.name
      g2m.genes <- varslist[[8]]$pig.gene.name
    } else {
      warning("Unsupported species, defaulting to human cell cycle genes")
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
    }
    
    missing_genes <- map(g2m.genes, ~ if(!(.x %in% rownames(seurat_obj))) .x else NULL) %>% compact()
    if (length(missing_genes) > 0) {
      warning("Missing genes: ", paste(missing_genes, collapse = ", "))
    }
    
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
    # Removed future calls, run directly
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    pdf(file.path(path, "pca_before_cc_regression.pdf"), width = 8, height = 6)
    print(DimPlot(seurat_obj, group.by = "Phase"))
    dev.off()

    seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score","G2M.Score"), features = all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    pdf(file.path(path,"pca_after_cc_regression.pdf"),width=8,height=6)
    print(DimPlot(seurat_obj, group.by="Phase"))
    dev.off()
  } else {
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  }
  
  return(seurat_obj)
}


run_and_visualize_pca <- function(seurat_obj, path) {
  top_n_dims <- config$run_and_visualize_pca$top_n_dims
  heatmap_dims <- 1:config$run_and_visualize_pca$heatmap_dims
  num_cells <- config$run_and_visualize_pca$num_cells
  dims <- 1:config$run_and_visualize_pca$dims
  num_replicate <- config$run_and_visualize_pca$num.replicate
  
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  pdf(file.path(path,"top_n_dims_with_genes.pdf"),width=8,height=6)
  print(VizDimLoadings(seurat_obj,dims=1:top_n_dims,reduction="pca"))
  dev.off()
  
  pdf(file.path(path,"pca_scatter_plot.pdf"),width=8,height=6)
  print(DimPlot(seurat_obj,reduction="pca"))
  dev.off()
  
  pdf(file.path(path,"pca_heat_map.pdf"),width=8,height=6)
  print(DimHeatmap(seurat_obj,nfeatures=5,dims=heatmap_dims,cells=num_cells,balanced=TRUE,fast=FALSE,combine=TRUE))
  dev.off()
  
  pdf(file.path(path,"elbow_pca.pdf"),width=8,height=6)
  elbow_pca <- ElbowPlot(seurat_obj,reduction="pca")
  print(elbow_pca)
  dev.off()
  
  if (num_replicate == "NA") {
    print("jackstraw not run")
  } else {
    seurat_obj <- JackStraw(seurat_obj,num.replicate=num_replicate)
    seurat_obj <- ScoreJackStraw(seurat_obj,dims=dims)
    
    pdf(file.path(path,"jack_straw.pdf"),width=8,height=6)
    jack_straw <- JackStrawPlot(seurat_obj,dims=dims)
    print(jack_straw)
    dev.off()
  }
  
  return(seurat_obj)
}


perform_batch_correction <- function(seurat_obj,path) {
  dims.use <- 1:config$perform_batch_correction$dims.use
  max_iter <- config$perform_batch_correction$max_iter
  
  pdf(file.path(path,"batch_uncorrected_pca.pdf"),width=8,height=6)
  p1_pre <- DimPlot(seurat_obj,reduction="pca",pt.size=.1,group.by="orig.ident")
  p2_pre <- VlnPlot(seurat_obj,features="PC_1",group.by="orig.ident",pt.size=.1)
  print(p1_pre+p2_pre)
  dev.off()
  
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars="orig.ident",dims.use=dims.use,max.iter.harmony=max_iter)
  harmony_embeddings <- Embeddings(seurat_obj,"harmony")
  
  pdf(file.path(path,"batch_corrected_pca.pdf"),width=8,height=6)
  p1_post <- DimPlot(seurat_obj,reduction="harmony",pt.size=.1,group.by="orig.ident")
  p2_post <- VlnPlot(seurat_obj,features="harmony_1",group.by="orig.ident",pt.size=.1)
  print(p1_post+p2_post)
  dev.off()
  
  return(list(seurat_obj=seurat_obj,harmony_embeddings=harmony_embeddings))
}

run_umap <- function(seurat_obj, path) {
  # Debug: Print initial function entry
  message("DEBUG: Entering run_umap function")
  
  # Debug: Print config object structure
  message("DEBUG: Config structure:")
  message(str(config))
  
  # Debug: Print raw values before conversion
  message("DEBUG: Raw config values:")
  message("dims_umap: ", paste(capture.output(str(config$run_umap$dims_umap)), collapse="\n"))
  message("umap.method: ", paste(capture.output(str(config$run_umap$umap.method)), collapse="\n"))
  
  # Try getting values with error catching
  tryCatch({
    dims_umap <- as.numeric(1:config$run_umap$dims_umap)
    message("DEBUG: dims_umap converted successfully: ", paste(dims_umap, collapse=", "))
  }, error = function(e) {
    message("ERROR in dims_umap conversion: ", e$message)
    stop(e)
  })
  
  tryCatch({
    umap.method <- as.character(config$run_umap$umap.method)
    message("DEBUG: umap.method converted successfully: ", umap.method)
  }, error = function(e) {
    message("ERROR in umap.method conversion: ", e$message)
    stop(e)
  })
  
  # Debug: Print object types
  message("DEBUG: Object types:")
  message("dims_umap type: ", typeof(dims_umap))
  message("umap.method type: ", typeof(umap.method))
  
  # Determine reduction method with debug
  message("DEBUG: Checking for harmony embeddings")
  umap.red <- if ("harmony" %in% names(Embeddings(seurat_obj))) {
    message("DEBUG: Using harmony reduction")
    "harmony"
  } else {
    message("DEBUG: Using pca reduction")
    "pca"
  }
  
  # Run UMAP based on method with debug
  tryCatch({
    message("DEBUG: Starting UMAP with method: ", umap.method)
    
    if (umap.method == "uwot") {
      message("DEBUG: Using uwot method")
      seurat_obj <- RunUMAP(seurat_obj, 
                           dims = dims_umap, 
                           umap.method = "uwot",
                           reduction = umap.red,
                           group.by = "orig.ident")
    } else {
      message("DEBUG: Using umap-learn method")
      library(reticulate)
      
      # Debug Python environment
      message("DEBUG: Python configuration:")
      message(str(py_config()))
      
      message("DEBUG: Importing UMAP")
      py_run_string("import umap.umap_ as umap")
      py_run_string("from umap import UMAP")
      
      seurat_obj <- RunUMAP(seurat_obj, 
                           dims = dims_umap, 
                           umap.method = "umap-learn",
                           reduction = umap.red,
                           group.by = "orig.ident",
                           metric = "correlation")
    }
    
    message("DEBUG: UMAP completed successfully")
    
  }, error = function(e) {
    message("DEBUG: Error in UMAP execution: ", e$message)
    warning("UMAP failed, falling back to uwot method")
    seurat_obj <<- RunUMAP(seurat_obj, 
                          dims = dims_umap, 
                          umap.method = "uwot",
                          reduction = umap.red,
                          group.by = "orig.ident")
  })
  
  # Generate plot with debug
  message("DEBUG: Generating UMAP plot")
  tryCatch({
    pdf(paste0(path, "umap_plot.pdf"), width = 8, height = 6)
    print(DimPlot(seurat_obj, reduction = "umap"))
    dev.off()
    message("DEBUG: Plot generated successfully")
  }, error = function(e) {
    message("ERROR in plot generation: ", e$message)
  })
  
  message("DEBUG: Exiting run_umap function")
  return(seurat_obj)
}

# Wrapper function for debugging
debug_run_umap <- function(seurat_obj, path) {
  tryCatch({
    message("DEBUG: Starting UMAP processing")
    result <- run_umap(seurat_obj, path)
    message("DEBUG: UMAP processing completed")
    return(result)
  }, error = function(e) {
    message("ERROR: ", e$message)
    message("Stack trace:")
    print(sys.calls())
    stop(e)
  })
}

perform_clustering <- function(seurat_obj, path) {
  resolution <- config$perform_clustering$resolution
  algorithm <- config$perform_clustering$algorithm
  reduction <- config$perform_clustering$reduction
  dims_snn <- 1:config$perform_clustering$dims_snn

  # Check if Harmony embeddings exist
  batch_corrected <- "harmony" %in% names(Reductions(seurat_obj))
  if (!batch_corrected && reduction == "harmony") {
    message("Batch correction was skipped. Updating reduction to 'pca'.")
    reduction <- "pca"
  }

  # Perform KNN
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_snn, reduction = reduction)
  # Import Python's leidenalg
  leidenalg <- reticulate::import("leidenalg")
  igraph <- reticulate::import("igraph")
  # Cluster using appropriate method
  # if (algorithm == "leiden") {
  #   message("Using Python leidenalg for clustering")
  #   seurat_obj <- run_leiden_clustering(seurat_obj, resolution = resolution)
  #   # Set the active clusters to leiden_clusters
  #   Idents(seurat_obj) <- "leiden_clusters"
  # } else {
  message("Using Seurat clustering algorithm: ", algorithm)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = algorithm)
  #}

  # Generate plots
  pdf(paste0(path, "umap_lanes.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = .5))
  dev.off()

  pdf(paste0(path, "umap_clusters.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",label = TRUE, pt.size = .5))
  dev.off()

  return(seurat_obj)
}



process_sample <- function(sample_name, sample_data, output_base_dir, config) {
  # Initialize an environment for storing intermediate results
  env <- new.env()
  
  # Make output directory for each sample
  sample_output_dir <- file.path(output_base_dir, sample_name)
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  sample_output_dir <- paste0(sample_output_dir, "/")
  
  # 1. Doublet finder and lane merging
  env$lane_and_merged_seurat_obj <- run_scDblFinder_and_merge(sample_data, sample_output_dir)
  
  # 2. Mitochondrial gene filtering
  env$seurat_obj_mt_filtered <- filter_cells(env$lane_and_merged_seurat_obj, sample_name, sample_output_dir)
  rm(list = "lane_and_merged_seurat_obj", envir = env); env$lane_and_merged_seurat_obj <- NULL
  
  # Perform garbage collection after each major step to manage memory
  gc(full = TRUE)
  
  # Normalize, select features, scale, and run PCA
  env$normalized_seurat_obj <- normalize_data(env$seurat_obj_mt_filtered, sample_output_dir)
  rm(list = "seurat_obj_mt_filtered", envir = env); env$seurat_obj_mt_filtered <- NULL
  env$feature_selected_seurat_obj <- feature_selection(env$normalized_seurat_obj)
  rm(list = "normalized_seurat_obj", envir = env); env$normalized_seurat_obj <- NULL
  env$scaled_seurat_obj <- scale_data(env$feature_selected_seurat_obj, sample_output_dir)
  rm(list = "feature_selected_seurat_obj", envir = env); env$feature_selected_seurat_obj <- NULL
  env$dim_reduced_seurat_obj <- run_and_visualize_pca(env$scaled_seurat_obj, sample_output_dir)
  rm(list = "scaled_seurat_obj", envir = env); env$scaled_seurat_obj <- NULL
  
  # Batch correction if needed
  if (length(unique(env$dim_reduced_seurat_obj$orig.ident)) > 1) {
    batchList <- perform_batch_correction(env$dim_reduced_seurat_obj, sample_output_dir)
    env$batch_corrected_obj <- batchList[["seurat_obj"]]
    saveRDS(env$batch_corrected_obj, file = paste0(sample_output_dir, sample_name, "_batchcorr_seurat_obj.rds"))
  } else {
    message("Skipping batch correction as 'orig.ident' has only one level.")
    env$batch_corrected_obj <- env$dim_reduced_seurat_obj
  }
  rm(list = "dim_reduced_seurat_obj", envir = env); env$dim_reduced_seurat_obj <- NULL
  
  # Process based on species consistency
  if (species_are_all_same(config)) {
    message("Species are consistent across all samples.")
    
    # Determine whether the batch_corrected_obj is a list and run UMAP accordingly
    if (is.list(env$batch_corrected_obj)) {
      final_obj <- process_consistent_species(env$batch_corrected_obj$seurat_obj, sample_output_dir, config, sample_name)
    } else {
      final_obj <- process_consistent_species(env$batch_corrected_obj, sample_output_dir, config, sample_name)
    }
    
    # Return the final processed object
    return(final_obj)
  } else {
    message("Species are not consistent. Returning batch corrected object for orthologous analysis.")
    # Return the batch corrected object for later orthologous analysis
    return(env$batch_corrected_obj)
  }
}



find_differentially_expressed_features <- function(seurat_obj, path) {
  # Get parameters from the config file
  min_pct <- config$find_differentially_expressed_features$min_pct
  logfc_threshold <- config$find_differentially_expressed_features$logfc_threshold
  top_n <- config$find_differentially_expressed_features$top_n

  # Find all markers
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)

  # Extract top_n markers
  markers %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC) -> topMarkers

  # Write out top_n markers
  write.table(topMarkers, file = paste0(path, "/seurat_obj.DE.markers.top", top_n, "_S2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

  # Define the top 4 features for FeaturePlot based on topMarkers
  top_features <- topMarkers$gene[1:4]

  # UMAP plots
  plot1 <- UMAPPlot(seurat_obj, group.by = "orig.ident")
  plot2 <- UMAPPlot(seurat_obj, label = T)
  plot3 <- FeaturePlot(seurat_obj, features = top_features, ncol = 2, pt.size = 0.1)

  # Combine plots
  combined_plot <- ((plot1 / plot2) | plot3) + plot_layout(width = c(1, 2))

  # Save plot
  pdf(file = paste0(path, "/combined_plot.pdf"), width = 8, height = 6)
  print(combined_plot)
  dev.off()

  # Heatmap of top_n markers
  heatmap <- DoHeatmap(seurat_obj, features = topMarkers$gene) + NoLegend()

  # Save heatmap
  pdf(file = paste0(path, "/heatmap_top_n_markers.pdf"))
  print(heatmap)
  dev.off()

  # Write out all markers
  write.table(markers, file = paste0(path, "/seurat_obj.DE.markers_S2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

  return(list(markers = markers, topMarkers = topMarkers))
}

analyze_known_markers <- function(seurat_obj, de_results, output_path) {
  # read in known GAMM retinoid markers
  known_markers_path <- config$score_and_plot_markers$known_markers_path
  known.markers <- read.csv2(known_markers_path, sep = "\t", header = TRUE)
  de.markers <- de_results[[1]]
  # match with any DE markers from data by merging dataframes
  marker_df <- merge(de.markers, known.markers, by = "gene")
  # write out marker df with known DE markers
  write.table(marker_df, file = paste0(output_path, "seurat_obj.knownDE.markers_S2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

  # first get unique cell type vector
  cell.types <- unique(marker_df$Cell.type)
  print(cell.types)

  # check gene number between known markers and DE known markers
  genesk <- unique(known.markers$gene)
  k <- length(genesk)
  genesDEk <- unique(marker_df$gene)
  DEk <- length(genesDEk)
  percent.markers.de <- (DEk / k) * 100
  percent.markers.de

  # check a specific marker
  df <- subset(marker_df, Cell.type == "Synaptic marker", select = c("gene"))
  length(df$gene)
  df2 <- subset(known.markers, Cell.type == "Synaptic marker", select = c("gene"))
  length(df2$gene)

  # subset all and plot using for loop
  for (i in 1:length(cell.types)) {
    new_df <- subset(marker_df, Cell.type == cell.types[i], select = c("gene", "Cell.type", "cluster"))
    new_vec <- unique(as.vector(new_df$gene))

    pdf(paste0(output_path, cell.types[i], "_featureplot.pdf"), width = 8, height = 6, bg = "white")
    # umap plot highlighting gene expression
    print(FeaturePlot(seurat_obj, features = new_vec))
    dev.off()

    pdf(paste0(output_path, cell.types[i], "_dotplot.pdf"), width = 8, height = 6, bg = "white")
    # expression dot plot
    dot.plot <- DotPlot(object = seurat_obj, features = new_vec)
    print(dot.plot + labs(title = cell.types[i]))
    dev.off()
  }
}

score_and_plot_markers <- function(seurat_obj, output_path) {
  known_markers_path <- config$score_and_plot_markers$known_markers_path
  known_markers <- config$score_and_plot_markers$known_markers
  top_n_markers <- config$score_and_plot_markers$top_n_markers
  cluster_type <- config$score_and_plot_markers$cluster_type
  pairwise <- config$score_and_plot_markers$pairwise
  logFC_thresh <- config$score_and_plot_markers$logFC_thresh
  auc_thresh <- config$score_and_plot_markers$auc_thresh


  sce_obj <- as.SingleCellExperiment(seurat_obj)

  # Score markers
  marker_info <- score_markers(sce_obj, cluster_type)

  # Read in known markers
  known.markers.df <- if (known_markers) {
    read.csv2(known_markers_path, sep = "\t", header = TRUE, row.names = 1)
  } else {
    NULL
  }

  annot_df <- data.frame(Cluster = integer(), Cell.type = character())

  # Determine clusters based on cluster_type
  clusters <- get_clusters(seurat_obj, cluster_type)

  # Process each cluster
  for (i in 1:length(clusters)) {
    # Process pairwise comparisons
    if (pairwise) {
      process_pairwise_comparisons(clusters, i, marker_info, output_path, logFC_thresh, auc_thresh, known.markers.df, seurat_obj)
    } else {
      print("No pairwise comparison")
    }

    # Process top genes for each cluster
    clust <- as.data.frame(marker_info[[clusters[i]]])
    annot_df <- process_top_genes(clust, clusters, i, known.markers.df, output_path, seurat_obj, annot_df)
  }

  return(annot_df)
}

score_markers <- function(sce_obj, cluster_type) {
  # Score markers based on cluster type
  if (cluster_type %in% c("seurat_clusters", "orig.ident")) {
    marker_field <- cluster_type
    marker_info <- scoreMarkers(sce_obj, sce_obj@colData@listData[[marker_field]], full.stats = TRUE)
  } else {
    stop("Invalid cluster_type. Please choose 'seurat_clusters' or 'orig.ident'.")
  }
  return(marker_info)
}

get_clusters <- function(seurat_obj, cluster_type) {
  if (cluster_type %in% c("seurat_clusters", "orig.ident")) {
    return(unique(seurat_obj@meta.data[[cluster_type]]))
  } else {
    stop("Invalid cluster_type. Please choose 'seurat_clusters' or 'orig.ident'.")
  }
}

process_top_genes <- function(clust, clusters, i, known.markers.df, output_path, seurat_obj, annot_df) {
  # Extract relevant configuration settings
  top_n_markers <- config$score_and_plot_markers$top_n_markers
  logFC_thresh <- config$score_and_plot_markers$logFC_thresh
  auc_thresh <- config$score_and_plot_markers$auc_thresh
  known_markers <- config$score_and_plot_markers$known_markers

  # Order by median Cohen's d and subset
  ordered <- subset(clust[order(clust$median.logFC.cohen, decreasing = TRUE), ], median.logFC.cohen > logFC_thresh & median.AUC > auc_thresh)
  top100 <- head(ordered, n = top_n_markers)

  # Create a subdirectory for Top100 DE genes
  top100_dir <- file.path(output_path, "Top100_DE_Genes")
  if (!dir.exists(top100_dir)) {
    dir.create(top100_dir, recursive = TRUE)
  }

  # Write out top100 genes
  write.table(top100,
    file = file.path(top100_dir, paste0("Top100DEgenes_clust_", clusters[i], ".txt")),
    quote = FALSE, sep = "\t", row.names = TRUE
  )

  # Process known markers
  return(process_known_markers(top100, known_markers, known.markers.df, clusters, i, output_path, top_n_markers, seurat_obj, annot_df))
}


process_pairwise_comparisons <- function(clusters, i, marker.info, output_path, logFC_thresh, auc_thresh, known_markers_df, seurat_obj) {
  # Extract cluster marker information
  clust <- as.data.frame(marker.info[[clusters[i]]])

  # Prepare the data for pairwise comparison
  clust_tbl <- rownames_to_column(clust, var = "gene") %>% as_tibble()
  df_clust0 <- clust_tbl %>% dplyr::select(starts_with("full.logFC.cohen"))
  df_clust <- clust_tbl %>% dplyr::select(c(
    "gene", starts_with("full.logFC.cohen"),
    starts_with("full.AUC")
  ))

  # Create a subdirectory for pairwise DE genes
  pairwise_de_path <- file.path(output_path, paste0("Pairwise_DE_Cluster_", clusters[i]))
  if (!dir.exists(pairwise_de_path)) {
    dir.create(pairwise_de_path, recursive = TRUE)
  }

  # Iterate through each pairwise comparison
  for (j in colnames(df_clust0)) {
    print(j)
    k <- str_split_1(j, "cohen.")
    l <- paste0("full.AUC.", k[2])

    # Subset based on logFC threshold and AUC threshold
    df_clust1 <- subset(df_clust, df_clust[j] > logFC_thresh & df_clust[l] > auc_thresh)
    df_clust1 <- df_clust1[order(df_clust1[j], decreasing = TRUE), ]

    # Write table of all DE genes for the comparison
    write.table(df_clust1,
      file = file.path(pairwise_de_path, paste0("DEgenes_", clusters[i], "_vs_", j, ".txt")),
      quote = FALSE, sep = "\t", row.names = TRUE
    )

    # Merge with known markers and check if empty
    df_clust2 <- merge(df_clust1, known_markers_df, by.x = "gene", by.y = "row.names")
    if (nrow(df_clust2) == 0) {
      print(paste0("This data frame is empty: ", clusters[i], ".vs_", j))
    } else {
      # Write table with known DE markers
      write.table(df_clust2,
        file = file.path(pairwise_de_path, paste0("KnownDEgenes_", clusters[i], "_vs_", j, ".txt")),
        quote = FALSE, sep = "\t", row.names = TRUE
      )

      # Make feature plot of genes
      new_vec0 <- unique(as.vector(df_clust2$gene))
      pdf(file.path(pairwise_de_path, paste0(clusters[i], "_vs_", j, "_featureplot.pdf")), bg = "white")
      print(FeaturePlot(seurat_obj, features = new_vec0), label = TRUE)
      dev.off()
    }
  }
}

process_known_markers <- function(top100, known_markers_flag, known_markers_df, clusters, i, output_path, top_n_markers, seurat_obj, annot_df) {
  annot_type <- config$process_known_markers$annot_type
  n_rank <- config$process_known_markers$n_rank
  if (known_markers_flag) {
    marker_df <- merge(top100, known_markers_df, by = "row.names")

    if (nrow(marker_df) == 0) {
      print(paste0("This data frame is empty: ", clusters[i]))
      new_row <- data.frame(Cluster = clusters[i], Cell.type = "unknown")
      annot_df <- rbind(annot_df, new_row)
    } else {
      subdirectory_path <- file.path(output_path, "Known_DE_Markers")
      if (!dir.exists(subdirectory_path)) {
        dir.create(subdirectory_path, recursive = TRUE)
      }
      # Write out marker dataframe with known DE markers
      write.table(marker_df,
        file = paste0(subdirectory_path, "/KnownDE.markers_clust_", clusters[i], ".txt"),
        quote = FALSE, sep = "\t", row.names = FALSE
      )

      # Subset data
      new_df <- marker_df[, c("Row.names", "rank.logFC.cohen", "Cell.type")]
      new_vec <- unique(as.vector(new_df$Row.names))

      # Get top ranked
      rank <- n_rank + 1
      new_df.ordered <- new_df[order(new_df$rank.logFC.cohen), ]

      if (annot_type == "manual"){
        new_df.ordered <- subset(new_df.ordered, rank.logFC.cohen < rank)
        new_vec2 <- unique(as.vector(new_df.ordered$Row.names))

        if (identical(new_vec2, character(0))) {
          print(paste0("This vector does not have any ranks in top ", top_n_markers, ": ", clusters[i]))
          new_row <- data.frame(Cluster = clusters[i], Cell.type = "unknown")
          annot_df <- rbind(annot_df, new_row)
        } else {
        # UMAP plot highlighting gene expression
          pdf(paste0(output_path, clusters[i], "_featureplot_top", top_n_markers, "ranks.pdf"), bg = "white")
          print(FeaturePlot(seurat_obj, features = new_vec2), label = TRUE)
         dev.off()
          allcelltypes <- unique(as.vector(new_df.ordered$Cell.type))
          result_string <- paste(allcelltypes, collapse = "-")
          new_row <- data.frame(Cluster = clusters[i], Cell.type = result_string)
          annot_df <- rbind(annot_df, new_row)
        }
      } else if (annot_type == "d120"| annot_type == "d40"){
        new_vec2 <- unique(as.vector(new_df.ordered$Row.names))
        cell_types <- unique(as.vector(known_markers_df[,"Cell.type"]))
        # create empty list to store cell types
        cell_type_list <- c()
        cell_type_list1 <- c()
        for (j in 1:length(cell_types)){
          cell_type <- cell_types[j]
          #print(cell_type)
          genes_df <- subset(known_markers_df, Cell.type == cell_type)
          #print(colnames(genes_df))
          genes <- unique(rownames(genes_df))
          #print(genes)
          count = 0
          for (k in 1:length(new_vec2)){
            gene <- new_vec2[k]
            if (gene %in% genes){
              count <- count + 1
            } else {
                next
            }
          }
          if (count >= 2){
            new_cell_type <- cell_type
            #print(new_cell_type)
          } else if (count >= 1 && cell_type == "Cone") {
            # check Pan PRs
            count2 <- 0
            genes_df <- subset(known_markers_df, Cell.type == "Pan PR")
            genes <- unique(rownames(genes_df))
            for (k in 1:length(new_vec2)){
              gene <- new_vec2[k]
              if (gene %in% genes){
                count2 <- count2 + 1
              } else {
                next
              }
            }
            if (count2 >= 2){
              new_cell_type <- "Cone"
            } else {
              new_cell_type <- "NA"
            }
          } else if (count >= 1 && cell_type == "Ganglion Cell") {
            # check Amacrine-Ganglion
            count3 <- 0
            genes_df <- subset(known_markers_df, Cell.type == "Amacrine-Ganglion") # nolint
            genes <- unique(rownames(genes_df))
            for (k in 1:length(new_vec2)){
              gene <- new_vec2[k]
              if (gene %in% genes){
                count3 <- count3 + 1
              } else {
                next
              }
            }
            if (count3 >= 2){
              new_cell_type <- "Ganglion Cell"
            } else if (annot_type == "d40" && count3 >= 1){
              new_cell_type <- "Ganglion Cell"
            } else {
              new_cell_type <- "NA"
            } 
          } else if (count >= 1 && cell_type == "Amacrine Cell") {
            # check Amacrine-Ganglion
            count4 <- 0
            genes_df <- subset(known_markers_df, Cell.type == "Amacrine-Ganglion") # nolint
            genes <- unique(rownames(genes_df))
            for (k in 1:length(new_vec2)){
              gene <- new_vec2[k]
              if (gene %in% genes){
                count4 <- count4 + 1
              } else {
                next
              }
            }
            if (count4 >= 2){
              new_cell_type <- "Amacrine Cell"
            } else {
              new_cell_type <- "NA"
            }
          } else {
            new_cell_type <- "NA"
          }
          if (count >= 1) {
            new_cell_type1 <- cell_type
            cell_type_list1 <- c(cell_type_list1, new_cell_type1)
          }
          cell_type_list <- c(cell_type_list, new_cell_type)
        }
        if ("Cone" %in% cell_type_list && "Pan PR" %in% cell_type_list) {
          cell_type_list <- c("Cone")
          print("Cone!")
        } else if ("Rod" %in% cell_type_list && "Pan PR" %in% cell_type_list) {
          print("Rod!")
          cell_type_list <- c("Rod")
        } else if ("Amacrine Cell" %in% cell_type_list && 
                    "Amacrine-Ganglion" %in% cell_type_list) {
          cell_type_list <- c("Amacrine Cell")
          print("Amacrine Cell!")
        } else if ("Ganglion Cell" %in% cell_type_list && 
                    "Amacrine-Ganglion" %in% cell_type_list) {
          cell_type_list <- c("Ganglion Cell")
          print("Ganglion Cell!")
        } else if (all(is.na(c("NA")) %in% names(cell_type_list))) {
          cell_type_list <- c("unknown")
        } else {
          cell_type_list <- cell_type_list[cell_type_list != "NA"]
        }
        if ("unknown" %in% cell_type_list | length(cell_type_list) == 0 | 
            length(cell_type_list) > 1) {
          if (length(unique(cell_type_list1)) > 2 && annot_type == "d40") {
            cell_type_list <- c("Retinal Prog")
          } else if (length(unique(cell_type_list)) == 2 && annot_type == "d120") {
            cell_type_list <- cell_type_list
          } else {
            cell_type_list <- c("unknown")
          }
        } else {
          cell_type_list <- cell_type_list
        }
        result_string <- paste(cell_type_list, collapse = "-")
        print(paste0("final cell type ", clusters[i], " ",result_string))
        new_row <- data.frame(Cluster = clusters[i], Cell.type = result_string)
        annot_df <- rbind(annot_df, new_row)
        } else {
          print("Need to set annotation type in config")
        }
      }
  } else {
    print("No known marker set")
  }
  return(annot_df)
}


annotate_clusters_and_save <- function(seurat_obj, new_cluster_ids, output_path = output) {
  # make sure idents are set to seurat_clusters
  Idents(seurat_obj) <- "seurat_clusters"
  # Rename the clusters based on the new IDs
  names(new_cluster_ids) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
  # store names as celltype
  seurat_obj$CellType <- Idents(seurat_obj)
  # Generate and plot the UMAP plot

  pdf(paste0(output_path, "labeled-clusters.pdf"), bg = "white")
  print(DimPlot(seurat_obj, reduction = "umap", group.by = 'CellType', label = TRUE, pt.size = 0.5))
  dev.off()
  # Save the Seurat object
  saveRDS(seurat_obj, file = paste0(output_path, "seurat_obj_labeled.rds"))

  return(seurat_obj)
}

create_feature_scatter_plot <- function(obj, feature1, feature2, save = TRUE, file_name = NULL, path = output) {
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

combine_feature_plots <- function(seurat_objs_list, path, feature_set1, feature_set2 = NULL, same_feature_set = TRUE, file_name_base = "post_SoupX_plot") {
  # If same_feature_set is TRUE, use feature_set1 for both plots
  if (same_feature_set) {
    feature_set2 <- feature_set1
  }

  # Create individual plots and save them as PDFs
  plots_list <- lapply(seq_along(seurat_objs_list), function(i) {
    if (i %% 2 == 1) { # if i is odd
      create_feature_scatter_plot(seurat_objs_list[[i]], feature_set1$feature1, feature_set1$feature2, save = TRUE, file_name = paste0(file_name_base, i), path = path)
    } else { # if i is even
      create_feature_scatter_plot(seurat_objs_list[[i]], feature_set2$feature1, feature_set2$feature2, save = TRUE, file_name = paste0(file_name_base, i), path = path)
    }
  })

  # Add plots together and save as a new PDF
  pdf(paste0(path, file_name_base, "_combined.pdf"), width = 8 * length(plots_list), height = 6)
  combined_plot <- cowplot::plot_grid(plotlist = plots_list)
  print(combined_plot)
  dev.off()

  # Delete the individual plot files
  for (i in seq_along(seurat_objs_list)) {
    file.remove(paste0(path, file_name_base, i, ".pdf"))
  }

  return(combined_plot)
}

annotate_with_clustifyR <- function(clustered_seurat_obj, output) {
  # Access the markers path
  markers_path <- config$score_and_plot_markers$known_markers_path
  markers <- read.csv2(markers_path, sep = "\t", header = TRUE)
  markers_df <- data.frame(markers$gene, markers$Cell.type)
  colnames(markers_df) <- c("gene", "cluster")

  # Convert Seurat object to expression matrix
  expr_matrix <- GetAssayData(clustered_seurat_obj, slot = "data")

  # Clustify lists with explicit matrix input
  list_res <- clustify_lists(
    input = expr_matrix,  # Use expression matrix directly
    metadata = clustered_seurat_obj@meta.data,  # Pass metadata separately
    cluster_col = "seurat_clusters",
    marker = markers_df,
    metric = "pct",
    marker_inmatrix = FALSE,
    obj_out = FALSE
  )

  if (length(unique(list_res)) > 1) {
    p1 <- plot_cor_heatmap(
      cor_mat = list_res,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      legend_title = "pct"
    )
    pdf(paste0(output, "correlation_heatmap.pdf"), width = 8, height = 6)
    print(p1)
    dev.off()
  } else {
    print(paste0("Insufficient distinct values in list_res for heatmap plotting: ", unique(list_res)))
    message("Insufficient distinct values in list_res for heatmap plotting.")
  }

  # Call cell types
  list_res2 <- cor_to_call(
    cor_mat = list_res,
    cluster_col = "seurat_clusters"
  )

  # Add clustifyr calls as metadata
  clust_call <- call_to_metadata(
    res = list_res2,
    metadata = clustered_seurat_obj@meta.data,
    cluster_col = "seurat_clusters",
    rename_prefix = "clustifyr_call"
  )
  clustered_seurat_obj <- AddMetaData(clustered_seurat_obj, metadata = clust_call)

  # Plot with clustifyr annotations
  pc <- DimPlot(clustered_seurat_obj,
    reduction = "umap", group.by = "clustifyr_call_type", label = TRUE,
    label.size = 3, repel = TRUE
  ) + ggtitle("Clustifyr annotated labels") +
    guides(fill = guide_legend(label.theme = element_text(size = 8)))
  pdf(paste0(output, "clustifyr_marker_annotation_umap.pdf"), width = 11, height = 6)
  print(pc)
  dev.off()

  # Save object with clustifyr annotation
  saveRDS(clustered_seurat_obj, file = paste0(output, "seurat_obj_clustifyr.rds"))

  return(clustered_seurat_obj)
}

clean_environment <- function(list_to_remove) {
  lapply(list_to_remove, function(x) rm(list = x, envir = .GlobalEnv))
  gc(full = TRUE)
}

create_sample_output_dir <- function(base_dir, sample_name) {
  sample_dir <- file.path(base_dir, sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(sample_dir)
}




process_consistent_species <- function(batch_corrected_obj, sample_output_dir, config, sample_name) {
  message("Processing consistent species data...")
  
  # Run UMAP
  umap_seurat_obj <- run_umap(batch_corrected_obj, sample_output_dir)
  
  # Perform clustering
  clustered_seurat_obj <- perform_clustering(umap_seurat_obj, sample_output_dir)
  
  # Save intermediate object
  saveRDS(clustered_seurat_obj, file = file.path(sample_output_dir, "clustered_seurat_obj.rds"))
  
  # Process based on DE method
  if (config$DE_method == "Seurat") {
    de_results <- find_differentially_expressed_features(clustered_seurat_obj, sample_output_dir)
    analyze_known_markers(clustered_seurat_obj, de_results, sample_output_dir)
    final_obj <- clustered_seurat_obj  # Store for return
  } else if (config$DE_method == "Scran") {
    annot_df <- score_and_plot_markers(clustered_seurat_obj, sample_output_dir)
    
    if (config$score_and_plot_markers$known_markers) {
      new_df_ordered <- annot_df[order(as.numeric(annot_df$Cluster)), ]
      new_cluster_ids <- new_df_ordered$Cell.type
      labeled_seurat_obj <- annotate_clusters_and_save(clustered_seurat_obj, new_cluster_ids, sample_output_dir)
      clustifyR_obj <- annotate_with_clustifyR(clustered_seurat_obj, sample_output_dir)
      final_obj <- list(labeled_seurat_obj = labeled_seurat_obj, clustifyR_obj = clustifyR_obj)
    } else {
      final_obj <- clustered_seurat_obj
    }
  }
  
  # Save final processed object
  saveRDS(final_obj, file = file.path(sample_output_dir, "final_processed_obj.rds"))
  
  return(final_obj)
}

perform_orthologous_gene_analysis <- function(processed_seurat_objs, config, output_dir) {
  # Add debug logging
  message("Starting orthologous gene analysis...")
  
  if (!species_are_all_same(config)) {
    # Categorize samples by species
    message("Categorizing samples by species...")
    categorized_samples <- categorize_samples_by_species(processed_seurat_objs, config)
    ref_name <- categorized_samples$ref_name
    query_name <- categorized_samples$query_name
    
    # Debug logging
    message("Reference name: ", ref_name)
    message("Query name: ", query_name)
    
    # Check if we have valid reference and query objects
    if (length(categorized_samples$ref_objects) > 0 && length(categorized_samples$query_objects) > 0) {
      ref_obj <- categorized_samples$ref_objects[[1]]
      query_obj <- categorized_samples$query_objects[[1]]
      
      # Verify objects have required data
      if (is.null(ref_obj) || is.null(query_obj)) {
        stop("Reference or query object is NULL")
      }
      
      # Check for required reductions
      if (!"harmony" %in% names(Reductions(ref_obj))) {
        message("Warning: Harmony reduction not found in reference object. Running harmony...")
        ref_obj <- RunHarmony(ref_obj, group.by.vars = "orig.ident")
      }
      
      if (!"harmony" %in% names(Reductions(query_obj))) {
        message("Warning: Harmony reduction not found in query object. Running harmony...")
        query_obj <- RunHarmony(query_obj, group.by.vars = "orig.ident")
      }
      
      # Get feature lists
      message("Getting feature lists...")
      feature_list_Q <- VariableFeatures(query_obj)
      feature_list_R <- VariableFeatures(ref_obj)
      
      # Get scaled data
      message("Getting scaled data...")
      scaled_matrix_Q <- GetAssayData(query_obj, layer = "scale.data")
      scaled_matrix_R <- GetAssayData(ref_obj, layer = "scale.data")
      
      # Get harmony embeddings
      message("Getting harmony embeddings...")
      harmony_embeddings_Q <- Embeddings(query_obj, reduction = "harmony")
      harmony_embeddings_R <- Embeddings(ref_obj, reduction = "harmony")
      
      # Get project names
      message("Getting project names...")
      ref_project <- unique(ref_obj$orig.ident)[1]
      query_project <- unique(query_obj$orig.ident)[1]
      ref_project <- sub("_lane.*", "", ref_project)
      query_project <- sub("_lane.*", "", query_project)
      project_names <- c(ref_project, query_project)
      
      # Subset orthologs
      message("Subsetting orthologs...")
      objs.list <- ortholog_subset(ref_obj, query_obj, project_names)
      
      ref.seurat <- objs.list[[1]]
      query.seurat <- objs.list[[2]]
      orthologs <- objs.list[[3]]
      
      # Save intermediate objects
      saveRDS(ref.seurat, file = file.path(output_dir, "ref_ortho-subset_seurat.rds"))
      saveRDS(query.seurat, file = file.path(output_dir, "query_ortho-subset_seurat.rds"))
      
      # Get metadata
      message("Getting metadata...")
      ref.seurat <- get_metadata(ref.seurat, "ref")
      query.seurat <- get_metadata(query.seurat, "query")
      
      # Add variable features back
      message("Adding variable features...")
      VariableFeatures(query.seurat) <- feature_list_Q
      VariableFeatures(ref.seurat) <- feature_list_R
      
      # Get the cells and features that are in the subsetted object
      query_cells <- colnames(query.seurat)
      ref_cells <- colnames(ref.seurat)
      query_features <- rownames(query.seurat)
      ref_features <- rownames(ref.seurat)

      # Debug: Print dimensions and check features/cells before subsetting
      message("Debug: Checking dimensions before subsetting")
      message("Original scaled matrix R dimensions: ", nrow(scaled_matrix_R), " x ", ncol(scaled_matrix_R))
      message("Requested features: ", length(ref_features))
      message("Requested cells: ", length(ref_cells))

      # Check which features/cells are missing
      missing_features <- ref_features[!ref_features %in% rownames(scaled_matrix_R)]
      missing_cells <- ref_cells[!ref_cells %in% colnames(scaled_matrix_R)]

      if(length(missing_features) > 0) {
        message("Missing features in scaled matrix: ", paste(head(missing_features, 5), collapse=", "), "...")
      }
      if(length(missing_cells) > 0) {
        message("Missing cells in scaled matrix: ", paste(head(missing_cells, 5), collapse=", "), "...")
      }

      # Only subset with features and cells that exist
      valid_features <- ref_features[ref_features %in% rownames(scaled_matrix_R)]
      valid_cells <- ref_cells[ref_cells %in% colnames(scaled_matrix_R)]

      # Subset using only valid features and cells
      scaled_matrix_R <- scaled_matrix_R[valid_features, valid_cells, drop = FALSE]

      # Do the same for query matrix
      message("Debug: Checking query matrix dimensions")
      message("Original scaled matrix Q dimensions: ", nrow(scaled_matrix_Q), " x ", ncol(scaled_matrix_Q))

      missing_features_Q <- query_features[!query_features %in% rownames(scaled_matrix_Q)]
      missing_cells_Q <- query_cells[!query_cells %in% colnames(scaled_matrix_Q)]

      if(length(missing_features_Q) > 0) {
        message("Missing features in query matrix: ", paste(head(missing_features_Q, 5), collapse=", "), "...")
      }
      if(length(missing_cells_Q) > 0) {
        message("Missing cells in query matrix: ", paste(head(missing_cells_Q, 5), collapse=", "), "...")
      }

      valid_features_Q <- query_features[query_features %in% rownames(scaled_matrix_Q)]
      valid_cells_Q <- query_cells[query_cells %in% colnames(scaled_matrix_Q)]

      scaled_matrix_Q <- scaled_matrix_Q[valid_features_Q, valid_cells_Q, drop = FALSE]

      # Update the Seurat objects to match the available scaled data
      ref.seurat <- subset(ref.seurat, features = valid_features, cells = valid_cells)
      query.seurat <- subset(query.seurat, features = valid_features_Q, cells = valid_cells_Q)

      # Add scaled data back
      LayerData(query.seurat, assay = "RNA", layer = "scale.data") <- scaled_matrix_Q
      LayerData(ref.seurat, assay = "RNA", layer = "scale.data") <- scaled_matrix_R
      # Subset the harmony embeddings to match the cells in the Seurat objects
      harmony_embeddings_Q <- harmony_embeddings_Q[query_cells, , drop = FALSE]
      harmony_embeddings_R <- harmony_embeddings_R[ref_cells, , drop = FALSE]
      # Create dimension reduction objects with matched embeddings
      query.seurat[["harmony"]] <- CreateDimReducObject(
        embeddings = harmony_embeddings_Q,
        key = "harmony_",
        assay = DefaultAssay(query.seurat)
      )
      ref.seurat[["harmony"]] <- CreateDimReducObject(
        embeddings = harmony_embeddings_R,
        key = "harmony_",
        assay = DefaultAssay(ref.seurat)
      )

      # Add debug messages to verify dimensions
      message("Query cells: ", ncol(query.seurat), " Query embeddings: ", nrow(harmony_embeddings_Q))
      message("Ref cells: ", ncol(ref.seurat), " Ref embeddings: ", nrow(harmony_embeddings_R))
      # Add debug messages
      message("Query scaled data dimensions: ", nrow(scaled_matrix_Q), " x ", ncol(scaled_matrix_Q))
      message("Query Seurat object dimensions: ", nrow(query.seurat), " x ", ncol(query.seurat))
      message("Ref scaled data dimensions: ", nrow(scaled_matrix_R), " x ", ncol(scaled_matrix_R))
      message("Ref Seurat object dimensions: ", nrow(ref.seurat), " x ", ncol(ref.seurat))
      
      # Return the list of objects
      obj.list2 <- list()
      obj.list2[[ref_name]] <- ref.seurat
      obj.list2[[query_name]] <- query.seurat
      
      # Save final results
      saveRDS(obj.list2, file = file.path(output_dir, "ortholog_objs_list.rds"))
      
      return(obj.list2)
    } else {
      stop("No reference and query objects found for orthologous gene analysis.")
    }
  } else {
    message("All samples are of the same species. Skipping orthologous gene analysis.")
    return(NULL)
  }
}

get_metadata <- function(seurat_obj, type) {
  # Get parameters from config
  metadata_file_ref <- config$get_metadata$metadata_file_ref
  metadata_file_query <- config$get_metadata$metadata_file_query
  metadata_subset1 <- config$get_metadata$metadata_subset1
  metadata_subset2 <- config$get_metadata$metadata_subset2
  # check type
  if (type == "ref") {
    metadata_file <- metadata_file_ref
    metadata_subset <- metadata_subset1
  } else if (type == "query") {
    metadata_file <- metadata_file_query
    metadata_subset <- metadata_subset2
  } else {
    stop("Invalid type. Please choose 'ref' or 'query'.")
  }
  # read in metadata file
  metadata <- read.csv2(metadata_file, sep="\t", header=TRUE, row.names=1)
  # subset metadata if needed
  if (metadata_subset != "NA") {
    metadata <- subset(metadata, metadata$source==metadata_subset)
  }
  # get complete metadata
  metadata<- as.data.frame(metadata)
  complete.cases(metadata)
  metadata<-na.omit(metadata)
  # add metadata to seurat object
  seurat_obj <- AddMetaData(object=seurat_obj, metadata=metadata)

  return(seurat_obj)
}

process_orthologous_objects <- function(seurat_obj, output_dir, config, sample_name){
  sample_output_dir <- file.path(output_dir, sample_name)
  sample_output_dir <- paste0(sample_output_dir, "/")
  seurat_obj <- run_and_visualize_pca(seurat_obj, sample_output_dir)
  umap_seurat_obj <- run_umap(seurat_obj, sample_output_dir)
  clustered_seurat_obj <- perform_clustering(umap_seurat_obj, sample_output_dir)
  if (config$DE_method == "Seurat") {
    de_results <- find_differentially_expressed_features(clustered_seurat_obj, sample_output_dir)
    analyze_known_markers(clustered_seurat_obj, de_results, sample_output_dir)
  } else if (config$DE_method == "Scran") {
    annot_df <- score_and_plot_markers(clustered_seurat_obj, sample_output_dir)
  }
  # Optionally annotate clusters and save if known markers are to be scored and plotted
  if (config$score_and_plot_markers$known_markers) {
    new_df_ordered <- annot_df[order(as.numeric(annot_df$Cluster)), ]
    # Get new cluster names ordered by cluster number
    new_cluster_ids <- new_df_ordered$Cell.type
    labeled_seurat_obj <- annotate_clusters_and_save(clustered_seurat_obj, new_cluster_ids, sample_output_dir)

    clustifyR_obj <- annotate_with_clustifyR(clustered_seurat_obj, sample_output_dir)
    # Append both annotated objects to the return list
    return(list(labeled_seurat_obj = labeled_seurat_obj, clustifyR_obj = clustifyR_obj))
  } else {
    # If no known markers scoring and plotting, return the clustered object only
    return(list(clustered_seurat_obj = clustered_seurat_obj))
  }
}

ortholog_subset <- function(ref_seurat, query_seurat, project_names, path = output) {
  # Get parameters from config
  ortholog_file <- config$ortholog_subset$ortholog_file
  # read in orrtholog file
  orthologs <- read.csv2(ortholog_file, sep = "\t", header = TRUE, row.names = 1)
  # ortholog gene list
  orthos <- as.vector(orthologs$human.gene.name)
  # turn seurat objects into matrices
  ref.matrix <- GetAssayData(ref_seurat, slot = "data")
  query.matrix <- GetAssayData(query_seurat, slot = "data")

  # subset matrix so they match only orthologous genes
  ref.matrix.flt <- ref.matrix[rownames(ref.matrix) %in% orthos, ]
  # subset orthologs to they only match expr matrix
  orthologs <- orthologs[orthologs$human.gene.name %in% rownames(ref.matrix.flt), ]
  # switch reference genes to query genes
  # reorder based on exp matrix genes
  ind_reorder <- match(rownames(ref.matrix.flt), orthologs$human.gene.name)
  orthologs.reorder <- orthologs[ind_reorder, ]
  # now replace rownames of exp matrix with pig genes
  pig.orthos <- as.vector(orthologs.reorder$pig.gene.name)
  row.names(ref.matrix.flt) <- pig.orthos
  # subset pig data so that there is same genes
  query.matrix.flt <- query.matrix[rownames(query.matrix) %in% pig.orthos, ]
  dim(query.matrix.flt)
  # subset human refernce based on gamm data
  ref.matrix.flt <- ref.matrix.flt[rownames(ref.matrix.flt) %in% rownames(query.matrix.flt), ]
  dim(ref.matrix.flt)
  # get rid of duplicates in reference by averaging them
  # Aggregate rows based on row names
  ref.matrix.flt <- aggregate(ref.matrix.flt, by = list(rownames(ref.matrix.flt)), FUN = mean)
  # use column 1 as rownames
  rownames(ref.matrix.flt) <- ref.matrix.flt$Group.1
  ref.matrix.flt <- select(ref.matrix.flt, -c("Group.1"))
  # turn back to sparse matrix
  library(Matrix)
  ref.matrix.flt <- as.matrix(ref.matrix.flt)
  ref.matrix.flt <- Matrix(ref.matrix.flt, sparse = TRUE)
  # create seurat object
  ref.seurat <- CreateSeuratObject(ref.matrix.flt, project = project_names[[1]], assay = "RNA")
  query.seurat <- CreateSeuratObject(query.matrix.flt, project = project_names[[2]], assay = "RNA")
  # put in a list to export
  obj.list <- list(ref.seurat, query.seurat, orthologs)
  return(obj.list)
}

categorize_samples_by_species <- function(sample_specific_list, config) {
  # Initialize empty lists to hold ref and query objects
  ref_objects <- list()
  query_objects <- list()

  # Iterate over each sample in the list
  for (i in seq_along(sample_specific_list)) {
    # Assuming each element in sample_specific_list has a unique identifier or name
    # This part may need adjustment depending on the structure of your sample_specific_list
    sample_name <- names(sample_specific_list)[i]
    species <- get_species_by_sample_name(sample_name, config)

    # Check if species is "Human" and categorize the sample accordingly
    if (tolower(species) == "human") {
      # Append to ref_objects
      ref_objects[[sample_name]] <- sample_specific_list[[i]]
      ref_name <- sample_name
    } else {
      # Append to query_objects
      query_objects[[sample_name]] <- sample_specific_list[[i]]
      query_name <- sample_name
    }
  }

  # Return a list containing both lists of objects
  return(list(ref_objects = ref_objects, query_objects = query_objects, ref_name = ref_name, query_name = query_name))
}

species_are_all_same <- function(config) {
  species_values <- future_map_chr(config$fastq_alignment, "species")
  all(species_values == species_values[1])
}

get_species_by_sample_name <- function(sample_name, config) {
  for (sample in names(config$fastq_alignment)) {
    sample_details <- config$fastq_alignment[[sample]]
    if (sample_details$NAME == sample_name) {
      return(sample_details$species)
    }
  }
  return(NA) # Return NA if no species is found for the sample name
}

run_leiden_clustering <- function(seurat_obj, resolution) {
  # Import Python's leidenalg
  leidenalg <- reticulate::import("leidenalg")
  igraph <- reticulate::import("igraph")
  
  # Get the SNN graph from Seurat object
  snn_graph <- seurat_obj@graphs$RNA_snn
  
  # Convert to igraph format
  edges <- which(snn_graph != 0, arr.ind = TRUE)
  weights <- snn_graph[edges]
  
  # Create Python igraph object
  g <- igraph$Graph(edges = edges - 1,  # Python uses 0-based indexing
                   directed = FALSE,
                   weights = weights)
  
  # Run Leiden clustering
  partition <- leidenalg$find_partition(g,
                                      leidenalg$RBConfigurationVertexPartition,
                                      resolution_parameter = resolution)
  
  # Convert results back to R
  clusters <- as.factor(partition$membership + 1)  # Convert back to 1-based indexing
  names(clusters) <- colnames(seurat_obj)
  
  # Add clusters to Seurat object
  seurat_obj$leiden_clusters <- clusters
  
  return(seurat_obj)
}