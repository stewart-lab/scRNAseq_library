read_aligned_data <- function(base_directory, project_name, output_base_dir) {
  # Set up directories in parallel
  dirs <- list(
    filtered = list(
      source = paste0(base_directory, "/filtered/"),
      target = paste0(output_base_dir, "/alignment_", project_name, "/filtered/")
    ),
    raw = list(
      source = paste0(base_directory, "/raw/"),
      target = paste0(output_base_dir, "/alignment_", project_name, "/raw/")
    )
  )
  
  # Create directories in parallel
  future_map(dirs, function(d) {
    dir.create(d$target, recursive = TRUE, showWarnings = FALSE)
  })
  
  # Copy files in parallel
  future_map(dirs, function(d) {
    gz_files <- list.files(d$source, pattern = "\\.gz$", full.names = TRUE)
    future_map(gz_files, function(file) {
      file.copy(file, file.path(d$target, basename(file)), overwrite = TRUE)
    })
  })
  
  # Copy summary file
  tryCatch({
    file.copy(
      paste0(base_directory, "/Summary.csv"),
      paste0(output_base_dir, "/alignment_", project_name, "/Alignment_Summary.csv"),
      overwrite = TRUE
    )
  }, error = function(e) message("Summary.csv copy failed: ", e$message))
  
  # Read data in parallel
  data <- future_map(dirs, function(d) Read10X(d$source))
  
  list(
    filtered = data$filtered,
    raw = data$raw,
    project = project_name
  )
}

prep_seurat_and_soupX <- function(data.raw, data, project) {
  dims_umap <- 1:config$prep_seurat_and_soupX$dims
  umap.method <- config$prep_seurat_and_soupX$umap.method
  tfidfMin <- config$prep_seurat_and_soupX$tfidfMin
  min.cells <- config$prep_seurat_and_soupX$min.cells
  
  # Create objects in parallel
  objects <- future_map(list(
    sc = list(raw = data.raw, filtered = data),
    seurat = list(data = data, project = project, min.cells = min.cells)
  ), function(x) {
    if ("raw" %in% names(x)) {
      SoupChannel(x$raw, x$filtered)
    } else {
      CreateSeuratObject(counts = x$data, project = x$project, min.cells = x$min.cells)
    }
  })
  
  # Clean up memory
  rm(data.raw, data)
  gc(full = TRUE)
  
  # Process Seurat object sequentially but with parallel internals
  seurat_obj <- objects$seurat
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- debug_run_umap(seurat_obj, path)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_umap, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, verbose = TRUE)
  
  # Process SoupX in parallel
  meta <- seurat_obj@meta.data
  umap <- seurat_obj@reductions$umap@cell.embeddings
  
  sc <- objects$sc %>%
    setClusters(setNames(meta$seurat_clusters, rownames(meta))) %>%
    autoEstCont(tfidfMin = tfidfMin, forceAccept = TRUE)
  
  out <- adjustCounts(sc, roundToInt = TRUE)
  
  list(seurat_obj = seurat_obj, meta = meta, umap = umap, out = out)
}

process_lane <- function(lane) {
  options(future.globals.maxSize = 131072 * 1024^2)
  
  # Process in parallel
  aligned_data <- read_aligned_data(lane$base_directory, lane$name, output_base_dir)
  
  soupX_obj <- prep_seurat_and_soupX(
    data.raw = aligned_data$raw,
    data = aligned_data$filtered,
    project = aligned_data$project
  )
  
  rm(aligned_data)
  gc(full = TRUE)
  
  # Create Seurat and SCE objects
  feature_set1 <- list(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  sce_obj <- create_seurat_and_sce(
    out = soupX_obj$out,
    project = lane$name,
    feature_set = feature_set1,
    sample_name = lane$name  # Now the function accepts this parameter
  )
  
  rm(soupX_obj)
  gc(full = TRUE)
  
  return(sce_obj)
}

# Helper function to clean environment
clean_environment <- function(list_to_remove) {
  rm(list = list_to_remove)
  gc(full = TRUE)
}

create_seurat_and_sce <- function(out, project, feature_set, sample_name = NULL) {
  # Use project as sample_name if not provided
  sample_name <- if (is.null(sample_name)) project else sample_name
  
  # Create singular Seurat object
  seu <- CreateSeuratObject(counts = out, project = project, data = NULL)
  
  # Add sample ID to metadata
  seu$sample_id <- sample_name
  
  # Calculate basic metrics
  seu[["nCount_RNA"]] <- Seurat::GetAssayData(seu, slot = "counts") %>% colSums()
  seu[["nFeature_RNA"]] <- Seurat::GetAssayData(seu, slot = "counts") %>% 
    apply(2, function(x) sum(x > 0))
  
  # Process features in parallel if needed
  if (length(feature_set) > 0) {
    future_map(names(feature_set), function(fname) {
      feature <- feature_set[[fname]]
      if (!feature %in% names(seu@meta.data)) {
        if (feature == "percent.mt") {
          seu[[feature]] <<- PercentageFeatureSet(seu, pattern = "^MT-", assay = "RNA")
        }
      }
    })
  }
  
  # Convert to SCE after all features are calculated
  sce <- as.SingleCellExperiment(seu)
  
  # Add metadata to SCE
  sce$sample_id <- sample_name
  
  # Add all features from Seurat to SCE
  for (feature in names(seu@meta.data)) {
    sce[[feature]] <- seu[[feature]]
  }
  
  # Clean memory
  gc(full = TRUE)
  
  return(list(seu = seu, sce = sce))
}