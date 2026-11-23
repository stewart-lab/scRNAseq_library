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
  tryCatch(
    {
      file.copy(
        paste0(base_directory, "/Summary.csv"),
        paste0(output_base_dir, "/alignment_", project_name, "/Alignment_Summary.csv"),
        overwrite = TRUE
      )
    },
    error = function(e) message("Summary.csv copy failed: ", e$message)
  )

  # Read data in parallel
  data <- future_map(dirs, function(d) Read10X(d$source))

  list(
    filtered = data$filtered,
    raw = data$raw,
    project = project_name
  )
}

read_aligned_data2 <- function(base_directory, project_name, output_base_dir) {
  print("reading in data")
  raw <- list(
    source = base_directory,
    target = paste0(output_base_dir, "/alignment_", project_name, "/raw/")
  )
  # create
  dir.create(raw$target, recursive = TRUE, showWarnings = FALSE)

  gz_files <- list.files(raw$source, pattern = "\\.gz$", full.names = TRUE)
  future_map(gz_files, function(file) {
    file.copy(file, file.path(raw$target, basename(file)), overwrite = TRUE)
  })
  # Read data (no need for future_map here)
  print(file.path(raw$source))
  data_raw <- Read10X(file.path(raw$source))
  print(data_raw)
  # Return in the same structure as read_aligned_data for compatibility
  list(
    filtered = NULL, # or data_raw if you want to use raw as filtered
    raw = data_raw,
    project = project_name
  )
}

# Parse Biosciences (split-pipe) DGE output: count_matrix.mtx is cells (rows) x
# genes (cols) - the transpose of the genes x cells layout Seurat/SoupX expect -
# with barcodes and gene names stored separately in cell_metadata.csv/all_genes.csv.
read_dge_matrix <- function(dge_dir, target_dir) {
  # Stage a copy of the input files alongside the rest of the run's output
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dge_files <- list.files(dge_dir, full.names = TRUE)
  file.copy(dge_files, file.path(target_dir, basename(dge_files)), overwrite = TRUE)

  mat <- Matrix::readMM(file.path(dge_dir, "count_matrix.mtx"))
  mat <- Matrix::t(mat)
  mat <- methods::as(mat, "CsparseMatrix")

  genes <- read.csv(file.path(dge_dir, "all_genes.csv"), stringsAsFactors = FALSE)
  cells <- read.csv(file.path(dge_dir, "cell_metadata.csv"), stringsAsFactors = FALSE)

  rownames(mat) <- make.unique(genes$gene_name)
  colnames(mat) <- cells$bc_wells

  mat
}

read_aligned_data3 <- function(base_directory, project_name, output_base_dir) {
  list(
    filtered = read_dge_matrix(
      file.path(base_directory, "DGE_filtered"),
      paste0(output_base_dir, "/alignment_", project_name, "/filtered/")
    ),
    raw = read_dge_matrix(
      file.path(base_directory, "DGE_unfiltered"),
      paste0(output_base_dir, "/alignment_", project_name, "/raw/")
    ),
    project = project_name
  )
}


# ---------------------------------------------------------------------------
# Species subsetting for human+pig chimeric libraries
# (see README_chimeric_analysis.md).  Optional per-sample keys in
# config$fastq_alignment$<SAMPLE>:
#   SPECIES_BARCODES  barcodes_<species>.txt from 03_classify_species.py:
#                     a path, or {"lane1": path, "lane2": path}
#   SOUPX_RHO         SoupX contamination fraction (rho_pig / rho_human from 03):
#                     a number, or {"lane1": 0.17, "lane2": 0.17}
# Lanes are matched by the suffix after "<NAME>_" (e.g. NAME_lane1 -> "lane1").
# ---------------------------------------------------------------------------
species_subset_for_lane <- function(lane_name, config) {
  samples <- Filter(function(s) is.list(s) && !is.null(s$NAME), config$fastq_alignment)
  hit <- NULL
  for (s in samples) {
    if (lane_name == s$NAME || startsWith(lane_name, paste0(s$NAME, "_"))) {
      if (is.null(hit) || nchar(s$NAME) > nchar(hit$NAME)) hit <- s
    }
  }
  if (is.null(hit) || (is.null(hit$SPECIES_BARCODES) && is.null(hit$SOUPX_RHO))) {
    return(NULL)
  }
  lane_key <- sub("^_", "", substring(lane_name, nchar(hit$NAME) + 1))
  pick <- function(x, what) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) {
      if (!nzchar(lane_key)) {
        stop(sprintf(paste0("%s for sample %s is a per-lane map, but this sample has no ",
                            "lane suffix (COMBINE_LANES / merged lanes): give a single value"),
                     what, hit$NAME))
      }
      v <- x[[lane_key]]
      if (is.null(v)) {
        stop(sprintf("%s for sample %s has no entry for lane '%s' (has: %s)",
                     what, hit$NAME, lane_key, paste(names(x), collapse = ", ")))
      }
      return(v)
    }
    x
  }
  ss <- list(
    sample = hit$NAME, lane_key = lane_key, species = hit$species,
    barcodes = pick(hit$SPECIES_BARCODES, "SPECIES_BARCODES"),
    rho = pick(hit$SOUPX_RHO, "SOUPX_RHO")
  )
  if (!is.null(ss$rho)) {
    ss$rho <- as.numeric(ss$rho)
    if (is.na(ss$rho) || ss$rho <= 0 || ss$rho >= 1) {
      stop(sprintf("SOUPX_RHO for %s lane %s must be between 0 and 1", hit$NAME, lane_key))
    }
  }
  if (!is.null(ss$barcodes) && !(ss$species %in% c("pig", "human"))) {
    stop(sprintf(paste0("Sample %s uses SPECIES_BARCODES, so its 'species' must be ",
                        "\"pig\" or \"human\" (the species kept), not \"%s\""),
                 hit$NAME, ss$species))
  }
  ss
}

# STARsolo/03 barcodes are bare sequences; drop any "-1"-style suffix so the
# lists match Read10X column names either way.
norm_barcodes <- function(x) sub("-[0-9]+$", "", trimws(x))

read_species_barcodes <- function(path) {
  if (!file.exists(path)) stop("SPECIES_BARCODES file not found: ", path)
  bc <- norm_barcodes(readLines(path))
  bc <- unique(bc[nzchar(bc)])
  if (length(bc) == 0) stop("SPECIES_BARCODES file is empty: ", path)
  bc
}

prep_seurat_and_soupX <- function(data.raw, data, project, rho = NULL, keep_cells = NULL) {
  dims_umap <- 1:config$prep_seurat_and_soupX$dims
  umap.method <- config$prep_seurat_and_soupX$umap.method
  tfidfMin <- config$prep_seurat_and_soupX$tfidfMin
  min.cells <- config$prep_seurat_and_soupX$min.cells

  # Species subsetting: SoupX runs on the whole library (soup from all empty
  # droplets, cells = STARsolo-filtered cells PLUS any species barcodes
  # EmptyDrops missed); the corrected matrix is subset to keep_cells at the end.
  subset_log <- NULL
  if (!is.null(keep_cells)) {
    raw_bc <- norm_barcodes(colnames(data.raw))
    filt_bc <- norm_barcodes(colnames(data))
    in_raw <- keep_cells %in% raw_bc
    if (!any(in_raw)) {
      stop("None of the SPECIES_BARCODES are in the raw matrix for ", project,
           " -- wrong lane or barcode file?")
    }
    toc_bc <- union(filt_bc, keep_cells[in_raw])
    data <- data.raw[, match(toc_bc, raw_bc)]
    subset_log <- data.frame(
      lane = project, rho = if (is.null(rho)) NA else rho,
      species_barcodes = length(keep_cells),
      species_barcodes_not_in_raw = sum(!in_raw),
      starsolo_filtered_cells = length(filt_bc),
      species_barcodes_not_in_filtered = sum(in_raw & !(keep_cells %in% filt_bc)),
      soupx_cells = length(toc_bc)
    )
    message(sprintf(
      "[%s] species subset: %d barcodes (%d not in raw); SoupX on %d cells (%d STARsolo-filtered + %d added)",
      project, length(keep_cells), sum(!in_raw), length(toc_bc), length(filt_bc),
      length(toc_bc) - length(filt_bc)))
  }

  # Create objects in parallel
  objects <- future_map(list(
    sc = list(raw = data.raw, filtered = data),
    seurat = list(data = data, project = project, min.cells = min.cells)
  ), function(x) {
    if ("raw" %in% names(x)) {
      SoupChannel(x$raw, x$filtered)
    } else {
      # names.field = 0 disables Seurat's default "_"-delimited barcode parsing for
      # orig.ident, which otherwise clobbers `project` for underscore-delimited barcodes
      # (e.g. Parse Biosciences bc_wells like "26_01_01").
      CreateSeuratObject(counts = x$data, project = x$project, min.cells = x$min.cells, names.field = 0)
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
    setClusters(setNames(meta$seurat_clusters, rownames(meta)))
  if (is.null(rho)) {
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE)
  } else {
    # rho measured from the other species (03_classify_species.py) replaces
    # SoupX's estimate; adjustCounts still removes soup gene by gene
    message(sprintf("[%s] SoupX: using configured rho = %.3f (autoEstCont skipped)", project, rho))
    sc <- setContaminationFraction(sc, rho, forceAccept = TRUE)
  }

  out <- adjustCounts(sc, roundToInt = TRUE)

  if (!is.null(keep_cells)) {
    keep <- norm_barcodes(colnames(out)) %in% keep_cells
    out <- out[, keep]
    meta <- meta[keep, , drop = FALSE]
    umap <- umap[keep, , drop = FALSE]
    subset_log$cells_kept <- ncol(out)
    subset_log$median_umi_before_soupx <- median(Matrix::colSums(sc$toc[, keep]))
    subset_log$median_umi_after_soupx <- median(Matrix::colSums(out))
    write.table(subset_log, file.path(output_base_dir, paste0("species_subset_", project, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[%s] kept %d species cells after SoupX", project, ncol(out)))
  }

  list(seurat_obj = seurat_obj, meta = meta, umap = umap, out = out)
}

process_lane <- function(lane, parse) {
  options(future.globals.maxSize = 131072 * 1024^2)
  # make sure directory is formatted right
  lane$base_directory <- normalizePath(lane$base_directory, mustWork = FALSE)
  # Process in parallel
  ss <- species_subset_for_lane(lane$name, config)
  if (!is.null(ss) && parse == TRUE) {
    stop("SPECIES_BARCODES / SOUPX_RHO are only supported for 10x (STARsolo) data")
  }
  if (parse == TRUE) {
    aligned_data <- read_aligned_data3(lane$base_directory, lane$name, output_base_dir)
  } else {
    aligned_data <- read_aligned_data(lane$base_directory, lane$name, output_base_dir)
  }

  soupX_obj <- prep_seurat_and_soupX(
    data.raw = aligned_data$raw,
    data = aligned_data$filtered,
    project = aligned_data$project,
    rho = ss$rho,
    keep_cells = if (!is.null(ss$barcodes)) read_species_barcodes(ss$barcodes) else NULL
  )

  rm(aligned_data)
  gc(full = TRUE)

  # Create Seurat and SCE objects
  feature_set1 <- list(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  sce_obj <- create_seurat_and_sce(
    out = soupX_obj$out,
    project = lane$name,
    feature_set = feature_set1,
    sample_name = lane$name # Now the function accepts this parameter
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
  # names.field = 0: see note in prep_seurat_and_soupX() above.
  seu <- CreateSeuratObject(counts = out, project = project, data = NULL, names.field = 0)

  # Add sample ID to metadata
  seu$sample_id <- sample_name

  # Calculate basic metrics
  seu[["nCount_RNA"]] <- seu[["RNA"]]$counts %>% colSums()
  seu[["nFeature_RNA"]] <- seu[["RNA"]]$counts %>%
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

filter_empty_droplets <- function(data.raw) {
  # filter out empty droplets
  sce <- SingleCellExperiment(list(counts = data.raw))
  tryCatch(
    {
      e.out <- emptyDrops(counts(sce))
      data.filtered <- sce[, which(e.out$FDR <= 0.001)]
      return(data.filtered)
    },
    error = function(e) {
      print(e)
      data.filtered <- NULL
      return(data.filtered)
    }
  )
}

process_lane2 <- function(lane) {
  options(future.globals.maxSize = 131072 * 1024^2)
  if (!is.null(species_subset_for_lane(lane$name, config))) {
    stop("SPECIES_BARCODES / SOUPX_RHO need STARsolo raw + filtered matrices; ",
         "set process_matrix to false")
  }
  # make sure directory is formatted right
  lane$base_directory <- normalizePath(lane$base_directory, mustWork = FALSE)
  # Process in parallel
  aligned_data <- read_aligned_data2(lane$base_directory, lane$name, output_base_dir)
  data.filtered <- filter_empty_droplets(aligned_data$raw)

  if (is.null(data.filtered) == TRUE) {
    # Create Seurat and SCE objects
    feature_set1 <- list(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    sce_obj <- create_seurat_and_sce(
      out = aligned_data$raw,
      project = lane$name,
      feature_set = feature_set1,
      sample_name = lane$name
    )
  } else {
    soupX_obj <- prep_seurat_and_soupX(
      data.raw = aligned_data$raw,
      data = data.filtered,
      project = aligned_data$project
    )
    # Create Seurat and SCE objects
    feature_set1 <- list(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    sce_obj <- create_seurat_and_sce(
      out = soupX_obj$out,
      project = lane$name,
      feature_set = feature_set1,
      sample_name = lane$name # Now the function accepts this parameter
    )
  }

  rm(aligned_data)
  rm(soupX_obj)
  gc(full = TRUE)

  return(sce_obj) #
}
