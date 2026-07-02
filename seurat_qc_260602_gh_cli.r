#!/usr/bin/env Rscript
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# https://github.com/BaderLab/MALAT1_threshold
# https://bioconductor.org/books/3.15/OSCA.basic/quality-control.html#quality-control-outlier

# Set libraries (seurat workflow)
suppressPackageStartupMessages({
  needed <- c(
    "future",
    "data.table",
    "Matrix",
    "Seurat", 
    "tools",
    "patchwork",
    "dplyr",
    "tidyr",
    "ggplot2", # version 3.5.1
    "DoubletFinder",
    "SoupX"
  )
  to_install <- needed[!needed %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) {
    suppressMessages(install.packages(to_install, repos = "https://cloud.r-project.org"))
  }
  invisible(lapply(needed, function(pkg) suppressMessages(require(pkg, character.only = TRUE))))
})

# increase ram limit 
options(future.globals.maxSize = 80 * 1024^3) # 1024^3 = 1Gb

# args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript seurat_droplet_12-3-2026_gh_cli.R <[cellranger|starsolo]> <indir> <outdir>")
}

# ---- SETTINGS ----
aligner <- args[1]       # cellranger or starsolo    
FindNeighbors.dims <- 1:15    # Check elbow plot
FindClusters.res <- 0.4       # Turn up to find more clusters, down to find fewer clusters

# droplet defaults
MAD_devs <- 2.5              # number of deviations (captures ~99% if normally distributed)
percent.mt.max <- 20

# ---- Specify paths ----
path <- here::here()
projdir <- normalizePath(file.path(path, "../.."))
indir <- args[2]
outdir <- args[3]
if (!dir.exists(indir)) {stop("ERROR: Provided indir does not exist!\n")}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

metadir <- file.path(dirname(indir), "SraRunTable_used.csv")

# ---- Helpers ----
source(file.path(path,"utils/seurat_pipe_20-2-2026_gh.r"))
source(file.path(path,"utils/seurat_doubletfinder_20-2-2026_gh.r"))
source(file.path(path,"utils/seurat_soupx_22-2-2026_gh.r"))
source(file.path(path,"utils/seurat_integrate_24-2-2026_gh.r"))
source(file.path(path,"utils/malat1_function.R"))

# ---- Read in data (organise by experiment) ----
metadata <- fread(metadir, select = c("Run", "Experiment"))

# filtered counts
crdir_filt <- lapply(unique(metadata$Experiment), function(exp) {
  runs <- metadata[Experiment == exp, Run]
  
  # Find directories matching aligner
  if (aligner == "cellranger") {
    regex <- paste0("(", paste(runs, collapse = "|"), ")_counts")
    crdir <- list.files(indir, full.names = TRUE, pattern = regex)
    crdir_filt <- file.path(crdir, "outs", "filtered_feature_bc_matrix")
    
  } else if (aligner == "starsolo") {
    regex <- paste0("(", paste(runs, collapse = "|"), ")_Solo.out")
    crdir <- list.files(indir, full.names = TRUE, pattern = regex)
    crdir_filt <- file.path(crdir, "Gene", "filtered")
    
  } else {
    stop("ERROR: aligner must be cellranger or starsolo\n")
  }
  
  # If nothing matched at all
  if (length(crdir_filt) == 0) {
    cat("   WARNING: no directories found for experiment", exp, "- removing\n")
    return(NULL)
  }
  
  # Keep only existing directories
  exists <- dir.exists(crdir_filt)
  
  if (!all(exists)) {
    missing <- crdir_filt[!exists]
    cat("   WARNING: missing directories for experiment", exp, ":\n")
    cat("            ", paste(missing, collapse = "\n            "), "\n")
  }
  
  crdir_filt <- crdir_filt[exists]
  
  # If all directories missing → drop experiment
  if (length(crdir_filt) == 0) {
    cat("   WARNING: all directories missing for experiment", exp, "- removing\n")
    return(NULL)
  }
  
  return(crdir_filt)
})

# Remove NULL entries (experiments with no valid runs)
crdir_filt <- Filter(Negate(is.null), crdir_filt)

# raw counts
crdir_raw <- lapply(unique(metadata$Experiment), function(exp) {
  runs <- metadata[Experiment == exp, Run]
  if (aligner == "cellranger") {
    regex <- paste0("(", paste(runs, collapse = "|"), ")_counts")
    crdir <- list.files(indir, full.names = TRUE, pattern = regex)
    crdir_raw <- file.path(crdir, "outs", "raw_feature_bc_matrix")
  } else if (aligner == "starsolo") {
    regex <- paste0("(", paste(runs, collapse = "|"), ")_Solo.out")
    crdir <- list.files(indir, full.names = TRUE, pattern = regex)
    crdir_raw <- file.path(crdir, "Gene", "raw")
  } else {
    stop("ERROR: aligner not set to cellranger or starsolo, cannot find raw counts directory\n")
  }
  if (length(crdir_raw) == 0 || any(!dir.exists(crdir_raw))) {
    crdir_raw <- character(0)
  }
  return(crdir_raw)
})

# label counts with experiment
names(crdir_filt) <- unique(metadata$Experiment)
names(crdir_raw) <- unique(metadata$Experiment)

# Create name for global labelling 
nobj <- basename(outdir)

# Initialize list to store Seurat objects
sobj_list <- list()
metrics_rows <- list()

# ---- MALAT1, SoupX, MAD-QC, DoubletFinder ----
for (exp in names(crdir_filt)) {
  exp_metrics <- list(experiment = exp)
  
  # Create Seurat object name
  cat(">>> Processing", exp, "\n")
  
  qcdir <- file.path(outdir, paste0(exp, "_QC_out"))
  dir.create(qcdir, showWarnings = FALSE, recursive = TRUE)
  qcdir <- normalizePath(qcdir)
  
  # Create pdf for plots (add time to prevent rewrite crash)
  pdf_path <- file.path(qcdir, paste(nobj, exp, format(Sys.Date(), "%H-%M-%S"), "plots.pdf", sep = "_"))
  pdf(pdf_path, width = 10, height = 6)
  
  cat("  - Reading in data\n")
  sobj.filt <- Read10X(data.dir = crdir_filt[[exp]])
  
  # Create seurat object for malat1 and soupx, with no filtering
  cat("  - Creating Seurat object\n")
  sobj <- CreateSeuratObject(
    counts = sobj.filt,
    project = nobj
  )
  
  # use log-normalisation for malat1 thresholding
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(sobj, verbose = FALSE)
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = FindNeighbors.dims, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = FindClusters.res, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = FindNeighbors.dims, verbose = FALSE)

  # ---- MALAT1 thresholding ----
  # apply malat1 thresholding per experiment: 
  cat("  - Applying MALAT1 thresholding\n")
  norm_counts <- GetAssayData(sobj, assay = "RNA", layer = "data")["MALAT1",]
  threshold <- define_malat1_threshold_ggplot2(norm_counts)
  malat1_threshold <- norm_counts > threshold
  sobj$malat1_threshold <- malat1_threshold
  sobj$malat1_threshold <- factor(sobj$malat1_threshold, levels = c(TRUE, FALSE))
  print(DimPlot(sobj, reduction = "umap", group.by = "malat1_threshold"))
  good_cells <- colnames(sobj)[malat1_threshold]

  # report number of cells pre/post malat1 thresholding
  ncells_pre_malat1 <- ncol(sobj)
  sobj <- subset(sobj, cells = good_cells)
  ncells_post_malat1 <- ncol(sobj)
  percent_cells_malat1 <- (ncells_post_malat1/ncells_pre_malat1) * 100
  cat(sprintf("  - Number of cells after MALAT1 thresholding: %i (%.2f%% remaining)\n", ncells_post_malat1, percent_cells_malat1))
  
  exp_metrics$malat1_ncells_pre <- ncells_pre_malat1
  exp_metrics$malat1_ncells_post <- ncells_post_malat1
  exp_metrics$malat1_percent_cells_remaining <- percent_cells_malat1

  sobj.filt.malat1 <- GetAssayData(sobj, assay = "RNA", layer = "counts")

  # reset sobj object for SoupX
  sobj[["RNA"]]$scale.data <- NULL
  sobj <- FindVariableFeatures(sobj, verbose = FALSE)
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = FindNeighbors.dims, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = FindClusters.res, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = FindNeighbors.dims, verbose = FALSE)

  # check if raw exists, if so -> SoupX
  if (length(crdir_raw[[exp]]) > 0) {
    
    sobj.raw <- Read10X(data.dir = crdir_raw[[exp]])
    
    # ---- Run soupx helper function ----
    # sobj has to have umap
    sobj <- seurat_soupx_23.2.2026_gh(
      sobj = sobj, 
      nobj = nobj,
      sobj.raw = sobj.raw, 
      sobj.filt = sobj.filt.malat1,
      outdir = qcdir
    )
    
  } else {
    cat("   WARNING: no raw counts directory found for", exp, "continuing with filtered counts only\n")
    
    # Create seurat object with filtered counts only
    cat("  - Creating Seurat object\n")
    sobj <- CreateSeuratObject(
      counts = sobj.filt.malat1,
      min.cells = 3,
      min.features = 200,
      project = nobj
    )
  }

  # ---- Raw QC ----
  cat("  - Performing QC\n")
  
  # Refilter and add mito to new sobj
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  
  # Visualize QC metrics with violins and scatters
  print(VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts"))
  print(FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  
  # Calculate nFeature_RNA thresholds using logs and Mean Absolute Deviations (MADs)
  logFeature_RNA <- log1p(sobj$nFeature_RNA)
  logFeature_RNA.min <- median(logFeature_RNA) - MAD_devs * mad(logFeature_RNA)
  logFeature_RNA.max <- median(logFeature_RNA) + MAD_devs * mad(logFeature_RNA)
  nFeature_RNA.min <- expm1(logFeature_RNA.min)
  nFeature_RNA.max <- expm1(logFeature_RNA.max)
  
  # Calculate nCount_RNA thresholds using logs and Mean Absolute Deviations (MADs)
  logCount_RNA <- log1p(sobj$nCount_RNA)
  logCount_RNA.min <- median(logCount_RNA) - MAD_devs * mad(logCount_RNA)
  logCount_RNA.max <- median(logCount_RNA) + MAD_devs * mad(logCount_RNA)
  nCount_RNA.min <- expm1(logCount_RNA.min)
  nCount_RNA.max <- expm1(logCount_RNA.max)
  
  # plot calculated thresholds on violins
  p.feat <- VlnPlot(sobj, features = "nFeature_RNA", layer = "counts") +
    geom_hline(yintercept = nFeature_RNA.min, linetype = "dashed", color = "darkblue") +
    geom_hline(yintercept = nFeature_RNA.max, linetype = "dashed", color = "tomato") +
    ggtitle("nFeature_RNA") +
    NoLegend()
  
  p.count <- VlnPlot(sobj, features = "nCount_RNA", layer = "counts") +
    geom_hline(yintercept = nCount_RNA.min, linetype = "dashed", color = "darkblue") +
    geom_hline(yintercept = nCount_RNA.max, linetype = "dashed", color = "tomato") +
    ggtitle("nCount_RNA") +
    NoLegend()
  
  p.mito <- VlnPlot(sobj, features = "percent.mt", layer = "counts") +
    geom_hline(yintercept = percent.mt.max, linetype = "dashed", color = "tomato") +
    ggtitle("percent.mt") +
    NoLegend()
  
  print(p.feat | p.count | p.mito)
  
  # Remove cells that fail QC
  num_cells_preQC <- ncol(sobj)
  sobj <- subset(sobj, subset = 
      nFeature_RNA > nFeature_RNA.min & 
      nFeature_RNA < nFeature_RNA.max &
      nCount_RNA > nCount_RNA.min & 
      nCount_RNA < nCount_RNA.max &
      percent.mt < percent.mt.max
  )

  num_cells_postQC <- ncol(sobj)
  percent_cells_kept <- (num_cells_postQC / num_cells_preQC) * 100
  cat(sprintf("  - Number of cells after QC: %i (%.2f%% remaining)\n", num_cells_postQC, percent_cells_kept))

  exp_metrics$qc_nFeature_RNA_min <- nFeature_RNA.min
  exp_metrics$qc_nFeature_RNA_max <- nFeature_RNA.max
  exp_metrics$qc_nCount_RNA_min <- nCount_RNA.min
  exp_metrics$qc_nCount_RNA_max <- nCount_RNA.max
  exp_metrics$qc_num_cells_pre <- num_cells_preQC
  exp_metrics$qc_num_cells_post <- num_cells_postQC
  exp_metrics$qc_percent_cells_kept <- percent_cells_kept
  
  # reset sobj object
  sobj[["RNA"]]$scale.data <- NULL
  sobj <- FindVariableFeatures(sobj, verbose = FALSE)
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = FindNeighbors.dims, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = FindClusters.res, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = FindNeighbors.dims, verbose = FALSE)
  
  # Plot highly variable features per experiment
  # (avoids multi model clashes after integration)
  top10 <- head(VariableFeatures(sobj), 10)
  VFplot <- VariableFeaturePlot(sobj)
  print(LabelPoints(plot = VFplot, points = top10, repel = TRUE))
  
  # ---- Run doubletfinder helper function ----
  doubletfinder_res <- seurat_doubletfinder_20.2.2026_gh(
    sobj = sobj, 
    nobj = nobj,
    outdir = qcdir,
    FindNeighbors.dims = FindNeighbors.dims,
    seq_method = "droplet"
  )
  sobj <- doubletfinder_res$sobj

  exp_metrics$doubletfinder_seq_method <- doubletfinder_res$params$seq_method
  exp_metrics$doubletfinder_doublet_rate <- doubletfinder_res$params$doublet_rate
  exp_metrics$doubletfinder_best_pK <- doubletfinder_res$params$best.pK
  exp_metrics$doubletfinder_nExp_poi <- doubletfinder_res$params$nExp_poi
  exp_metrics$doubletfinder_nExp_poi_adj <- doubletfinder_res$params$nExp_poi.adj
  exp_metrics$doubletfinder_doublet_count <- doubletfinder_res$params$doublet_count
  
  # Store sobj in list
  sobj_list[[exp]] <- sobj
  metrics_rows[[exp]] <- as.data.frame(exp_metrics, stringsAsFactors = FALSE)
  cat("  - Saved", exp, "to sobj_list\n")
  
  # Close pdf
  dev.off()
  
}

# save qc metrics to a table
if (length(metrics_rows) > 0) {
  metrics_table <- bind_rows(metrics_rows)
  write.table(
    metrics_table,
    file = file.path(outdir, paste0(nobj, "_per_experiment_metrics.tsv")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

# ---- create merged object to save ----
if (length(sobj_list) > 1) {

  merged <- merge(
    x = sobj_list[[1]], 
    y = sobj_list[-1], 
    add.cell.ids = names(sobj_list), 
    project = nobj
  )

  # collapse sequencing runs into one layer
  merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

} else if (length(sobj_list) == 1) {
  
  # already one layer, just take the first object
  merged <- sobj_list[[1]]

} else {
  stop("ERROR: no Seurat objects found")
}

# ---- Save clustered and annotated sobj ----
cat(">>> Saving final Seurat object\n")
saveRDS(merged, file = file.path(outdir, paste0(nobj,"_merged.rds")))

# Save parameters to a text file
params <- list(
  script = "seurat_qc_260602_gh_cli.r",
  indir = indir,
  outdir = outdir,
  nobj = nobj,
  aligner = aligner,
  MAD_devs = MAD_devs,
  percent.mt.max = percent.mt.max,
  FindNeighbors.dims = max(FindNeighbors.dims), 
  FindClusters.res = FindClusters.res,
  final_cell_count = ncol(merged)
)

df <- data.frame(
  name  = names(params),
  value = unlist(params),
  row.names = NULL
)

write.table(
  df,
  file = file.path(outdir, paste0(nobj, "_parameters.txt")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

