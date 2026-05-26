#!/usr/bin/env Rscript
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

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

# ---- SoupX, QC, DoubletFinder ----
for (exp in names(crdir_filt)) {
  
  # Create Seurat object name
  cat(">>> Processing", exp, "\n")
  
  expdir <- file.path(outdir, paste0(exp, "_QC_out"))
  dir.create(expdir, showWarnings = FALSE, recursive = TRUE)
  expdir <- normalizePath(expdir)
  
  # Create pdf for plots (add time to prevent rewrite crash)
  pdf_path <- file.path(expdir, paste(nobj, exp, format(Sys.Date(), "%H-%M-%S"), "plots.pdf", sep = "_"))
  pdf(pdf_path, width = 10, height = 6)
  
  cat("  - Reading in data\n")
  sobj.filt <- Read10X(data.dir = crdir_filt[[exp]])
  
  if (length(crdir_raw[[exp]]) > 0) {
    
    sobj.raw <- Read10X(data.dir = crdir_raw[[exp]])
    
    # Create seurat object for soupx, with no filtering
    cat("  - Creating Seurat object\n")
    sobj <- CreateSeuratObject(
      counts = sobj.filt,
      project = nobj
    )
    
    # runs SCT, PCA, neighbours, clusters, UMAP (requires percent.mt)
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    sobj <- seurat_pipe_20.2.2026_gh(
      sobj = sobj, 
      FindNeighbors.dims = FindNeighbors.dims,
      FindClusters.res = FindClusters.res
    )
    
    # ---- Run soupx helper function ----
    # sobj has to have umap
    sobj <- seurat_soupx_23.2.2026_gh(
      sobj = sobj, 
      nobj = nobj,
      FindNeighbors.dims = FindNeighbors.dims,
      FindClusters.res = FindClusters.res,
      sobj.raw = sobj.raw, 
      sobj.filt = sobj.filt,
      outdir = expdir
    )
    
  } else {
    cat("   WARNING: no raw counts directory found for", exp, "continuing with filtered counts only\n")
    
    # Create seurat object with filtered counts only
    cat("  - Creating Seurat object\n")
    sobj <- CreateSeuratObject(
      counts = sobj.filt,
      min.cells = 3,
      min.features = 200,
      project = nobj
    )
    
  }
  
  # ---- Raw QC ----
  cat("  - Performing QC\n")
  
  # Refilter and add mito to new sobj
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  
  # Calculate thresholds using Mean Absolute Deviations (MADs)
  # https://bioconductor.org/books/3.15/OSCA.basic/quality-control.html#quality-control-outlier
  MAD_feats.min <- median(sobj$nFeature_RNA) - MAD_devs * mad(sobj$nFeature_RNA)
  nFeature_RNA.min <- max(200, MAD_feats.min) # clamp lower threshold (match to seurat object creation)
  nFeature_RNA.max <- median(sobj$nFeature_RNA) + MAD_devs * mad(sobj$nFeature_RNA)
  
  MAD_count.min <- median(sobj$nCount_RNA) - MAD_devs * mad(sobj$nCount_RNA)
  nCount_RNA.min <- max(0, MAD_count.min) # clamp lower threshold (feature filtering should capture cells with really low counts)
  nCount_RNA.max <- median(sobj$nCount_RNA) + MAD_devs * mad(sobj$nCount_RNA)
  
  #Visualize QC metrics with violins and scatters
  print(VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts"))
  
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
  
  p.mt <- VlnPlot(sobj, features = "percent.mt", layer = "counts") +
    geom_hline(yintercept = percent.mt.max, linetype = "dashed", color = "tomato") +
    ggtitle("percent.mt") +
    NoLegend()
  
  print(p.feat | p.count | p.mt)
  
  print(FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  
  # Remove cells that fail QC
  num_cells_preQC <- ncol(sobj)
  sobj <- subset(
    sobj, 
    subset = 
      nFeature_RNA > nFeature_RNA.min & 
      nFeature_RNA < nFeature_RNA.max &
      nCount_RNA > nCount_RNA.min &
      nCount_RNA < nCount_RNA.max &
      percent.mt < percent.mt.max
  )
  num_cells_postQC <- ncol(sobj)
  percent_cells_kept <- (num_cells_postQC/num_cells_preQC)*100
  cat(sprintf("  - Number of cells after QC: %i (%.2f%% remaining)\n", num_cells_postQC, percent_cells_kept))
  
  # save qc stats per experiment
  params <- list(
    nFeature_RNA.min = nFeature_RNA.min,
    nFeature_RNA.max = nFeature_RNA.max,
    nCount_RNA.min = nCount_RNA.min,
    nCount_RNA.max = nCount_RNA.max,
    num_cells_preQC = num_cells_preQC,
    num_cells_postQC = num_cells_postQC,
    percent_cells_kept = percent_cells_kept
  )
  
  df <- data.frame(name = names(params), value = unlist(params), row.names = NULL)
  write.table(df, file = file.path(expdir, paste(nobj,exp,"QC_stats.txt", sep="_")), quote = FALSE, sep = "\t", row.names = FALSE)
  
  # runs SCT, PCA, neighbours, clusters, UMAP
  sobj <- seurat_pipe_20.2.2026_gh(
    sobj = sobj, 
    FindNeighbors.dims = FindNeighbors.dims,
    FindClusters.res = FindClusters.res
  )
  
  # Plot highly variable features per experiment
  # (avoids multi model clashes after integration)
  top10 <- head(VariableFeatures(sobj), 10)
  VFplot <- VariableFeaturePlot(sobj)
  print(LabelPoints(plot = VFplot, points = top10, repel = TRUE))
  
  # ---- Run doubletfinder helper function ----
  sobj <- seurat_doubletfinder_20.2.2026_gh(
    sobj = sobj, 
    nobj = nobj,
    outdir = expdir,
    FindNeighbors.dims = FindNeighbors.dims,
    seq_method = "droplet"
  )
  
  # Store sobj in list
  sobj_list[[exp]] <- sobj
  cat("  - Saved", exp, "to sobj_list\n")
  
  # Close pdf
  dev.off()
  
}

# Create pdf for plots (add time to prevent rewrite crash)
pdf_path <- file.path(outdir, paste(nobj, format(Sys.Date(), "%H-%M-%S"), "plots.pdf", sep = "_"))
pdf(pdf_path, width = 10, height = 6)

# ---- Integrate data ----
if (length(sobj_list) > 1) {
  # runs merging, seurat pipeline, integration, and pipeline again
  merged <- seurat_integrate_24.2.2026_gh(
    sobj_list = sobj_list,
    nobj = nobj,
    FindClusters.res = FindClusters.res,
    FindNeighbors.dims = FindNeighbors.dims
  )
} else if (length(sobj_list) == 1) {
  # runs sctransform, pca, neighbours, clusters, umap
  merged <- seurat_pipe_20.2.2026_gh(
    sobj = sobj_list[[1]],
    FindNeighbors.dims = FindNeighbors.dims,
    FindClusters.res = FindClusters.res
  )
} else {
  stop("ERROR: no Seurat objects found")
}

# Report variable features and plot dimensionality reduction plots
print(VizDimLoadings(merged, dims = 1:2, reduction = "pca"))
print(DimPlot(merged, reduction = "pca", raster = FALSE) + NoLegend())
print(DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE))
print(DimPlot(merged, reduction = "umap", raster = FALSE))
print(DimPlot(merged, reduction = "umap", raster = FALSE, label = TRUE)+ NoLegend())

# ---- Plotting markers ----
# Highlight genes of interest/markers
all_markers_list  <- list(
  duct =  c("SOX9", "KRT19", "CFTR", "MUC1"), 
  acin =  c("CPA1", "PRSS1", "CTRC", "AMY2B"),
  vasc =  c("PDGFRB", "PECAM1", "VWF", "RGS5"), # endothelial, pericytes
  imm =   c("MS4A1", "LYZ", "PTPRC", "CD68"), # macrophages, T/B cells
  strom = c("COL1A1", "PDGFRA", "DCN", "COL3A1"), # fibroblasts, stellate
  endoc =  c("ISL1", "NEUROD1", "PDX1", "CHGA"),
  horm = c("INS", "GCG", "SST", "PPY", "GHRL"),
  empd =  c("MALAT1", "nFeature_RNA", "nCount_RNA", "percent.mt")
)

# filter all markers to remove warning messages
allowed <- union(rownames(merged), colnames(merged@meta.data))
all_markers_filt <- lapply(
  all_markers_list, 
  function(g) intersect(g, allowed)
)
all_markers_filt <- Filter(length, all_markers_filt)

# show marker expression per cluster
for (markers in all_markers_filt) {
  print(FeaturePlot(merged, features = markers, reduction = "umap", raster = FALSE))
  print(VlnPlot(merged, features = markers))
}

# dotplot, fix gene names
print(
  DotPlot(merged, features = all_markers_filt) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
)

# ---- Calculate cluster stats ----
# Number and percentage of cells per cluster
per_cluster_stats <- 
  table(Idents(merged)) %>% 
  as.data.frame() %>% 
  rename(cluster = Var1, cell_count = Freq) %>% 
  mutate(percentage = round((cell_count/ncol(merged))*100, 2))
write.table(
  per_cluster_stats,
  file = file.path(outdir, paste0(nobj, "_per_cluster_stats.txt")),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# Finding differentially expressed features
if (length(merged[["SCT"]]@SCTModel.list) > 1) {
  cat("  - Detected multi-model SCT object, prepping before FindAllMarkers()\n")
  merged <- PrepSCTFindMarkers(merged, verbose = FALSE)
}

merged.markers <- FindAllMarkers(merged, only.pos = TRUE, verbose = FALSE)
write.table(
  merged.markers, 
  file = file.path(outdir, paste0(nobj,"_allmarkers.txt")),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# heatmap top 5 DEGs per cluster 
top5 <- 
  merged.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()
print(DoHeatmap(merged, features = top5$gene) + NoLegend())

# Save significantly expressed genes
signif_markers <- 
  merged.markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05)
write.table(
  signif_markers,
  file = file.path(outdir, paste0(nobj,"_adjpval0.05_logfc1_markers.txt")),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# Save top 50 markers per cluster for GSEA/GO
ranked_top50 <- merged.markers %>%
  filter(p_val_adj < 0.05) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  mutate(rank = row_number()) %>% 
  slice_head(n = 50) %>%
  ungroup() %>% 
  select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)
write.table(
  ranked_top50,
  file = file.path(outdir, paste0(nobj,"_ranked_top50_markers.txt")),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# ---- Save clustered and annotated sobj ----
cat(">>> Saving final Seurat object\n")
dev.off()
saveRDS(merged, file = file.path(outdir, paste0(nobj,"_merged.rds")))

# Save parameters to a text file
params <- list(
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

