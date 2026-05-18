library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(tibble)
library(data.table)
library(DESeq2)
library(future)
library(pysch)

source("utils/malat1_function.R")

# Threaded math 
options(future.globals.maxSize = 110 * 1024^3) # 1024^3 = 1Gb

# ---- settings ----
FindClusters.res <- 0.4
FindNeighbors.dims <- 1:30

# ---- read in data ----
path <- here::here()
indir <- normalizePath(file.path(path, "../../processed/seurat_droplet_11-5-2026_gh_cli"))
outdir <- normalizePath(file.path(path, "../../processed/seurat_SC-vs-Human_2-5-2026_gh_integrated.v3"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data.dirs <- list.dirs(indir, recursive = FALSE)
rds.dirs <- list.files(data.dirs, pattern = ".rds", full.names = TRUE)
data.list <- lapply(rds.dirs, readRDS)

name.dirs <- sapply(data.dirs, basename)
name.dirs <- setdiff(name.dirs, "other")
names(data.list) <- sapply(strsplit(name.dirs, "_"), `[`, 1)

# ---- Create combined seurat object ----
data.list <- lapply(data.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x[["SCT"]] <- NULL
  
  # apply malat1 thresholding per object
  # https://github.com/BaderLab/MALAT1_threshold
  x <- NormalizeData(x, verbose = FALSE)
  norm_counts <- GetAssayData(x, assay = "RNA", layer = "data")["MALAT1",]
  threshold <- define_malat1_threshold_ggplot2(norm_counts)
  malat1_threshold <- norm_counts > threshold
  x$malat1_threshold <- malat1_threshold
  x$malat1_threshold <- factor(x$malat1_threshold, levels = c(TRUE, FALSE))
  DimPlot(x, group.by = "malat1_threshold")
  good_cells <- colnames(x)[malat1_threshold]
  x <- subset(x, cells = good_cells)
  
  return(x)
})

# Merge data for comparison 
merged <- merge(
  x = data.list[[1]],
  y = data.list[-1],
  add.cell.ids = names(data.list),
  project = "SC-vs-Human"
)

# extracting experiment from cell id
merged@meta.data$experiment <- sapply(strsplit(rownames(merged@meta.data), "_"), "[", 1)

# extracting condition
parts <- strsplit(rownames(merged@meta.data), "_")
merged@meta.data$condition <- ifelse(
  merged@meta.data$experiment == "akerman",
  ifelse(
    sapply(parts, "[", 2) == "D45",
    sapply(parts, function(x) paste(x[2:3], collapse = "_")),
    sapply(parts, function(x) paste(x[2:4], collapse = "_"))
  ),
  merged@meta.data$experiment
)

# extracting source
merged@meta.data <- merged@meta.data %>% 
  mutate(source = case_when(
    condition %in% c("D45_H1", "D39_ZKSCAN1_WT37") ~ "H1",
    condition %in% c("D39_SIM1_WT8", "D39_SIM1_WT18") ~ "H3",
    TRUE ~ condition
  ))
merged$source <- factor(merged$source, levels = c("bandesh", "fasolino", "kang", "H1", "H3"))

## create combined matrix then split by source
merged <- JoinLayers(merged)
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$source)

# ---- Start analysis ----
gc()
merged <- NormalizeData(merged, assay = "RNA", verbose = FALSE)
merged <- FindVariableFeatures(merged, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
gc()

merged <- RunPCA(merged, verbose = FALSE)
print(ElbowPlot(merged, ndims = 30))
merged <- FindNeighbors(merged, dims = FindNeighbors.dims, verbose = FALSE)
merged <- FindClusters(merged, resolution = 1, cluster.name = "unintegrated_clusters", verbose = FALSE)
merged <- RunUMAP(merged, dims = FindNeighbors.dims, reduction.name = "umap.unintegrated", verbose = FALSE)

# ---- Plotting ----
gene_markers <- c(
  "ISL1", "INS", "IAPP", "GCG", "ARX", "SST", "GHRL", "PPY","TPH1",
  "KRT19", # ductal
  "PRSS1", # acinar
  "PECAM1", # vasc
  "PTPRC", # all immune
  "COL1A1" # stromal
)
qc_markers <- c("MALAT1", "nFeature_RNA", "nCount_RNA", "percent.mt")
markers <- c(gene_markers, qc_markers)

for (m in markers) {
  p1 <- FeaturePlot(merged, reduction = "umap.unintegrated", features = m, raster = FALSE)
  p2 <- VlnPlot(merged, features = m) + NoLegend()
  print(p1 | p2)
}

print(DotPlot(merged, features = gene_markers))
print(RidgePlot(merged, features = qc_markers, ncol = length(qc_markers)))

print(DimPlot(merged, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", raster = FALSE))
print(DimPlot(merged, reduction = "umap.unintegrated", group.by = "source", raster = FALSE))
print(DimPlot(merged, reduction = "umap.unintegrated", label = TRUE, raster = FALSE) + NoLegend())

# ---- Removing empty droplets ----
empty_clusters <- c(0, 10, 20, 27, 32, 34, 35, 38)
DimPlot(
  merged,
  reduction = "umap.unintegrated",
  raster = FALSE,
  cells.highlight = WhichCells(merged, idents = empty_clusters),
  cols = "grey",
  cols.highlight = "tomato",
  pt.size = 0.1,
  sizes.highlight = 0.05   # <-- make highlighted points smaller
)

## remove empty droplets and try UMAPing from fresh (avoid crash)
integrated <- subset(merged, !unintegrated_clusters %in% empty_clusters)
integrated[["RNA"]]$scale.data <- NULL
gc()

# re run the pipline
integrated <- NormalizeData(integrated, assay = "RNA", verbose = FALSE)
integrated <- FindVariableFeatures(integrated, verbose = FALSE)
integrated <- ScaleData(integrated, verbose = FALSE)
gc()

integrated <- RunPCA(integrated, verbose = FALSE)
print(ElbowPlot(integrated, ndims = 30))

# harmony integration, two different resolutions for picking out smaller clusters
integrated <- IntegrateLayers(integrated, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = FindNeighbors.dims, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = FindClusters.res, cluster.name = "harmony_clusters", verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 8, cluster.name = "harmony_clusters_8", verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = FindNeighbors.dims, reduction.name = "umap.harmony", verbose = FALSE)

# # rpca integration (squashes too much)
# integrated <- IntegrateLayers(integrated, method = RPCAIntegration, new.reduction = "rpca", verbose = FALSE) 
# integrated <- FindNeighbors(integrated, reduction = "rpca", dims = FindNeighbors.dims, verbose = FALSE)
# integrated <- FindClusters(integrated, resolution = FindClusters.res, cluster.name = "rpca_clusters", verbose = FALSE)
# integrated <- RunUMAP(integrated, reduction = "rpca", dims = FindNeighbors.dims, reduction.name = "umap.rpca", verbose = FALSE)

for (reduc in c("harmony")) {
  
  reduc_name <- paste("umap", reduc, sep = ".")
  clust_name <- paste(reduc, "clusters", sep = "_")
  
  # plot new data
  for (m in markers) {
    p1 <- FeaturePlot(integrated, reduction = reduc_name, features = m, raster = FALSE)
    ggsave(file.path(outdir, paste0("FeaturePlots_", reduc), paste0(m,".png")), plot = p1, width = 11, height = 10)
    p2 <- VlnPlot(integrated, group.by = clust_name, features = m) + NoLegend()
    print(p1 | p2)
    ggsave(file.path(outdir, paste0("FeatureVlnPlots_", reduc), paste0(m, ".png")), width = 14, height = 9)
  }

  print(DotPlot(integrated, group.by = clust_name, features = gene_markers))
  print(RidgePlot(integrated, group.by = clust_name, features = qc_markers, ncol = length(qc_markers)))

  print(DimPlot(integrated, reduction = reduc_name, group.by = clust_name, raster = FALSE))
  print(DimPlot(integrated, reduction = reduc_name, group.by = "source", raster = FALSE))
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_source.png", clust_name)), height = 10, width = 11)
  print(DimPlot(integrated, reduction = reduc_name, group.by = clust_name, raster = FALSE, label = TRUE) + NoLegend())
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_numbered.png", clust_name)), height = 10, width = 11)
  
}
  
# ---- labelling clusters ----
integrated@meta.data <- integrated@meta.data %>%
  mutate(cell_type = case_when(
    
    # high-res first, to overide 
    harmony_clusters_8 %in% c(118) ~ "delta",
    harmony_clusters_8 %in% c(108) ~ "epsilon",
    
    # coarse clusters last, for anything left 
    harmony_clusters %in% c(0, 3, 8, 17) ~ "beta",
    harmony_clusters %in% c(1, 4, 6, 14) ~ "alpha",
    harmony_clusters %in% c(2) ~ "ductal",
    harmony_clusters %in% c(9) ~ "ec",
    harmony_clusters %in% c(7) ~ "acinar",
    harmony_clusters %in% c(11) ~ "vasc",
    harmony_clusters %in% c(10) ~ "stromal",
    harmony_clusters %in% c(12) ~ "ppy",
    harmony_clusters %in% c(5) ~ "delta",
    harmony_clusters %in% c(15 , 16) ~ "immune",
    
    TRUE ~ "other"
  ))
integrated$cell_type <- factor(integrated$cell_type)

DimPlot(integrated, reduction = "umap.harmony", group.by = "cell_type", label = TRUE, raster = FALSE)
DotPlot(integrated, features = gene_markers, group.by = "cell_type")

# annotate cell type origin
primary <- c("bandesh", "fasolino", "kang")
integrated$beta_cell <- case_when(
  integrated$cell_type == "beta" & integrated$experiment == "akerman" ~ "sc-beta",
  integrated$cell_type == "beta" & integrated$experiment %in% primary ~ "primary-beta",
  TRUE ~ "other"
)
integrated$alpha_cell <- case_when(
  integrated$cell_type == "alpha" & integrated$experiment == "akerman" ~ "sc-alpha",
  integrated$cell_type == "alpha" & integrated$experiment %in% primary ~ "primary-alpha",
  TRUE ~ "other"
)
integrated$delta_cell <- case_when(
  integrated$cell_type == "delta" & integrated$experiment == "akerman" ~ "sc-delta",
  integrated$cell_type == "delta" & integrated$experiment %in% primary ~ "primary-delta",
  TRUE ~ "other"
)

# find sc-alpha cells
DimPlot(
  integrated,
  reduction = "umap.harmony",
  raster = FALSE,
  cells.highlight = WhichCells(
    integrated, 
    expression = alpha_cell == "sc-alpha"
  )
)

# find sc-delta cells
DimPlot(
  integrated,
  reduction = "umap.harmony",
  raster = FALSE,
  cells.highlight = WhichCells(
    integrated, 
    expression = delta_cell == "sc-delta"
  )
)

# group H1 and H3 beta cell sub populations for cleaner plotting
integrated$beta_cell_condition <- ifelse(integrated$cell_type == "beta", as.character(integrated$source), NA)

# factorise beta_cell_type column for cleaner plotting
integrated$beta_cell_condition <- factor(
  integrated$beta_cell_condition,
  levels = c("bandesh", "fasolino", "kang", "H1", "H3")
)

print(DimPlot(
  integrated,
  reduction = "umap.harmony",
  group.by = "beta_cell_condition",
  raster = FALSE
))

# create human vs h1 vs h3 for pairwise analysis and lickert contamination
integrated$cell_origin <- ifelse(
  integrated$source %in% primary, "primary", 
  as.character(integrated$source)
)
integrated$cell_origin <- factor(integrated$cell_origin)

# ---- beta cell specific analysis ----
# subset beta cells first to remove NA
beta <- subset(integrated, cell_type == "beta")
beta[["RNA"]]$scale.data <- NULL
gc()

# renormalise after subsetting, safe to join because we never scale
beta[["RNA"]] <- JoinLayers(beta[["RNA"]])
beta <- NormalizeData(beta, assay = "RNA", verbose = FALSE)
beta <- FindVariableFeatures(beta, nfeatures = 3000, verbose = FALSE)
gc()

beta <- SCTransform(beta, assay = "RNA", verbose = FALSE)

print(VlnPlot(
    beta, 
    assay = "RNA", 
    layer = "data", 
    group.by = "source", 
    features = "INS",
    pt.size = 0,
    cols = c(
      "bandesh" = "#990000",
      "fasolino" = "#CC0000",
      "kang" = "#FF6666",
      "H1" = "#38BDF8",
      "H3" = "#003e86"
    )
  ) +
  labs(title = t) +
  NoLegend()
)

# calculate pairwise significance values inbetween each group
pairs <- combn(levels(beta$beta_cell_origin), 2, simplify = FALSE)
done <- lapply(pairs, function(p) {
  FindMarkers(
    beta,
    assay = "RNA",
    layer = "data",
    features = "INS",
    verbose = FALSE,
    group.by = "beta_cell_origin",
    ident.1 = p[1],
    ident.2 = p[2],
    logfc.threshold = 0
  )
})
names(done) <- vapply(pairs, function(p) paste(p, collapse = "-vs-"), "")

# Actual differences between INS expr across conditions
avg.ins.expr <- AverageExpression(
  beta, 
  group.by = "beta_cell_condition", 
  assays = "RNA", 
  layer = "data",
  features = "INS", 
  verbose = FALSE
)$RNA

# AverageExpression un-logs values before averaging, re-log
avg.ins.expr_df <- avg.ins.expr %>%
  as.matrix() %>%              
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "condition", values_to = "avg.INS") %>% 
  select(-(gene)) %>% 
  mutate(avg.INS_ln = log1p(avg.INS))

write.table(
  avg.ins.expr_df,
  file.path(outdir, "avg-ins-expr_source.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# plot genes of interest
beta_markers <- c(
  "PDX1", "NKX6-1", "MAFA", "MNX1", "GLIS3", "HNF1A",  # beta cell
  "NEUROD1", "INSM1", "ISL1", "NKX2-2", "PAX6", "RFX6", "MYT1", "MYT1L", # beta cell/ neuronal
  "MAFB", "PAX4", "SIX3", "FOXO1", "FOXA2", "SOX9" # alpha cell / other 
)

# reverse order for plotting
beta$beta_cell_condition_rev <- factor(
  beta$beta_cell_condition,
  levels = rev(levels(beta$beta_cell_condition))
)

DotPlot(beta, assay = "RNA", features = beta_markers, group.by = "beta_cell_condition_rev") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(
  file.path(outdir, "DotPlots", "DotPlot_beta-markers_combined.rev.png"),
  width = 14, height = 4
)

VlnPlot(beta, assay = "RNA", features = beta_markers, group.by = "beta_cell_condition_rev", stack = TRUE) +
  NoLegend() 
ggsave(
  file.path(outdir, "VlnPlots", "VlnPlot_beta-markers_combined.rev.png"),
  width = 14, height = 4
)

# Load in lickert contaminants for 
lickert_contaminants <- fread(
  file.path(indir, "../..", "raw_data", "Lickert_contaminants.csv"),
  select = c("gene", "category"),
  verbose = FALSE
)

# for SBCs and PBCs, calculate pct expressed for each contaminant gene
for (bc in unique(beta$beta_cell_origin)) {
  
  beta.bc <- subset(beta, beta_cell_origin == bc)
  genes <- intersect(lickert_contaminants$gene, rownames(beta.bc))
  lc_expr <- GetAssayData(beta.bc[genes, ], assay = "RNA", layer = "data")
  
  # set pct thresholds
  for (pct in c(0, 0.5, 1, 2)) {
    
    lc_expr.pct <- rowMeans(lc_expr > pct) * 100
    col_name <- paste("pct", pct, bc, sep = "_")
    lickert_contaminants[[col_name]] <- lc_expr.pct[lickert_contaminants$gene]
    
  }
  remove(beta.bc, genes, lc_expr, lc_expr.pct, col_name)
}
lickert_contaminants <- lickert_contaminants[order(lickert_contaminants$category),]

# save and plot results
write.table(
  lickert_contaminants,
  file.path(outdir, "pheatmaps", "lickert-contaminants_SBC-vs-PBC.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

mat <- lickert_contaminants |>
  column_to_rownames("gene") |>
  select(starts_with("pct")) |>
  as.matrix()

annotation_row <- data.frame(Category = lickert_contaminants$category)
rownames(annotation_row) <- rownames(mat)

annotation_col <- data.frame(beta_cell = sub(".*\\_", "", colnames(mat)))
rownames(annotation_col) <- colnames(mat)
  
ann_colours <- list(
  beta_cell = c(
    primary = "#CC0000",
    H1 = "#38BDF8",
    H3 = "#003e86" 
  )
)

pheatmap::pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = ann_colours,
  display_numbers = TRUE,
  #fontsize = 15,
  filename = file.path(outdir, "pheatmaps", "pheatmap_pct.expr.lickert_contaminants.png"),
  width = 7,
  height = 12
)

# ---- Pearson correlation coefficients (pcc) ----
# create new subdivisions for pcc analysis
integrated$pcc_clusters <- case_when(
  integrated$cell_type == "beta" ~ integrated$beta_cell_condition,
  integrated$cell_type == "alpha" ~ integrated$alpha_cell,
  integrated$cell_type == "delta" ~ integrated$delta_cell,
  TRUE ~ integrated$cell_type
)

integrated@meta.data <- integrated@meta.data %>%
  mutate(pcc_clusters.reps = case_when(
      cell_type == "beta" ~ paste("beta", condition, sep = "_"),
      cell_type == "alpha" ~ paste("alpha", as.character(cell_origin), sep = "_"),
      cell_type == "delta" ~ delta_cell,
      TRUE ~ cell_type
  ))

# factorise pcc_cluster column for cleaner plotting
integrated$pcc_clusters <- factor(
  integrated$pcc_clusters,
  levels = c(
    "bandesh", "fasolino", "kang", 
    "H1", "H3", 
    "primary-alpha", "sc-alpha", 
    "primary-delta", "sc-delta",
    "ductal"
  )
)

integrated$pcc_clusters.reps <- factor(
  integrated$pcc_clusters.reps,
  levels = c(
    "beta_bandesh", "beta_fasolino", "beta_kang", 
    "beta_D39_ZKSCAN1_WT37", "beta_D45_H1", 
    "beta_D39_SIM1_WT8", "beta_D39_SIM1_WT18", 
    "alpha_primary", "alpha_H1", "alpha_H3", 
    "primary-delta", "sc-delta",
    "ductal", "acinar"
  )
)

DimPlot(integrated, reduction = "umap.harmony", group.by = "pcc_clusters", raster = FALSE)
DimPlot(integrated, reduction = "umap.harmony", group.by = "pcc_clusters.reps", raster = FALSE)

# create bar plot to show counts per population
ggplot(integrated@meta.data, aes(
    x = factor(pcc_clusters.reps, levels = rev(levels(pcc_clusters.reps))), fill = pcc_clusters.reps)
  ) +
    geom_bar(show.legend = FALSE) +
    coord_flip() +
    theme_classic() + 
    expand_limits(y = max(table(integrated@meta.data$pcc_clusters.reps)) * 1.1) +
    geom_text(
      stat = "count", 
      aes(label = after_stat(count)), 
      hjust = -0.2
    )

# subset endocrine and get SCT HVFs
endocrine <- subset(integrated, cell_type %in% c("beta", "alpha", "delta", "ppy", "epsilon"))
endocrine <- SCTransform(endocrine, assay = "RNA", verbose = FALSE)
gc()

# SCTransform for HVFs on a per model basis
integrated <- SCTransform(integrated, assay = "RNA", verbose = FALSE)
gc()

# join integrated layers for expr matrix agg/avg
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

# build a function for data generation
make_pcc_step <- function(step, pcc, pcc_clust, beta) {
  
  if (step == "top3000-beta-vfs") {
    
    vfs <- VariableFeatures(beta, method = "sctransform", nfeatures = 3000)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else if (step == "top3000-vfs") {
    
    vfs <- VariableFeatures(pcc, method = "sctransform", nfeatures = 3000)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else if (step == "all-genes") {
    
    vfs <- rownames(pcc)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE) 
    
  } else if (step == "top3000-endo-vfs") {
    
    vfs <- VariableFeatures(endocrine, method = "sctransform", nfeatures = 3000)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else if (step == "top1000-vfs") {
    
    vfs <- VariableFeatures(pcc, method = "sctransform", nfeatures = 1000)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else {
    stop("Unknown step")
  }
  
  ## colours
  groups <- factor(colnames(agg$RNA), levels = levels(pcc@meta.data[[pcc_clust]]))
  annotation_col <- data.frame(Group = groups)
  rownames(annotation_col) <- colnames(agg$RNA)
  
  ## correlations
  corr.res.agg <- corr.test(as.matrix(agg$RNA), method = "pearson", adjust = "BH")
  corr.res.avg <- corr.test(as.matrix(avg$RNA), method = "pearson", adjust = "BH")
  
  agg.corr2 <- (corr.res.agg$r)^2
  avg.corr2 <- (corr.res.avg$r)^2
  
  list(
    title = step,
    features = vfs,
    agg = agg,
    avg = avg,
    agg.corr2 = agg.corr2,
    avg.corr2 = avg.corr2,
    agg.p_adj = corr.res.agg$p,
    avg.p_adj = corr.res.avg$p,
    annotation_col = annotation_col
  )
}

# build a function for colour matching
make_pcc_colors <- function(pcc, pcc_clust) {
  
  lvls <- levels(pcc@meta.data[[pcc_clust]])
  gg_cols <- scales::hue_pal()(length(lvls))
  
  list(
    levels = lvls,
    annotation_colors = list(
      Group = setNames(gg_cols, lvls)
    )
  )
}

# build data for heatmaps
pcc_results <- list()

steps <- c(
  "top3000-beta-vfs",
  "top3000-vfs",
  "all-genes",
  "top3000-endo-vfs",
  "top1000-vfs"
)

for (pcc_clust in c("pcc_clusters", "pcc_clusters.reps")) {
  
  outdir_pheatmap <- file.path(outdir, "pheatmaps", pcc_clust)
  dir.create(outdir_pheatmap, recursive = TRUE, showWarnings = FALSE)
  
  colors <- make_pcc_colors(integrated, pcc_clust)
  
  pcc_results[[pcc_clust]] <- list(
    outdir = outdir_pheatmap,
    colors = colors,
    steps  = lapply(
      steps, make_pcc_step,
      pcc = integrated,
      pcc_clust = pcc_clust,
      beta = beta
    )
  )
  
  names(pcc_results[[pcc_clust]]$steps) <- steps
}

# plot heatmaps
for (pcc_clust in names(pcc_results)) {
  
  if (pcc_clust == "pcc_clusters") {fon <- 15; wid <- 9; hei <- 7} 
  else {fon <- 15; wid <- 11; hei <- 9}
  
  res <- pcc_results[[pcc_clust]]
  
  for (step in names(res$steps)) {
    
    obj <- res$steps[[step]]
    
    obj$agg.corr2[lower.tri(obj$agg.corr2)] <- NA 
    pheatmap::pheatmap(
      obj$agg.corr2,
      annotation_col = obj$annotation_col,
      annotation_colors = res$colors$annotation_colors,
      color = colorRampPalette(c("white", "red"))(100),
      breaks = seq(0, 1, length.out = 101),   # fixed scale,
      main = sprintf("pearson r2 agg expr (%s)", obj$title),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      fontsize = fon,
      filename = file.path(res$outdir, sprintf("pheatmap_agg-corr2_%s.png", obj$title)),
      width = wid, height = hei,
      na_col = "white"
    )
    
    #obj$agg.p_adj[lower.tri(obj$agg.p_adj)] <- NA 
    pheatmap::pheatmap(
      obj$agg.p_adj,
      annotation_col = obj$annotation_col,
      annotation_colors = res$colors$annotation_colors,
      color = colorRampPalette(c("white", "blue"))(100),
      breaks = seq(0, 0.05, length.out = 101),
      main = sprintf("bh p-adj agg expr (%s)", obj$title),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      fontsize = fon,
      filename = file.path(res$outdir, sprintf("pheatmap_agg-padj_%s.png", obj$title)),
      width = wid, height = hei,
      na_col = "white"
    )
    
    obj$avg.corr2[lower.tri(obj$avg.corr2)] <- NA
    pheatmap::pheatmap(
      obj$avg.corr2,
      annotation_col = obj$annotation_col,
      annotation_colors = res$colors$annotation_colors,
      color = colorRampPalette(c("white", "red"))(100),
      breaks = seq(0, 1, length.out = 101),   # fixed scale
      main = sprintf("pearson r2 avg expr (%s)", obj$title),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      fontsize = fon,
      filename = file.path(res$outdir, sprintf("pheatmap_avg-corr2_%s.png", obj$title)),
      width = wid, height = hei,
      na_col = "white"
    )
    
    #obj$avg.p_adj[lower.tri(obj$avg.p_adj)] <- NA
    pheatmap::pheatmap(
      obj$avg.p_adj,
      annotation_col = obj$annotation_col,
      annotation_colors = res$colors$annotation_colors,
      color = colorRampPalette(c("white", "blue"))(100),
      breaks = seq(0, 0.05, length.out = 101),
      main = sprintf("bh p-adj avg expr (%s)", obj$title),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      fontsize = fon,
      filename = file.path(res$outdir, sprintf("pheatmap_avg-padj_%s.png", obj$title)),
      width = wid, height = hei,
      na_col = "white"
    )
    
  }
}

# helper: flatten upper triangle
flatten_corr <- function(mat) {mat[upper.tri(mat, diag = FALSE)]}

# main function: build meta-correlation matrix
make_meta_corr <- function(step_list, slot = c("agg.corr2", "avg.corr2")) {
  
  slot <- match.arg(slot)
  
  # extract + vectorise
  vecs <- lapply(step_list, function(x) {
    flatten_corr(x[[slot]])
  })
  
  # combine into matrix (columns = steps)
  mat <- do.call(cbind, vecs)
  colnames(mat) <- names(step_list)
  
  # compute correlation between steps
  meta_corr2 <- cor(mat, method = "pearson", use = "pairwise.complete.obs")^2
  
  return(meta_corr2)
}


# calculate correlation matrices
meta_results <- list()

for (pcc_clust in names(pcc_results)) {
  
  steps_list <- pcc_results[[pcc_clust]]$steps
  
  meta_results[[pcc_clust]] <- list(
    agg_meta = make_meta_corr(steps_list, "agg.corr2"),
    avg_meta = make_meta_corr(steps_list, "avg.corr2")
  )
}

# plot results
for (pcc_clust in names(meta_results)) {
  
  pheatmap::pheatmap(
    meta_results[[pcc_clust]]$agg_meta,
    main = paste(pcc_clust, "- AGG meta-correlation"),
    display_numbers = TRUE,
    color = colorRampPalette(c("white", "red"))(100),
    filename = file.path(outdir, "pheatmaps", sprintf("pheatmap_meta-agg_%s.png", pcc_clust)),
    width = 6, height = 5
  )
  
  pheatmap::pheatmap(
    meta_results[[pcc_clust]]$avg_meta,
    main = paste(pcc_clust, "- AVG meta-correlation"),
    display_numbers = TRUE,
    color = colorRampPalette(c("white", "red"))(100),
    filename = file.path(outdir, "pheatmaps", sprintf("pheatmap_meta-avg_%s.png", pcc_clust)),
    width = 6, height = 5
  )

}


# calculate expression sets for each pcc cluster
integrated$expr.sets <- case_when(
  pcc$cell_type == "beta" ~ pcc$beta_cell,
  pcc$cell_type == "alpha" ~ pcc$alpha_cell,
  pcc$cell_type == "delta" ~ pcc$delta_cell,
  TRUE ~ NA
)

avg.expr.sets <- AverageExpression(pcc, group.by = "expr.sets", assays = "RNA", layer = "data", verbose = FALSE)$RNA
colnames(avg.expr.sets) <- sub("-", "_", colnames(avg.expr.sets))

for (coln in colnames(avg.expr.sets)) {
  avg.expr.sets %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    select(gene, coln) %>% 
    filter(.data[[coln]] > 0) %>% 
    arrange(desc(.data[[coln]])) %>% 
    write.table(
      file.path(outdir, "pcc-clusters_avg-expr-sets", paste0("avg-expr-sets_linear_", coln, ".rnk")),
      quote = FALSE,
      col.names = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
}

# ---- Perform DEA ----
if (length(integrated[["SCT"]]@SCTModel.list) > 1) {
  integrated <- PrepSCTFindMarkers(integrated, verbose = FALSE)
}

integrated.markers <- FindAllMarkers(integrated, group.by = "cell_type", only.pos = TRUE, verbose = FALSE)
write.table(
  integrated.markers,
  file.path(outdir, "FindAllMarkers_cell-type_output.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# find DEGs for stem-cell beta cells vs primary beta cells
integrated$beta_cell <- factor(integrated$beta_cell)
sbc.vs.pbc <- FindMarkers(
  integrated, 
  group.by = "beta_cell",
  ident.1 = "sc-beta", 
  ident.2 = "primary-beta"
)

# rename columns for clarity and create gene column
sbc.vs.pbc <- dplyr::rename(sbc.vs.pbc, pct.sbc = pct.1, pct.pbc = pct.2)
sbc.vs.pbc$gene <- row.names(sbc.vs.pbc)

write.table(
  sbc.vs.pbc,
  file.path(outdir, "FindMarkers_SBC-vs-PBC_output.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# create ranked gene set for gsea
sbc.vs.pbc %>% 
  select(gene, avg_log2FC) %>% 
  arrange(desc(avg_log2FC)) %>% 
  write.table(
    file.path(outdir, "FindMarkers_SBC-vs-PBC_output.rnk"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

# find markers for sc and primary beta cells
integrated.markers.beta_cell <- FindAllMarkers(integrated, group.by = "beta_cell", verbose = FALSE)
integrated.markers.beta_cell <- dplyr::rename(integrated.markers.beta_cell, beta_cell = cluster)

write.table(
  integrated.markers.beta_cell,
  file.path(outdir, "FindAllMarkers_beta-cell_output.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# create preranked gene set lists for markers of SBCs and PBCs
for (bc in unique(integrated.markers.beta_cell$beta_cell)) {
  if (bc == "other") next
  integrated.markers.beta_cell %>% 
    filter(beta_cell == bc) %>% 
    select(gene, avg_log2FC) %>% 
    arrange(desc(avg_log2FC)) %>% 
    write.table(
      file.path(outdir, sprintf("FindAllMarkers_beta-cell_output.%s.rnk", bc)),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
}

# create pseudo bulk count matrix for deseq2
condition.expr <- AggregateExpression(
  beta,
  group.by = "condition",
  assay = "RNA",
  layers = "counts",
  normalization.method = NULL,
  scale.factor = NULL,
  margin = NULL,
  verbose = FALSE
)$RNA

deseq2.info <- data.frame(
  sample = colnames(condition.expr),
  condition = c("primary-beta", replicate(4, "sc-beta"), replicate(2, "primary-beta")),
  replicate = c("pri1", "sc1", "sc2", "sc3", "sc4", "pri2", "pri3")
)
rownames(deseq2.info) <- deseq2.info$sample

dds <- DESeqDataSetFromMatrix(
  countData = round(condition.expr), 
  colData = deseq2.info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "sc-beta", "primary-beta")
)

vsd <- vst(dds)
plotPCA(vsd, intgroup = "condition")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df %>%
  select(gene, log2FoldChange, padj) %>%
  filter(!is.na(padj)) %>%
  mutate(
    neglog10p = -log10(pmax(padj, 1e-300)),
    neglog10p.sign = neglog10p * ifelse(log2FoldChange >= 0, 1, -1)
  ) %>%
  select(gene, neglog10p.sign) %>%
  arrange(desc(neglog10p.sign)) %>%
  write.table(
    file.path(outdir, "DESeq2_SBC-vs-PBC.neg-log10-sign.rnk"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

# create ranked linear expression of all genes in SBCs
## aggregate expression
sbc.agg.expr <- AggregateExpression(beta, group.by = "beta_cell", assay = "RNA", verbose = FALSE)$RNA
sbc.agg.expr_df <- sbc.agg.expr %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rename("sc-beta" = "sc_beta") %>% 
  rownames_to_column(var = "gene") %>%
  select(gene, sc_beta) %>% 
  arrange(desc(sc_beta))

write.table(
  sbc.agg.expr_df,
  file.path(outdir, "AggregateExpression_SBCs.rnk"),
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

sbc.agg.expr_df %>% 
  mutate(sc_beta_int = round(sc_beta)) %>% 
  select(gene, sc_beta_int) %>% 
  write.table(
    file.path(outdir, "AggregateExpression_SBCs.int.rnk"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

## average expression
sbc.avg.expr <- AverageExpression(beta, group.by = "beta_cell", assay = "RNA", verbose = FALSE)$RNA
sbc.avg.expr_df <- sbc.avg.expr %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rename("sc-beta" = "sc_beta") %>% 
  rownames_to_column(var = "gene") %>%
  select(gene, sc_beta) %>% 
  arrange(desc(sc_beta)) %>% 
  mutate(sc_beta_log = log1p(sc_beta))

sbc.avg.expr_df %>% 
  select(gene, sc_beta) %>% 
  write.table(
    file.path(outdir, "AverageExpression_SBCs.rnk"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

sbc.avg.expr_df %>% 
  select(gene, sc_beta_log) %>% 
  write.table(
    file.path(outdir, "AverageExpression_SBCs.log1p.rnk"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
