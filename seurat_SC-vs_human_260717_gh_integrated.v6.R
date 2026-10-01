# check conda env has been set up correctly
Sys.getenv("CONDA_PREFIX")
Sys.getenv("LD_LIBRARY_PATH")
Sys.getenv("LD_PRELOAD")
reticulate::py_config()

library(reticulate)
py_run_string("import pyexpat; print('success')")
system(sprintf("cat /proc/%d/maps | grep expat", Sys.getpid()))

# load in libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(DESeq2)
library(psych)
library(SeuratWrappers)
library(ggVennDiagram)
library(ggrepel)

# increase max worker size
options(future.globals.maxSize = 160 * 1024^3) # 1024^3 = 1Gb

# ---- default settings ----
FindClusters.res <- 1
FindNeighbors.dims <- 1:30
plot_cols <- c("lightgrey", "red")

# ---- data ----
# set paths
path <- normalizePath(file.path(here::here(), "../.."))
indir <- file.path(path, "processed/seurat_qc_260722_gh_cli")
outdir <- file.path(path, "processed/seurat_SC-vs-Human_260717_gh_integrated.v6")

# create outdirs (prevents ggsave interrupting console)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "VlnPlots"), showWarnings = FALSE)
dir.create(file.path(outdir, "FeaturePlots"), showWarnings = FALSE)
dir.create(file.path(outdir, "DimPlots"), showWarnings = FALSE)
dir.create(file.path(outdir, "DotPlots"), showWarnings = FALSE)
dir.create(file.path(outdir, "pheatmaps"), showWarnings = FALSE)
dir.create(file.path(outdir, "dendrograms"), showWarnings = FALSE)
dir.create(file.path(outdir, "pca"), showWarnings = FALSE)
dir.create(file.path(outdir, "volcanos"), showWarnings = FALSE)

# read in data
data.dirs <- list.dirs(indir, recursive = FALSE)
rds.dirs <- list.files(data.dirs, pattern = ".rds", full.names = TRUE)
data.list <- lapply(rds.dirs, readRDS)

# label data
name.dirs <- sapply(data.dirs, basename)
name.dirs <- setdiff(name.dirs, "other")
names(data.list) <- sapply(strsplit(name.dirs, "_"), `[`, 1)

# ---- Create combined seurat object ----
# Merge data for comparison 
merged <- merge(
  x = data.list[[1]],
  y = data.list[-1],
  add.cell.ids = names(data.list),
  project = "SC-vs-Human"
)

# clear space and unused variables
remove(data.list)
keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "pANN")
merged@meta.data <- merged@meta.data[, keep, drop = FALSE]
gc()

# overwrite orig.ident with lab
merged$orig.ident <- sapply(strsplit(rownames(merged@meta.data), "_"), "[", 1)
merged$orig.ident <- factor(
  merged$orig.ident, 
  levels = c(
    "bandesh",
    "fasolino",
    "kang",
    "akerman",
    "balboa",
    "dadheech",
    "kelley",
    "wu"
  ) 
)

orig.ident_cols <- c(
  "bandesh" = "grey30",
  "fasolino" = "grey50",
  "kang" = "grey70",
  "akerman" = "mediumpurple",
  "balboa" = "slateblue",
  "dadheech" = "thistle",
  "kelley" = "maroon",
  "wu" = "hotpink"
)

# mark source, primary or stem-cell
primary <- c("bandesh", "fasolino", "kang")
merged$source <- ifelse(
  merged$orig.ident %in% primary,
  "primary",
  "stem-cell"
)

# extracting replicates
merged$sample <- sapply(strsplit(rownames(merged@meta.data), "_"), '[', 2)

# combine then split to standardize RNA layers (eg SIX3)
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample)

# cell_line
merged@meta.data <- merged@meta.data %>% 
  mutate(cell_line = case_when(
    grepl("^(D39-Z|D45|SRX101888|SRX25096888)", sample) ~ "H1",
    grepl("^D39-S", sample) ~ "H3",
    grepl("^SRX25096890", sample) ~ "HS980",
    grepl("^(1|2)", sample) ~ "HUES8",
    grepl("^SRX2452304", sample) ~ "iPSC",
    TRUE ~ "primary"
  ))
merged$cell_line <- factor(
  merged$cell_line, 
  levels = c(
    "primary", "H1", "H3", "HS980", "HUES8", "iPSC"
  )
)

# stage.week
merged@meta.data <- merged@meta.data %>% 
  mutate(stage_week = case_when(
    orig.ident == "wu" ~ "s6_w0",
    orig.ident %in% c("dadheech", "kelley") ~ "s6_w1",
    sample %in% c("SRX10188829", "SRX10188830", "SRX10188831") ~ "s7_w0",
    orig.ident == "akerman" ~ "s7_w2",
    sample %in% c("SRX10188832", "SRX10188833", "SRX10188834") ~ "s7_w3",
    sample %in% c("SRX10188835", "SRX10188836") ~ "s7_w6",
    TRUE ~ "primary"
  ))
merged$stage_week <- factor(
  merged$stage_week,
  levels = c(
    "primary", "s7_w6", "s7_w3", "s7_w2", "s7_w0", "s6_w1", "s6_w0"
  )
)

# make grouped bar charts for cell types per lab
merged$id <- ifelse(
  merged$orig.ident == "akerman",
  paste(merged$orig.ident, merged$sample, sep = "."),
  ifelse(
    merged$orig.ident %in% primary,
    paste(merged$orig.ident, merged$source, sep = "."),
    paste(merged$orig.ident, merged$cell_line, merged$stage_week, sep = ".")
  )
)
merged$id <- factor(
  merged$id,
  levels = c(
    "bandesh.primary", "fasolino.primary", "kang.primary",
    "akerman.D39-ZKSCAN1-WT37", "akerman.D45-H1",
    "akerman.D39-SIM1-WT8", "akerman.D39-SIM1-WT18",
    "balboa.H1.s7_w0", "balboa.H1.s7_w3", "balboa.H1.s7_w6",
    "kelley.HUES8.s6_w1",
    "wu.H1.s6_w0", "wu.HS980.s6_w0",
    "dadheech.iPSC.s6_w1"
  )
)

# create reversed id for plotting
merged$id.rev <- factor(merged$id,levels = rev(levels(merged$id)))

# set colours to use for id
id_cols <- c(
  "bandesh.primary" = "grey30",
  "fasolino.primary" = "grey50",
  "kang.primary" = "grey70",
  "akerman.D39-ZKSCAN1-WT37" = "mediumpurple4",
  "akerman.D45-H1" = "mediumpurple3",
  "akerman.D39-SIM1-WT8" = "mediumpurple2",
  "akerman.D39-SIM1-WT18" = "mediumpurple1",
  "balboa.H1.s7_w0" = "slateblue4",
  "balboa.H1.s7_w3" = "slateblue3",
  "balboa.H1.s7_w6" = "slateblue2",
  "kelley.HUES8.s6_w1" = "maroon",
  "wu.H1.s6_w0" = "hotpink2",
  "wu.HS980.s6_w0" = "hotpink1",
  "dadheech.iPSC.s6_w1" = "thistle"
)


# create a column for plotting best (cant contain NAs)
# CHANGE SO PRIMARY ALL BECOME ONE, LIKE IN VOLCANOS 
merged@meta.data <- merged@meta.data %>% 
  mutate(id.best = case_when(
    id %in% paste0(primary, ".primary") ~ "primary",
    id == "akerman.D39-ZKSCAN1-WT37" ~ "akerman.H1.rep1",
    id == "akerman.D39-SIM1-WT18" ~ "akerman.H3.rep2",
    id == "balboa.H1.s7_w3" ~ "balboa.H1.s7_w3",
    id == "balboa.H1.s7_w6" ~ "balboa.H1.s7_w6",
    id == "kelley.HUES8.s6_w1" ~ "kelley.HUES8.s6_w1",
    id == "wu.H1.s6_w0" ~ "wu.H1.s6_w0",
    TRUE ~ "other"
  ))

merged$id.best <- factor(
  merged$id.best,
  levels = c(
    "other", "primary",
    "akerman.H1.rep1", "akerman.H3.rep2",
    "balboa.H1.s7_w3", "balboa.H1.s7_w6",
    "kelley.HUES8.s6_w1", "wu.H1.s6_w0"
  )
)

# reverse id.best for plotting
merged$id.best.rev <- factor(merged$id.best, levels = rev(levels(merged$id.best)))

# set colours to use for id
id.best_cols <- c(
  "primary" = "grey50",
  "akerman.H1.rep1" = "mediumpurple4",
  "akerman.H3.rep2" = "mediumpurple1",
  "balboa.H1.s7_w3" = "slateblue3",
  "balboa.H1.s7_w6" = "slateblue2",
  "kelley.HUES8.s6_w1" = "maroon",
  "wu.H1.s6_w0" = "hotpink2"
)

# label all replicates appropriately
merged@meta.data <- merged@meta.data %>% 
  mutate(replicates = case_when(
    id == "bandesh.primary" ~ "primary.rep1",
    id == "fasolino.primary" ~ "primary.rep2",
    id == "kang.primary" ~ "primary.rep3",
    id == "akerman.D39-ZKSCAN1-WT37" ~ "akerman.H1.rep1",
    id == "akerman.D45-H1" ~ "akerman.H1.rep2",
    id == "akerman.D39-SIM1-WT8" ~ "akerman.H3.rep1",
    id == "akerman.D39-SIM1-WT18" ~ "akerman.H3.rep2",
    sample == "SRX10188829" ~ "balboa.s7_w0.rep1",
    sample == "SRX10188830" ~ "balboa.s7_w0.rep2",
    sample == "SRX10188831" ~ "balboa.s7_w0.rep3",
    sample == "SRX10188832" ~ "balboa.s7_w3.rep1",
    sample == "SRX10188833" ~ "balboa.s7_w3.rep2",
    sample == "SRX10188834" ~ "balboa.s7_w3.rep3",
    sample == "SRX10188835" ~ "balboa.s7_w6.rep1",
    sample == "SRX10188836" ~ "balboa.s7_w6.rep2",
    sample == "1" ~ "kelley.rep1",
    sample == "2" ~ "kelley.rep2",
    sample == "SRX25096888" ~ "wu.H1.rep1",
    sample == "SRX25096890" ~ "wu.HS980.rep1",
    sample == "SRX24523042" ~ "dadheech.rep1",
    sample == "SRX24523043" ~ "dadheech.rep2"
  ))
merged$replicates <- factor(
  merged$replicates,
  levels = c(
    "primary.rep1", "primary.rep2", "primary.rep3",
    "akerman.H1.rep1", "akerman.H1.rep2",
    "akerman.H3.rep1", "akerman.H3.rep2",
    "balboa.s7_w0.rep1", "balboa.s7_w0.rep2", "balboa.s7_w0.rep3",
    "balboa.s7_w3.rep1", "balboa.s7_w3.rep2", "balboa.s7_w3.rep3",
    "balboa.s7_w6.rep1", "balboa.s7_w6.rep2",
    "kelley.rep1", "kelley.rep2",
    "wu.H1.rep1", 
    "wu.HS980.rep1",
    "dadheech.rep1", "dadheech.rep2"
  )
)

# ---- final qc ----
merged$logGenePerUMI <- log1p(merged$nFeature_RNA) / log1p(merged$nCount_RNA)
merged$percent.rp <- PercentageFeatureSet(merged, pattern = "^RP[LS]")

qc_markers <- c("nFeature_RNA", "nCount_RNA", "pANN", "logGenePerUMI", "percent.mt", "percent.rp")

for (q in qc_markers) {
  print(VlnPlot(merged, group.by = "id", features = q, cols = id_cols, pt.size = 0) + NoLegend())
  ggsave(file.path(outdir, "VlnPlots", paste0(q, ".png")), height = 8, width = 8)
  remove(q)
}

# # remove remaining poor quality cells
# merged <- subset(merged, logGenePerUMI < 0.95)

# ---- Start analysis ----
merged <- SCTransform(merged, verbose = FALSE)
gc()

merged <- RunPCA(merged, verbose = FALSE)
print(ElbowPlot(merged, ndims = 30))
#ggsave(file.path(outdir, "merged_elbowplot.png"), width = 11, height = 10)

merged <- FindNeighbors(merged, dims = FindNeighbors.dims, verbose = FALSE)
merged <- FindClusters(merged, resolution = FindClusters.res, cluster.name = "unintegrated_clusters", verbose = FALSE)
merged <- RunUMAP(merged, dims = FindNeighbors.dims, reduction.name = "umap.unintegrated", verbose = FALSE)

# ---- Plotting ----
DimPlot(merged, reduction = "umap.unintegrated", raster = FALSE, group.by = "id", cols = id_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_merged_id.png"), width = 12, height = 10)

RidgePlot(merged, features = c("MALAT1", "percent.mt", "percent.rp"), group.by = "id", cols = id_cols) + NoLegend()

for (q in qc_markers) {
  print(FeaturePlot(merged, raster = FALSE, reduction = "umap.unintegrated", features = q, cols = plot_cols))
  ggsave(file.path(outdir, "FeaturePlots_unintegrated", paste0(q, ".png")), height = 8, width = 8)
  remove(q)
}

FeaturePlot(merged, raster = FALSE, reduction = "umap.unintegrated", features = "MALAT1", cols = plot_cols, min.cutoff = 3)

FeaturePlot(merged, reduction = "umap.unintegrated", raster = FALSE, features = "logGenePerUMI", min.cutoff = 0.95, cols = plot_cols)
ggsave(file.path(outdir, "FeaturePlots", "logGenePerUMI.min95.png"), width = 11, height = 10)

for (g in c("INS", "GCG", "TPH1", "KRT19")) {
  FeaturePlot(merged, reduction = "umap.unintegrated", raster = FALSE, features = g, cols = plot_cols)
  ggsave(file.path(outdir, "FeaturePlots_unintegrated", paste0(g, ".png")), width = 11, height = 10)
}

for (g in c("unintegrated_clusters", "orig.ident", "sample")) {
  print(DimPlot(merged, reduction = "umap.unintegrated", group.by = g, raster = FALSE))
  ggsave(file.path(outdir, "DimPlots", paste0("DimPlot_unintegrated_", g, ".png")), width = 11, height =10)
  remove(g)
}

# experiments dont mix well, integrate

# ---- integration ----
integrated <- merged; remove(merged); gc()

# ## ---- harmony ----
# integrated <- IntegrateLayers(
#   integrated,
#   method = HarmonyIntegration,
#   normalization.method = "SCT",
#   new.reduction = "harmony",
#   verbose = FALSE
# )
# 
# integrated <- integrated %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>% 
#   FindClusters(cluster.name = "harmony_clusters", resolution = 1.5, verbose = FALSE) %>% 
#   RunUMAP(reduction = "harmony", reduction.name = "umap.harmony", dims = 1:30, verbose = FALSE)
# 
# DimPlot(integrated, reduction = "umap.harmony", raster = FALSE, group.by = "id.akerman_samples")
# DimPlot(integrated, reduction = "umap.harmony", raster = FALSE, group.by = "harmony_clusters")
# FeaturePlot(integrated, reduction = "umap.harmony", raster = FALSE, features = c("INS", "GCG", "TPH1", "GHRL"))

## ---- scvi integration ----
# pull SCT HVGs before swapping to RNA assay
sct_hvgs <- VariableFeatures(integrated, nfeatures = 5000)

# uses raw counts (needs RNA assay)
# https://satijalab-seurat-wrappers.mintlify.app/methods/scvi
# https://cbib.github.io/Seurat-Integrate/reference/scVIIntegration.html
DefaultAssay(integrated) <- "RNA"
integrated <- integrated %>% 
  FindVariableFeatures(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(verbose = FALSE)

integrated <- IntegrateLayers(
  integrated,
  method = scVIIntegration,
  features = head(sct_hvgs, 3000), #SCT
  # features = VariableFeatures(integrated), # RNA
  ndims = 20,
  conda_env = "~/anaconda3/envs/scvi-env",
  new.reduction = "scvi",
  verbose = FALSE
)

# check all embeddings are useful
png(filename = file.path(outdir, "barplot_scvi-embeddings.png"))
barplot(apply(Embeddings(integrated[["scvi"]]), 2, var))
dev.off()

integrated <- integrated %>% 
  FindNeighbors(reduction = "scvi", dims = 1:20, verbose = FALSE) %>% 
  FindClusters(resolution = c(0.1, 2, 4), cluster.name = "scvi_clusters", verbose = FALSE) %>% 
  RunUMAP(reduction = "scvi", dims = 1:20, reduction.name = "umap.scvi", verbose = FALSE)
  
DimPlot(integrated, reduction = "umap.scvi", raster = FALSE, group.by = "scvi_clusters_2", label = TRUE)
FeaturePlot(integrated, reduction = "umap.scvi", raster = FALSE, features = c("INS", "GCG", "TPH1", "GHRL"))

## ---- fastmnn ----
integrated <- IntegrateLayers(
  integrated,
  method = FastMNNIntegration,
  new.reduction = "fastmnn",
  verbose = FALSE
)

integrated <- integrated %>% 
  FindNeighbors(reduction = "fastmnn", dims = 1:30, verbose = FALSE) %>% 
  FindClusters(cluster.name = "fastmnn_clusters", resolution = 1.5, verbose = FALSE) %>% 
  RunUMAP(reduction = "fastmnn", dims = 1:30, reduction.name = "umap.fastmnn", verbose = FALSE)

# ---- plot markers ----
gene_markers <- c(
  "ISL1", "NKX6-1", "PDX1",
  "INS", "IAPP", 
  "GCG", "ARX", 
  "SST", "GHRL", "PPY",
  "TPH1", "LMX1A", # EC
  "KRT19", # ductal
  "PRSS1", # acinar
  "RGS5", # pericytes
  "PECAM1", # endothelial
  "PTPRC", "C1QC", # all immune
  "COL1A1", "COL3A1", # stellate/stromal
  "FABP4", # quiescent stellate
  "TOP2A", # cycling
  "GATA4", "GATA6", "FOXA1", "ONECUT1", # progenitor
  "NGFR", "SLIT2", # schwann
  "TPSB2" # mast
)
markers <- c(qc_markers, gene_markers)

for (reduc in c("scvi", "fastmnn")) {
  
  reduc_name <- paste("umap", reduc, sep = ".")
  clust_name <- paste(reduc, "clusters", sep = "_")
  
  dir.create(file.path(outdir, paste0("FeaturePlots_", reduc)), showWarnings = FALSE)
  dir.create(file.path(outdir, paste0("VlnPlots_", reduc)), showWarnings = FALSE)
  
  # plot new data
  for (m in markers) {
    FeaturePlot(integrated, reduction = reduc_name, features = m, raster = FALSE, cols = plot_cols)
    ggsave(file.path(outdir, paste0("FeaturePlots_", reduc), paste0(m,".png")), width = 11, height = 10)
    VlnPlot(integrated, group.by = clust_name, features = m) + NoLegend()
    ggsave(file.path(outdir, paste0("VlnPlots_", reduc), paste0(m, ".png")), width = 11, height = 10)

    remove(m)
  }

  print(DimPlot(integrated, reduction = reduc_name, group.by = "id", raster = FALSE, cols = id_cols))
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_orig-ident.png", clust_name)), height = 10, width = 11)
  print(DimPlot(integrated, reduction = reduc_name, group.by = clust_name, raster = FALSE, label = TRUE) + NoLegend())
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_numbered.png", clust_name)), height = 10, width = 11)
  print(DimPlot(integrated, reduction = reduc_name, group.by = "id", split.by = "source", raster = FALSE, cols = id_cols))
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_orig-ident_source.png", clust_name)), height = 10, width = 22)
  
  remove(reduc, reduc_name, clust_name)
  
}

other_markers <- list(
  ribo = c("FBL", "NCL", "NPM1", "DKC1", "SBDS"),
  #ribo = c("RPL24", "RPL10", "EFL1", "NMD3", "RPS3", "LTV1"),
  er = c("XBP1", "HSPA5", "HSP90B1", "P4HB", "PDIA4"),
  mito = c("PPARGC1A", "NRF1", "TFAM", "COX4I1", "ATP5F1A"),
  fat = c("PPARA", "CPT1A", "ACOX1", "ACADM", "HADH"),
  gluc = c("SLC2A1", "PDK4", "PFKM", "PKM", "LDHA")
)

DotPlot(integrated, assay = "RNA", group.by = "id.rev", features = other_markers, cols = plot_cols) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(outdir, "DotPlots", "DotPlot_ribo-er-mito-fat-gluc_markers.png"))

# ---- find markers ----
integrated <- PrepSCTFindMarkers(integrated, verbose = FALSE)

integrated.markers <- FindAllMarkers(
  integrated,
  assay = "SCT",
  group.by = "scvi_clusters_2",
  verbose = FALSE
)

# ---- labelling clusters ----
integrated@meta.data <- integrated@meta.data %>%
  mutate(cell_type = case_when(
    
    # micro res
    scvi_clusters_4 %in% c(75) ~ "delta_GCG+",
    
    # normal resolution
    scvi_clusters_2 %in% c(68) ~ "schwann", # S100B, NGFR
    scvi_clusters_2 %in% c(38) ~ "delta_INS+",
    scvi_clusters_2 %in% c(10) ~ "prog.epithelial", # KRT19, CLDN4
    scvi_clusters_2 %in% c(40) ~ "stellate_q", # FABP4
    scvi_clusters_2 %in% c(64) ~ "epsilon", # GHRL
    scvi_clusters_2 %in% c(58) ~ "pericyte", # RGS5
    scvi_clusters_2 %in% c(52) ~ "msc", # COL3A1, PDGFRA
    scvi_clusters_2 %in% c(34,53) ~ "poly",
    scvi_clusters_2 %in% c(67) ~ "mast", # TPSB2
    scvi_clusters_2 %in% c(47) ~ "prog.neuroepi", # ONCECUT1, SOX9, HES4
    scvi_clusters_2 %in% c(60) ~ "neuronal", # SLIT2
    scvi_clusters_2 %in% c(26,61) ~ "prog.enteroendo", # PITX1, ONECUT2, ARX
    scvi_clusters_2 %in% c(55) ~ "prog.foregut", # YAP1, GATAs
    scvi_clusters_2 %in% c(29,59) ~ "prog.hepato", # TBC
    scvi_clusters_2 %in% c(0,6,4,7,35,50,49,36,27,12,20,15,17,13,28) ~ "beta",  # INS
    scvi_clusters_2 %in% c(48,3,37,25,1,11,16,41,32,8,2,30) ~ "alpha", # GCG
    scvi_clusters_2 %in% c(46,14,42) ~ "delta", # SST
    scvi_clusters_2 %in% c(33,65,62) ~ "ppy", # PPY
    scvi_clusters_2 %in% c(5,24,31,71,21) ~ "ec", # TPH1
    scvi_clusters_2 %in% c(9,54,39,45,19,44,56,69) ~ "ductal", # KRT19
    scvi_clusters_2 %in% c(18) ~ "stellate", # RGS5, COL1A1
    scvi_clusters_2 %in% c(22,70) ~ "endothelial", # PECAM1
    scvi_clusters_2 %in% c(51) ~ "immune", # PTPRC
    scvi_clusters_2 %in% c(57,43,23,66,63) ~ "acinar", # PRSS1
    
    TRUE ~ scvi_clusters_2
  ))
integrated$cell_type <- factor(
  integrated$cell_type,
  levels = c(
    "beta", "alpha", "ec",
    "delta", "delta_INS+", "delta_GCG+", "ppy", "epsilon", "poly",
    "acinar", "ductal",
    "stellate", "stellate_q", "pericyte", "endothelial",
    "immune", "mast", "schwann", "neuronal",
    "msc", 
    "prog.epithelial", "prog.enteroendo", "prog.foregut", "prog.hepato", "prog.neuroepi"
  )
)
integrated$cell_type.rev <- factor(
  integrated$cell_type,
  levels = rev(levels(integrated$cell_type))
)

cell_type_cols <- DiscretePalette(length(levels(integrated$cell_type)))
names(cell_type_cols) <- levels(integrated$cell_type)

DimPlot(integrated, reduction = "umap.scvi", group.by = "cell_type", label = TRUE, raster = FALSE, cols = cell_type_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_integrated_cell-type.png"), width = 12, height = 10)

DimPlot(integrated, reduction = "umap.scvi", group.by = "cell_type", raster = FALSE, cols = cell_type_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_integrated_cell-type.unlabelled.png"), width = 12, height = 10)

DotPlot(integrated, features = gene_markers, group.by = "cell_type", cols = plot_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "DotPlots", "DotPlot_integrated_cell-type.png"), width = 12, height = 5) 

# label source
DimPlot(integrated, reduction = "umap.scvi", raster = FALSE, group.by = "cell_type", split.by = "source")
ggsave(file.path(outdir, "DimPlots", "DimPlot_integrated_cell-type_primary-stemcell.png"), width = 15, height = 8)

# create UMAPs for each lab
Idents(integrated) <- "id"
groups <- list(
  "primary" = grep("primary", unique(integrated$id), value = TRUE), 
  "akerman" = grep("akerman", unique(integrated$id), value = TRUE), 
  "balboa" = grep("balboa", unique(integrated$id), value = TRUE), 
  "dadheech" = grep("dadheech", unique(integrated$id), value = TRUE), 
  "wu" = grep("wu", unique(integrated$id), value = TRUE), 
  "kelley" = grep("kelley", unique(integrated$id), value = TRUE)
)

# loop through each lab
for (i in seq_along(groups)) {
  
  g <- groups[[i]]
  n <- names(groups)[[i]]

  DimPlot(
    integrated, 
    reduction = "umap.scvi", 
    raster = FALSE, 
    cells.highlight = WhichCells(integrated, idents = g), 
    cols.highlight = unname(id_cols[g]),
    cols = "grey90",
    sizes.highlight = 0.1
  )
  
  ggsave(file.path(outdir, "DimPlots", paste0("DimPlots_scvi_", n, ".png")), width = 11, height = 10)
  
  remove(i, g, n)
}

# plot percentages of cell types per id
ggplot(integrated@meta.data, aes(fill=cell_type.rev, x=id.rev)) + 
  geom_bar(position="fill", color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = cell_type_cols) +
  theme_classic()
ggsave(file.path(outdir, "geom-bar_stacked_replicate_cell-type.png"), width = 8, height = 6)

percent_celltypes <- integrated@meta.data %>%
  dplyr::count(id, cell_type) %>%          # number of cells per type per sample
  group_by(id) %>%
  mutate(per = n / sum(n) * 100) %>%
  ungroup()

write.table(
  percent_celltypes,
  file.path(outdir, "percent_celltypes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggplot(percent_celltypes, aes(fill = cell_type, y = per, x = id)) +
  geom_bar(position = "dodge", stat = "identity", colour = "grey30") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = percent_celltypes %>% filter(cell_type == "beta"),
    aes(label = round(per, 1)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    hjust = 1
  )
ggsave(file.path(outdir, "geom-bar_id-celltype_grouped.labelled.png"), width = 11, height = 6)

# ---- investigate EC cells ----
ec_markers <- c(
  "TPH1", "DDC", # serotonin enzymes
  "SLC18A1", # serotonin packing
  "AANAT", "ASMT", # make melatonin
  "HTR2B", "HTR1F", "HTR2C", "HTR3A", # serotonin receptors
  "MTNR1A", "MTNR1B", # melatonin receptors
  "PIEZO2" # pressure > action potential
)
ec_markers.clean <- setdiff(ec_markers, rownames(integrated[["RNA"]]))
ec_markers.clean <- setdiff(ec_markers, ec_markers.clean)

for (e in ec_markers.clean) {
  print(FeaturePlot(integrated, features = e, reduction = "umap.scvi", raster = FALSE, cols = plot_cols))
  ggsave(file.path(outdir, "FeaturePlots_scvi", paste0(e, ".png")), width = 11, height = 10)
  print(VlnPlot(integrated, features = e, assay = "RNA", group.by = "cell_type", split.by = "source", cols = pri_vs_sc))
  ggsave(file.path(outdir, "VlnPlots", paste0("VlnPlot_", e, "_celltype-source.png")), width = 15, height = 8)
  remove(e)
}

Idents(integrated) <- "source"

for (s in c("primary", "stem-cell")) {
  
  cols <- if (s == "primary") {
    c("lightgrey", "grey20")
  } else {
    c("lightgrey", "purple")
  } 
  
  DotPlot(
    integrated, 
    features = ec_markers.clean, 
    group.by = "cell_type", 
    idents = s,
    cols = cols
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(outdir, "DotPlots", paste0("DotPlot_integrated_ildem-ec-markers.", s, ".png")), width = 8, height = 10)

  remove(s, cols)
}
  
insulin_regulators <- c(
  "GIP",        # glucose-dependent insulinotropic polypeptide
  "GH1",        # growth hormone
  "CRH",        # corticotropin-releasing hormone
  "POMC",       # ACTH precursor
  "NR3C1",      # glucocorticoid receptor (cortisol signaling)
  "ADCYAP1",    # PACAP
  "VIP"         # vasoactive intestinal peptide
)

for (s in c("primary", "stem-cell")) {
  
  cols <- if (s == "primary") {
    c("lightgrey", "grey20")
  } else {
    c("lightgrey", "purple")
  } 

  DotPlot(
    integrated, 
    features = insulin_regulators, 
    group.by = "cell_type", 
    idents = s,
    cols = cols
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(outdir, "DotPlots", paste0("DotPlot_integrated_chatgpt-beta-hormones.", s, ".png")), width = 8, height = 10)

  remove(s, cols)
}
  
# now just plot EC cells
ec <- subset(
  integrated,
  cell_type %in% c("ec", "prog.epithelial") &
    source == "stem-cell"
)

DefaultAssay(ec) <- "RNA"
ec[["SCT"]] <- NULL
gc()

DotPlot(
  ec, 
  features = ec_markers.clean, 
  group.by = "id.rev",
  cols = c("lightgrey", "purple")
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "EC cell RNA expr - functionality markers")
ggsave(file.path(outdir, "DotPlots", "DotPlot_ec_ildem-ec-markers.png"), width = 8, height = 10)

DotPlot(
  ec, 
  features = insulin_regulators, 
  group.by = "id.rev",
  cols = c("lightgrey", "purple")
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "EC cell RNA expr - negulators of beta cells")
ggsave(file.path(outdir, "DotPlots", "DotPlot_ec_chatgpt-beta-hormones.png"), width = 8, height = 10)

# ---- beta cell specific analysis ----
# subset beta cells first to remove NA
beta <- subset(integrated, cell_type == "beta")
DefaultAssay(beta) <- "RNA"
beta[["SCT"]] <- NULL
gc()

# pipeline and reanalyse just beta cells
beta[["RNA"]] <- JoinLayers(beta[["RNA"]])
beta[["RNA"]] <- split(beta[["RNA"]], f = beta$sample)
beta <- SCTransform(beta, verbose = FALSE)
gc()

# cluster beta cells separately
beta <- RunPCA(beta, verbose = FALSE)
ElbowPlot(beta)

beta <- FindNeighbors(beta, dims = 1:15, verbose = FALSE)
beta <- FindClusters(beta, resolution = 0.8, cluster.name = "unintegrated_clusters", verbose = FALSE)
beta <- RunUMAP(beta, dims = 1:15, reduction.name = "umap.unintegrated", verbose = FALSE)

# plot new clusters
DimPlot(beta, reduction = "umap.unintegrated", group.by = "id", cols = id_cols, raster = FALSE)
ggsave(file.path(outdir, "DimPlots", paste0("DimPlot_beta_id.png")), width = 11, height = 10)

beta$split.by <- ifelse(
  beta$orig.ident %in% primary,
  "primary",
  as.character(beta$orig.ident)
)
beta$split.by <- factor(beta$split.by, levels = c("primary", "akerman", "balboa", "kelley", "wu", "dadheech"))

DimPlot(beta, reduction = "umap.unintegrated", group.by = "id", cols = id_cols, split.by = "split.by", raster = FALSE)
ggsave(file.path(outdir, "DimPlots", "DimPlot_beta_id-splitby.png"), width = 30, height = 8)

'---
rpca and harmony integration overcorrects beta cell clustering
only makes one large ball, we know they are all beta cells but different
i would argue integration is not needed now
integration probably wont work if there is just one cell type?
 
---'

# return to RNA assay for plotting
DefaultAssay(beta) <- "RNA"

for (q in qc_markers) {
  print(FeaturePlot(beta, reduction = "umap.unintegrated", raster = FALSE, features = q, cols = plot_cols))
  ggsave(file.path(outdir, "FeaturePlots", paste0("DimPlot_beta_", q, ".png")), width = 11, height = 10)
  
  print(VlnPlot(beta, group.by = "id", features = q, cols = id_cols, pt.size = 0) + NoLegend())
  ggsave(file.path(outdir, "VlnPlots", paste0("VlnPlot_beta_", q, ".png")), width = 8, height = 8)
  
  remove(q)
}

VlnPlot(
  beta, 
  layer = "data", 
  group.by = "replicates", 
  features = "INS",
  pt.size = 0,
  #cols = id_cols
) + NoLegend()
ggsave(file.path(outdir, "VlnPlots", "VlnPlot_beta-INS_id.png"), width = 10, height = 6)

# Actual differences between INS expr across conditions
avg.ins.expr <- AverageExpression(
  beta, 
  group.by = "id", 
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
  pivot_longer(cols = -gene, names_to = "source", values_to = "avg.logINS") %>% 
  select(-(gene)) %>% 
  mutate(avg.INS = expm1(avg.logINS))

write.table(
  avg.ins.expr_df,
  file.path(outdir, "avg-ins-expr_cell-line.tsv"),
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

DotPlot(beta, features = beta_markers, group.by = "id.best.rev", cols = plot_cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(
  file.path(outdir, "DotPlots", "DotPlot_beta-markers_id.best.rev.png"),
  width = 14, height = 6
)

VlnPlot(beta, features = beta_markers, group.by = "id.best.rev", cols = id.best_cols, stack = TRUE, fill.by = "ident") +
  NoLegend() +
  scale_fill_manual(values = id.best_cols)
ggsave(
  file.path(outdir, "VlnPlots", "VlnPlot_beta-markers_id.best.rev.png"),
  width = 14, height = 6
)

tremmel_markers <- list(
  cir = c("NR1D1", "ARNTL"),
  difr = c("MAFA", "SIX2", "SIX3", "ONECUT2"),
  tfs = c("KLF9", "HDAC9", "HOPX"),
  er = c("ERO1B"),
  gluc = c("SLC2A1", "SLC2A2", "GCK", "G6PC2"),
  rec = c("ITGA1", "SHISAL2B", "ENTPD3"),
  exo = c("IAPP", "SLC30A8", "CHGB", "SYT4")
)

DotPlot(beta, assay = "RNA", group.by = "id.rev", features = tremmel_markers, cols = plot_cols) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(outdir, "DotPlots", "DotPlot_tremmel_markers.png"), height = 6, width = 10)

# ---- beta dgea ----
beta <- PrepSCTFindMarkers(beta, verbose = FALSE)

# compare each groups beta cells to primary beta cells
beta.dgea.res <- list()
p <- paste0(primary, ".primary")

# loop through each group 
for (g in unique(beta$id)) {
  
  if (g %in% p) next # dont include primary
  
  dgea <- FindMarkers( # perform dgea
    beta,
    assay = "SCT",
    group.by = "id",
    ident.1 = g,
    ident.2 = p,
    min.pct = 0.25,
    verbose = FALSE
  )
  
  dgea$gene <- rownames(dgea) # add gene col
  
  # add grouping for counting and volcanos
  dgea <- dgea %>% 
    mutate(neg.log10.p_val_adj = -log10(p_val_adj + 1e-300)) %>% 
    mutate(expr = case_when(
      avg_log2FC > 1 & 
        neg.log10.p_val_adj > -log10(0.05 + 1e-300) ~ "sig_up",
      avg_log2FC < -1  &
        neg.log10.p_val_adj > -log10(0.05 + 1e-300) ~ "sig_down",
      TRUE ~ "no_diff"
    ))
  
  beta.dgea.res[[g]] <- dgea # save results
  remove(g, dgea) # clean up environment
}

saveRDS(beta.dgea.res, file.path(outdir, "FindMarkers_list_sbc-vs-pbc_per-id.RDS")) # save results
remove(p) # clean up

# count how many degs per group
beta.dgea.res.counts <- lapply(beta.dgea.res, function(df) {
  
  df <- df %>% 
    group_by(expr) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% 
    mutate(percent = count / sum(count) * 100)
})

# combine counts tables into one table and save
write.table(
  bind_rows(beta.dgea.res.counts, .id = "id"),
  file.path(outdir, "FindMarkers_dgea_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# plot dgea
for (i in seq_along(beta.dgea.res)) {
  
  df <- beta.dgea.res[[i]]
  nm <- names(beta.dgea.res[i])
  
  print(
    ggplot(df, aes(x = avg_log2FC, y = neg.log10.p_val_adj, colour = expr)) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(values = c(
        "sig_up" = "tomato",
        "sig_down" = "turquoise4",
        "no_diff" = "grey50"
      )) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) +
      geom_vline(xintercept = 1, linetype = 2) +
      geom_vline(xintercept = -1, linetype = 2) +
      geom_label_repel(
        data = subset(df, gene %in% c(beta_markers, "INS")),
        aes(label = gene),
        size = 6,
        max.overlaps = Inf
      ) +
      labs(title = nm) +
      theme_classic() +
      theme(legend.position = "none")
  )
  
  ggsave(file.path(outdir, "volcanos", paste0("sbc-vs-pbc.", nm, ".png")), width = 6, height = 6)
  
  remove(i, df, nm)
}

# plot corresponding percentages
for (i in seq_along(beta.dgea.res.counts)) {
  
  df <- beta.dgea.res.counts[[i]]
  nm <- names(beta.dgea.res.counts)[[i]]
  
  # match to order of volcanos
  df$expr <- factor(df$expr, levels = c("sig_down", "no_diff", "sig_up"))
  
  print(
    ggplot(df, aes(x = expr, y = percent, fill = expr)) +
      geom_col() +
      scale_fill_manual(values = c(
        "sig_down" = "turquoise4",
        "no_diff" = "grey50",
        "sig_up" = "tomato"
      )) +
      geom_text(
        aes(y = percent, label = paste0(round(percent,1), "%")),
        position = position_dodge(width = 0.9),
        vjust = -0.5,
        size = 8
      ) +
      labs(title = nm) +
      theme_classic() +
      theme(legend.position = "none") +
      scale_y_continuous(limits = c(0, 60))
  ) 
  
  ggsave(file.path(outdir, "geom-bar", paste0("geom-bar_dgea-counts_", nm, ".png")), height = 6, width = 6)
  
  remove(i, df, nm)
}

# ---- pseudobulk for cross comparison ----
integrated$cell_type.id <- ifelse(
  integrated$cell_type %in% c("beta", "alpha", "ec"),
  paste(integrated$cell_type, integrated$id, sep = "."),
  NA
)

# plot cell counts
ggplot(integrated@meta.data, aes(x = factor(
    cell_type.id, 
    levels = sort(unique(integrated$cell_type.id), decreasing = TRUE)
  ), 
  fill = cell_type.id)
) +
  geom_bar(show.legend = FALSE) +
  coord_flip() +
  theme_classic() + 
  expand_limits(y = max(table(integrated$cell_type.id)) * 1.2) +
  geom_text(
    stat = "count", 
    aes(label = after_stat(count)), 
    hjust = -0.2
  )
ggsave(file.path(outdir, "geom-bar_integrated_celltype-id.png"))

# plot percentages of total population
percent_celltypes.filt <- percent_celltypes %>% 
  filter(cell_type %in% c("beta", "alpha", "ec", "prog.epithelial")) %>% 
  mutate(cell_type = as.character(cell_type)) %>% 
  complete(id, cell_type, fill = list(per = 0))

ggplot(percent_celltypes.filt, aes(fill = id, y = per, x = cell_type)) +
  geom_bar(position = "dodge", stat = "identity", colour = "grey30") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = id_cols)
ggsave(file.path(outdir, "geom-bar_celltype-id_percentages.png"))

# identify groups with more than 50 cells in them
keep <- table(integrated$cell_type.id) > 50
keep <- names(keep)[keep]

# analyse rna instead, see differences with sct
expr_sets.agg_rna <- AggregateExpression(
  integrated,
  assays = "RNA",
  layers = "counts",
  group.by = "cell_type.id",
  features = head(sct_hvgs, 3000),
  normalization.method = NULL,
  scale.factor = NULL,
  margin = NULL,
  verbose = FALSE
)$RNA

# aggregateexpression changes "_" to "-", match keep to aggregateexpression
keep <- gsub("_", "-", keep)

# remove columns that have less than one cell
expr_sets.agg_rna <- expr_sets.agg_rna[, colnames(expr_sets.agg_rna) %in% keep]

# create deseq2.info
deseq2.info <- data.frame(
  sample = colnames(expr_sets.agg_rna),
  cell_type = sapply(strsplit(colnames(expr_sets.agg_rna), "\\."), `[`, 1),
  lab = sapply(strsplit(colnames(expr_sets.agg_rna), "\\."), `[`, 2),
  source = ifelse(
    sapply(strsplit(colnames(expr_sets.agg_rna), "\\."), `[`, 2) %in% primary,
    "primary", "stem-cell"
  ),
  id = sub("^[^.]*\\.", "", colnames(expr_sets.agg_rna))
)
rownames(deseq2.info) <- deseq2.info$sample

# create deseq object
dds <- DESeqDataSetFromMatrix(
  countData = round(expr_sets.agg_rna),
  colData = deseq2.info,
  design = ~1
)

# variance stabilise
vsd <- vst(dds, blind = TRUE)

## ---- euclidean distance ----
euc_dist <- dist(t(assay(vsd)))

# organise output
euc_dist <- as.matrix(euc_dist)
rownames(euc_dist) <- colnames(expr_sets.agg_rna)
colnames(euc_dist) <- colnames(expr_sets.agg_rna)

# plot
pheatmap::pheatmap(
  euc_dist,
  color = hcl.colors(100),
  main = "Euclidean distance",
  filename = file.path(outdir, "pheatmaps", "euc-dist_top3000hvg_alpha-beta-ec.png")
)

# hierarchical clustering
ed_hc <- hclust(as.dist(euc_dist), method = "average")

png(file.path(outdir, "dendrograms", "euc-dist_top3000hvg_alpha-beta-ec.png"))
plot(
  as.dendrogram(ed_hc), 
  horiz = TRUE, 
  xlim = c(max(euc_dist), -100),
  xlab = "Euclidean distance"
)
dev.off()

## ---- pearson ----
ps_cor <- cor(as.matrix(expr_sets.agg_rna), method = "pearson")^2

# add colouring
ps_cor_annot <- data.frame(
  cell_type = sapply(strsplit(colnames(ps_cor), "\\."), '[', 1)
  #id = sub("^([^\\.]*\\.)", "", colnames(ps_cor))
)

rownames(ps_cor_annot) <- colnames(ps_cor)

annot_cols <- list(
  cell_type = c(
    cell_type_cols["beta"],
    cell_type_cols["alpha"],
    cell_type_cols["ec"]
  )
)

# without values
pheatmap::pheatmap(
  ps_cor, 
  color = hcl.colors(100, rev = TRUE),
  main = "Pearson corr (r^2)",
  annotation_col = ps_cor_annot,
  annotation_row = ps_cor_annot,
  annotation_colors = annot_cols,
  show_colnames = FALSE,
  cutree_rows = 3,
  cutree_cols = 3,
  filename = file.path(outdir, "pheatmaps", "ps-cor_top3000hvg_alpha-beta-ec.png"),
  width = 9.5
); dev.off()

# with values
pheatmap::pheatmap(
  ps_cor, 
  color = hcl.colors(100, rev = TRUE),
  main = "Pearson corr (r^2)",
  cutree_rows = 3,
  cutree_cols = 3,
  display_numbers = TRUE,
  filename = file.path(outdir, "pheatmaps", "ps-cor_top3000hvg_alpha-beta-ec.values.png"),
  width = 16, height = 16
); dev.off()

# hierarchical clustering
ps_dist <- as.dist(1 - ps_cor)
ps_hc <- hclust(ps_dist, method = "average")

png(file.path(outdir, "dendrograms", "ps-cor_top3000hvg_alpha-beta-ec.png"))
plot(
  as.dendrogram(ps_hc), 
  horiz = TRUE, 
  xlim = c(1, -0.4), 
  xlab = "Pearson corr (r^2) distance"
)
dev.off()

## ---- spearman ----
sm_cor <- cor(as.matrix(expr_sets.agg_rna), method = "spearman")^2
pheatmap::pheatmap(
  sm_cor, 
  color = hcl.colors(100, rev = TRUE),
  main = "Spearman corr (r^2)",
  filename = file.path(outdir, "pheatmaps", "sm-cor_top3000hvg_alpha-beta-ec.png")
); dev.off()

# hierarchical clustering
sm_dist <- as.dist(1 - sm_cor)
sm_hc <- hclust(sm_dist, method = "average")

png(file.path(outdir, "dendrograms", "sm-cor_top3000hvg_alpha-beta-ec.png"))
plot(
  as.dendrogram(sm_hc), 
  horiz = TRUE, 
  xlim = c(0.5, -0.2),
  xlab = "Spearman corr (r^2) distance"
)
dev.off()

## ---- PCA ----
pca_data <- plotPCA(
  vsd, 
  ntop = 3000,
  intgroup = c("cell_type", "id"), 
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

# fix id_cols to match new names without _
id_cols.clean <- id_cols
names(id_cols.clean) <- gsub("_", "-", names(id_cols))

# fix order of id in pca_data
pca_data$id <- factor(pca_data$id, levels = names(id_cols.clean))

ggplot(pca_data, aes(PC1, PC2, color = id, shape = cell_type)) +
  geom_point(size = 3) +
  scale_color_manual(values = id_cols.clean) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_classic()

ggsave(file.path(outdir, "pca", "pca_top3000hvg_alpha-beta-ec.png"), width = 9, height = 9)

# # ---- Deseq2 ----
# dds <- DESeq(dds)
# 
# # create ranked list, ranked with p values and signed with fc
# res <- results(
#   dds,
#   contrast = c("condition", "sc-beta", "primary-beta")
# )
# 
# res_df <- as.data.frame(res)
# res_df$gene <- rownames(res_df)
# 
# res_df %>%
#   select(gene, log2FoldChange, padj) %>%
#   filter(!is.na(padj)) %>%
#   mutate(
#     neglog10p = -log10(pmax(padj, 1e-300)),
#     neglog10p.sign = neglog10p * ifelse(log2FoldChange >= 0, 1, -1)
#   ) %>%
#   select(gene, neglog10p.sign) %>%
#   arrange(desc(neglog10p.sign)) %>%
#   write.table(
#     file.path(outdir, "DESeq2_SBC-vs-PBC_neg-log10-sign.rnk"),
#     row.names = FALSE,
#     col.names = FALSE,
#     quote = FALSE,
#     sep = "\t"
#   )

# --- Lickert contaminants ----
lickert_contaminants <- fread(
  file.path(indir, "../..", "raw_data", "Lickert_contaminants.csv"),
  select = c("gene", "category"),
  verbose = FALSE
)

# join layers for GetAssayData
beta[["RNA"]] <- JoinLayers(beta[["RNA"]])

# for SBCs and PBCs, calculate pct expressed for each contaminant gene
for (bc in unique(beta$id)) {
  
  beta.bc <- subset(beta, id == bc)
  genes <- intersect(lickert_contaminants$gene, rownames(beta.bc))
  lc_expr <- GetAssayData(beta.bc[genes, ], assay = "RNA", layer = "data")
  
  # set pct thresholds
  for (pct in c(0, 0.5, 1, 2)) {
    
    lc_expr.pct <- rowMeans(lc_expr > pct) * 100
    col_name <- paste0(bc, "_", pct, "fc")
    lickert_contaminants[[col_name]] <- lc_expr.pct[lickert_contaminants$gene]
    
    remove(pct, lc_expr.pct, col_name)
  }
  remove(bc, beta.bc, genes, lc_expr)
}

# order the contaminants by category and columns alphabetically
lickert_contaminants <- as.data.frame(lickert_contaminants)
lickert_contaminants <- lickert_contaminants[order(lickert_contaminants$category),]
lickert_contaminants <- lickert_contaminants[,sort(names(lickert_contaminants))]
rownames(lickert_contaminants) <- NULL # reset rownames

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
  select(ends_with("fc")) |>
  as.matrix()

annotation_row <- data.frame(Category = lickert_contaminants$category)
rownames(annotation_row) <- rownames(mat)

annotation_col <- data.frame(beta_cell = sub("_\\d.*$", "", colnames(mat)))
rownames(annotation_col) <- colnames(mat)

pheatmap::pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  display_numbers = FALSE,
  #fontsize = 15,
  filename = file.path(outdir, "pheatmaps", "pheatmap_pct.expr.lickert_contaminants.png"),
  width = 12,
  height = 12,
  color = hcl.colors(100),
); dev.off()

remove(annotation_col, annotation_row)

# ---- metabolism ----
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

metab_pathways <- list(
  REACTOME_MITOCHONDRIAL_BIOGENESIS = file.path(path, "raw_data", "REACTOME_MITOCHONDRIAL_BIOGENESIS.v2026.1.Hs.grp"),
  GOBP_MITOCHONDRION_ORGANIZATION = file.path(path, "raw_data", "GOBP_MITOCHONDRION_ORGANIZATION.v2026.1.Hs.grp"),
  GOBP_RIBOSOME_BIOGENESIS = file.path(path, "raw_data", "GOBP_RIBOSOME_BIOGENESIS.v2026.1.Hs.grp"),
  GOBP_RRNA_PROCESSING = file.path(path, "raw_data", "GOBP_RRNA_PROCESSING.v2026.1.Hs.grp")
)

metab_pathways <- lapply(metab_pathways, scan, what = character(), skip = 2, quiet = TRUE)

metab_pathways <- lapply(metab_pathways, function(x) {
  AggregateExpression(
    integrated,
    assay = "RNA",
    features = intersect(x, rownames(integrated)),
    group.by = "id"
  )$RNA
})

lapply(seq_along(metab_pathways), function(i) {
  data <- metab_pathways[[i]]
  nm <- names(metab_pathways)[[i]]
  
  pheatmap::pheatmap(
    t(as.matrix(data)),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    color = colorRampPalette(plot_cols)(100),
    border_color = NA,
    gaps_row = c(3,7,10,11,13),
    scale = "column",
    main = nm,
    filename = file.path(outdir, "pheatmaps", paste0("pheatmap_", nm, "_integrated.png"))
  )
})


# ---- save output ----
saveRDS(integrated, file.path(outdir, "integrated.v4.RDS"))