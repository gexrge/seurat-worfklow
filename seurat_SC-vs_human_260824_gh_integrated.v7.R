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
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridis)

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
outdir <- file.path(path, "processed/seurat_SC-vs-Human_260824_gh_integrated.v7")

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
data.dirs <- data.dirs[grepl("(akerman|bandesh|fasolino|kang)", data.dirs)]
rds.dirs <- list.files(data.dirs, pattern = ".rds", full.names = TRUE)
data.list <- lapply(rds.dirs, readRDS)

# label data
name.dirs <- sapply(data.dirs, basename)
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
    "akerman"
  ) 
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
    grepl("^(D39-Z|D45)", sample) ~ "H1",
    grepl("^D39-S", sample) ~ "H3",
    TRUE ~ "primary"
  ))
merged$cell_line <- factor(
  merged$cell_line,
  levels = c("primary", "H1", "H3")
)

# make an id that groups samples appropriately
merged$id <- ifelse(
  merged$orig.ident %in% primary,
  paste(merged$orig.ident, "primary", sep = "."),
  paste(merged$orig.ident, merged$sample, sep = ".")
)
merged$id <- factor(
  merged$id,
  levels = c(
    "bandesh.primary", "fasolino.primary", "kang.primary",
    "akerman.D39-ZKSCAN1-WT37", "akerman.D45-H1",
    "akerman.D39-SIM1-WT8", "akerman.D39-SIM1-WT18"
  )
)

# reverse order for plotting
merged$id.rev <- factor(merged$id, levels = rev(levels(merged$id)))

id_cols <- c(
  "bandesh.primary" = "grey30",
  "fasolino.primary" = "grey50",
  "kang.primary" = "grey70",
  "akerman.D39-ZKSCAN1-WT37" = "mediumpurple4",
  "akerman.D45-H1" = "mediumpurple3",
  "akerman.D39-SIM1-WT8" = "mediumpurple2",
  "akerman.D39-SIM1-WT18" = "mediumpurple1"
)

# ---- final qc ----
merged$logGenePerUMI <- log1p(merged$nFeature_RNA) / log1p(merged$nCount_RNA)
merged$percent.rp <- PercentageFeatureSet(merged, pattern = "^RP[LS]")

qc_markers <- c("nFeature_RNA", "nCount_RNA", "pANN", "logGenePerUMI", "percent.mt", "percent.rp")

for (q in qc_markers) {
  print(VlnPlot(
    integrated, 
    group.by = "id", 
    features = q,
    cols = id_cols
  ) + NoLegend())
  ggsave(file.path(outdir, "VlnPlots", paste0(q, ".png")), height = 8, width = 8)
  remove(q)
}

# remove remaining poor quality cells
merged <- subset(merged, logGenePerUMI < 0.95)

# ---- Start analysis ----
merged <- SCTransform(merged, verbose = FALSE)
gc()

merged <- RunPCA(merged, verbose = FALSE)
print(ElbowPlot(merged, ndims = 30))
ggsave(file.path(outdir, "merged_elbowplot.png"), width = 11, height = 10)

merged <- FindNeighbors(merged, dims = FindNeighbors.dims, verbose = FALSE)
merged <- FindClusters(merged, resolution = FindClusters.res, cluster.name = "unintegrated_clusters", verbose = FALSE)
merged <- RunUMAP(merged, dims = FindNeighbors.dims, reduction.name = "umap.unintegrated", verbose = FALSE)

# ---- Plotting ----
for (g in c("unintegrated_clusters", "orig.ident", "sample")) {
  print(DimPlot(merged, reduction = "umap.unintegrated", group.by = g, raster = FALSE))
  ggsave(file.path(outdir, "DimPlots", paste0("DimPlot_unintegrated_", g, ".png")), width = 11, height =10)
  remove(g)
}

DimPlot(merged, reduction = "umap.unintegrated", raster = FALSE, group.by = "id", cols = id_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_unintegrated_id.png"), width = 11, height = 10)

# experiments dont mix well, integrate

# ---- integration ----
integrated <- merged; remove(merged); gc()

## ---- scvi ----
# pull SCT HVGs before swapping to RNA assay
sct_hvgs <- VariableFeatures(integrated, nfeatures = 5000)

# scvi integration - uses raw counts (needs RNA assay)
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
  features = head(sct_hvgs, 3000), # SCT
  # features = VariableFeatures(integrated),
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
  FindClusters(resolution = 2, cluster.name = "scvi_clusters_2", verbose = FALSE) %>% 
  FindClusters(resolution = 0.4, cluster.name = "scvi_clusters", verbose = FALSE) %>% 
  RunUMAP(reduction = "scvi", dims = 1:20, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(integrated, reduction = "umap.scvi", raster = FALSE, group.by = "id")
DimPlot(integrated, reduction = "umap.scvi", raster = FALSE, group.by = "scvi_clusters", label = TRUE)
FeaturePlot(integrated, reduction = "umap.scvi", raster = FALSE, features = c("INS", "GCG", "TPH1", "GHRL"))

## ---- fastmnn ----
integrated <- IntegrateLayers(
  integrated,
  method = FastMNNIntegration,
  new.reduction = "fastmnn",
  verbose = FALSE
)
integrated <- FindNeighbors(integrated, reduction = "fastmnn", dims = 1:20, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 1, cluster.name = "fastmnn_clusters", verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "fastmnn", dims = 1:20, reduction.name = "umap.fastmnn", verbose = FALSE)

DimPlot(integrated, reduction = "umap.fastmnn", raster = FALSE, group.by = "id")
FeaturePlot(integrated, reduction = "umap.fastmnn", raster = FALSE, features = c("INS", "GCG", "TPH1", "GHRL"))

## ---- rpca ----
integrated <- IntegrateLayers(
  integrated,
  method = RPCAIntegration,
  new.reduction = "rpca",
  verbose = FALSE
)

integrated <- integrated %>% 
  FindNeighbors(reduction = "rpca", dims = 1:20, verbose = FALSE) %>% 
  FindClusters(resolution = 1, cluster.name = "rpca_clusters", verbose = FALSE) %>% 
  RunUMAP(reduction = "rpca", dims = 1:20, reduction.name = "umap.rpca", verbose = FALSE)

DimPlot(integrated, reduction = "umap.rpca", raster = FALSE, group.by = "id")
FeaturePlot(integrated, reduction = "umap.rpca", raster = FALSE, features = c("INS", "GCG", "TPH1", "GHRL"))

# ---- plot markers ----
# return assay to SCT after integrating with RNA
DefaultAssay(integrated) <- "SCT"

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
  "NGFR", "SLIT2", # schwann
  "TPSB2" # mast
)
markers <- c(qc_markers, gene_markers)

for (reduc in c("scvi", "rpca", "fastmnn")) {
  
  reduc_name <- paste("umap", reduc, sep = ".")
  clust_name <- paste(reduc, "clusters", sep = "_")
  
  dir.create(file.path(outdir, paste0("FeaturePlots_", reduc)), showWarnings = FALSE)

  # plot new data
  for (m in markers) {
    FeaturePlot(integrated, reduction = reduc_name, features = m, raster = FALSE, cols = plot_cols)
    ggsave(file.path(outdir, paste0("FeaturePlots_", reduc), paste0(m,".png")), width = 11, height = 10)

    remove(m)
  }

  print(DimPlot(integrated, reduction = reduc_name, group.by = "id", raster = FALSE, cols = id_cols))
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_id.png", clust_name)), height = 10, width = 11)
  print(DimPlot(integrated, reduction = reduc_name, group.by = clust_name, raster = FALSE, label = TRUE) + NoLegend())
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_numbered.png", clust_name)), height = 10, width = 11)
  print(DimPlot(integrated, reduction = reduc_name, group.by = "id", split.by = "source", raster = FALSE, cols = id_cols))
  ggsave(file.path(outdir, "DimPlots", sprintf("DimPlot_%s_id-source.png", clust_name)), height = 10, width = 22)
  
  remove(reduc, reduc_name, clust_name)
  
}

DimPlot(integrated, reduction = "umap.scvi", group.by = "id", raster = FALSE, cols = id_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_scvi_id.png"), width = 11, height = 10)


# ---- labelling clusters ----
integrated@meta.data <- integrated@meta.data %>%
  mutate(cell_type = case_when(
    
    # micro resolution
    scvi_clusters_2 %in% c(34) ~ "delta_INS+",
    scvi_clusters_2 %in% c(41) ~ "delta_GCG+",
    
    # macro resolution
    scvi_clusters %in% c(1,2,9) ~ "beta",  # INS
    scvi_clusters %in% c(0,5,13) ~ "alpha", # GCG
    scvi_clusters %in% c(4) ~ "delta", # SST
    scvi_clusters %in% c(11) ~ "ppy", # PPY
    scvi_clusters %in% c(19)    ~ "epsilon", # GHRL
    scvi_clusters %in% c(6) ~ "ec", # TPH1
    scvi_clusters %in% c(3,20,21) ~ "ductal", # KRT19
    scvi_clusters %in% c(7,17) ~ "acinar", # PRSS1
    scvi_clusters %in% c(8) ~ "stellate", # COL1A1
    scvi_clusters %in% c(14) ~ "stellate_q", # FABP4
    scvi_clusters %in% c(16) ~ "pericytes", # RGS5
    scvi_clusters %in% c(10) ~ "endothelial", # PECAM1
    scvi_clusters %in% c(15) ~ "immune", # PTPRC, C1QC
    scvi_clusters %in% c(18) ~ "mast", # PTPRC, TPSB2
    scvi_clusters %in% c(22) ~ "schwann", # NGFR
    scvi_clusters %in% c(12) ~ "poly",
    
    TRUE ~ scvi_clusters
  ))
integrated$cell_type <- factor(
  integrated$cell_type,
  levels = c(
    "beta", "alpha", "ec",
    "delta", "delta_INS+", "delta_GCG+", "ppy", "epsilon", "poly",
    "acinar", "ductal",
    "stellate", "stellate_q", "pericytes", "endothelial",
    "immune", "mast", "schwann"
  )
)

# make cell type colours
cell_type_cols <- DiscretePalette(length(levels(integrated$cell_type)))
names(cell_type_cols) <- levels(integrated$cell_type)

# plot cell types
DimPlot(integrated, reduction = "umap.scvi", group.by = "cell_type", label = TRUE, raster = FALSE, cols = cell_type_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_scvi_cell-type.png"), width = 12, height = 10)

DimPlot(integrated, reduction = "umap.scvi", group.by = "cell_type", raster = FALSE, cols = cell_type_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_scvi_cell-type.unlabelled.png"), width = 12, height = 10)
  
DotPlot(integrated, features = gene_markers, group.by = "cell_type", cols = plot_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "DotPlots", "DotPlot_integrated_cell-type.png"), width = 12, height = 5) 

# label source
DimPlot(integrated, reduction = "umap.scvi", raster = FALSE, group.by = "cell_type", split.by = "source")
ggsave(file.path(outdir, "DimPlots", "DimPlot_integrated_cell-type_primary-stemcell.png"), width = 15, height = 8)

# plot percentages of cell type per id
ggplot(integrated@meta.data, aes(fill=cell_type, x=id.rev)) + 
  geom_bar(position="fill", color = "grey30") +
  scale_fill_manual(values = cell_type_cols) +
  coord_flip()
ggsave(file.path(outdir, "geom-bar_stacked_id_cell-type.png"), width = 8, height = 6)

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

# ---- beta cell specific analysis ----
# subset beta cells first to remove NA
beta <- subset(integrated, cell_type == "beta")
DefaultAssay(beta) <- "RNA"
beta[["SCT"]] <- NULL
gc()

# pipeline and reanalyse just beta cells
beta <- SCTransform(beta, verbose = FALSE)
gc()

# cluster beta cells separately
beta <- RunPCA(beta, verbose = FALSE)
ElbowPlot(beta)

beta <- beta %>% 
  FindNeighbors(dims = 1:15, verbose = FALSE) %>% 
  FindClusters(resolution = 0.8, cluster.name = "unintegrated_clusters", verbose = FALSE) %>% 
  RunUMAP(dims = 1:15, reduction.name = "umap.unintegrated", verbose = FALSE)

# plot new clusters
DimPlot(beta, reduction = "umap.unintegrated", group.by = "id", raster = FALSE, cols = id_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_beta_id.png"), width = 9, height = 7)

DimPlot(beta, reduction = "umap.unintegrated", group.by = "id", split.by = "source", raster = FALSE, cols = id_cols)
ggsave(file.path(outdir, "DimPlots", "DimPlot_beta_id-source.png"), width = 15, height = 8)

'###
rpca and harmony integration overcorrects beta cell clustering
only makes one large ball, we know they are all beta cells but different
i would argue integration is not needed now
integration probably wont work if there is just one cell type?
'###

# qc markers
for (q in qc_markers) {
  print(VlnPlot(
    beta, 
    group.by = "id", 
    features = q,
    cols = id_cols
  ) + NoLegend())
  ggsave(file.path(outdir, "VlnPlots", paste0("beta_", q, ".png")), height = 8, width = 8)
  remove(q)
}

VlnPlot(
  beta, 
  assay = "RNA", 
  layer = "data", 
  group.by = "id", 
  features = "INS",
  pt.size = 0,
  cols = id_cols
) + NoLegend()
ggsave(file.path(outdir, "VlnPlots", "VlnPlot_beta-INS_replicate.RNA.png"), width = 10, height = 6)

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
  file.path(outdir, "avg-ins-expr_cell-line.RNA.tsv"),
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

DotPlot(beta, assay = "RNA", features = beta_markers, group.by = "id.rev", cols = c("beige", "brown")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(
  file.path(outdir, "DotPlots", "DotPlot_beta-markers_id.rev.png"),
  width = 14, height = 6
)

VlnPlot(beta, assay = "RNA", features = beta_markers, group.by = "id.rev", cols = id_cols, stack = TRUE, fill.by = "ident") +
  NoLegend() +
  scale_fill_manual(values = id_cols)
ggsave(
  file.path(outdir, "VlnPlots", "VlnPlot_beta-markers_id.rev.png"),
  width = 14, height = 6
)

# ---- Deseq2 ----
# analyse rna instead, see differences with sct
beta.agg_rna <- AggregateExpression(
  beta,
  assays = "RNA",
  layers = "counts",
  group.by = "id",
  normalization.method = NULL,
  scale.factor = NULL,
  margin = NULL,
  verbose = FALSE
)$RNA

# create deseq2.info
deseq2.info <- data.frame(
  sample = colnames(beta.agg_rna),
  lab = sapply(strsplit(colnames(beta.agg_rna), "\\."), `[`, 1),
  source = ifelse(
    sapply(strsplit(colnames(beta.agg_rna), "\\."), `[`, 1) %in% primary,
    "primary", "stem-cell"
  ),
  cell_line = ifelse(
    sapply(strsplit(colnames(beta.agg_rna), "\\."), `[`, 1) == "akerman",
    ifelse(
      colnames(beta.agg_rna) %in% c("akerman.D39-ZKSCAN1-WT37", "akerman.D45-H1"),
      "H1", "H3"
    ),
    "primary"
  )
)
rownames(deseq2.info) <- deseq2.info$sample

# create deseq object
dds <- DESeqDataSetFromMatrix(
  countData = round(beta.agg_rna),
  colData = deseq2.info,
  design = ~ source
)

# prefilter
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
remove(keep)

# set reference level
dds$source <- relevel(dds$source, ref = "primary")

# DGEA
dds <- DESeq(dds)
res <- results(dds)

resLFC <- lfcShrink(dds, coef = resultsNames(dds)[[2]], type = "apeglm")

# plotting
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

## ---- euc dist ----
# variance stabilise
vsd <- vst(dds, blind = FALSE)

euc_dist <- dist(t(assay(vsd)))
ed_mat <- as.matrix(euc_dist)

rownames(ed_mat) <- rownames(deseq2.info)
colnames(ed_mat) <- rownames(deseq2.info)

pheatmap::pheatmap(
  ed_mat,
  color = hcl.colors(100, palette = "viridis")
)

## ---- PCA ----
pca_data <- plotPCA(
  vsd, 
  ntop = 3000,
  intgroup = c("cell_line", "lab"), 
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = lab, shape = cell_line)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "bandesh" = "grey30",
      "fasolino" = "grey50",
      "kang" = "grey70",
      "akerman" = "mediumpurple"
    )
  ) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_classic()

ggsave(file.path(outdir, "pca", "pca_top3000hvg_beta.png"), width = 6, height = 6)

remove(percent_var)

## ---- DGEA list ----
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
    file.path(outdir, "DESeq2_SBC-vs-PBC_neg-log10-sign.rnk"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

# plot gsea results (done externally)
gsea_res <- read.csv(file.path(outdir, "gsea_loop_23-4-2026_gh", "combined_gsea_results.xlsx.csv"))

gsea_res <- gsea_res %>% 
  arrange(desc(GROUP), desc(NES)) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME))) %>% 
  mutate(barplot_col = case_when(
    FDR.q.val < 0.05 & NES < 0 ~ "turquoise4",
    FDR.q.val < 0.05 & NES > 0 ~ "tomato2",
    TRUE ~ "grey50"
  ))

ggplot(gsea_res, aes(x = NAME, y = NES, fill = barplot_col)) +
  geom_col() +
  scale_fill_identity() +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme_classic() 

ggsave(file.path(outdir, "geom-col_deseq2_SBC-vs-PBC_gsea-res.png"), width = 10, height = 10)

## ---- volcano ----
res_df <- res_df %>% 
  mutate(neglog10.padj = -log10(padj + 1e-300)) %>% 
  mutate(volcano_col = case_when(
    log2FoldChange > 1 & 
      neglog10.padj > -log10(0.05 + 1e-300) ~ "tomato2",
    log2FoldChange < -1  &
      neglog10.padj > -log10(0.05 + 1e-300) ~ "turquoise4",
    TRUE ~ "grey50"
  ))

top10 <- res_df %>% 
  arrange(padj) %>% 
  slice_head(n = 10) %>% 
  pull(gene)

ggplot(res_df, aes(x = log2FoldChange, y = neglog10.padj, colour = volcano_col)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_identity() +
  geom_hline(yintercept = -log10(0.05 + 1e-300), linetype = 2) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_label_repel(
    data = subset(res_df, gene %in% c(beta_markers, top10, "INS")),
    aes(label = gene),
    size = 6,
    max.overlaps = Inf
  ) +
  theme_classic()

ggsave(file.path(outdir, "volcanos", "deseq2_SBC-vs-PBC.png"), width = 10, height = 10)

## ---- ClusterProfiler ----
# save sig up (SBC enriched) and sig down (pri enriched) genes
dgea_res <- list(
  "up" = dplyr::filter(res_df, volcano_col == "tomato2"),
  "down" = dplyr::filter(res_df, volcano_col == "turquoise4")
)

# create empty list to save results
ora_res <- list()

# perform ORA for up and down genes
for (i in seq_along(dgea_res)) {
  
  df <- dgea_res[[i]]
  nm <- names(dgea_res)[[i]]
  
  # add entrez ids for enrichGO
  entrez <- bitr(
    df$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # ORA
  ora <- enrichGO(
    gene = entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # save output
  ora_res[[nm]] <- ora
  
  remove(i, df, nm, entrez, ora)
}

options(enrichplot.colours = hcl.colors(100))

dotplot(ora_res$up, showCategory = 10)
ggsave(file.path(outdir, "DotPlots", "deseq2_SBC-vs-PBC_ORA_up.png"))

dotplot(ora_res$down, showCategory = 10)
ggsave(file.path(outdir, "DotPlots", "deseq2_SBC-vs-PBC_ORA_down.png"))

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
  file.path(outdir, "pheatmaps", "lickert-contaminants_SBC-vs-PBC.RNA.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

mat <- lickert_contaminants |>
  column_to_rownames("gene") |>
  dplyr::select(ends_with("fc")) |>
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
  color = hcl.colors(100, palette = "viridis", rev = TRUE)
); dev.off()

remove(annotation_col, annotation_row)

# ---- Perform DEA ----
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

# compare sc beta and pri beta to everything else
integrated$beta <- ifelse(
  integrated$cell_type == "beta",
  integrated$source,
  "other"
)
integrated$beta <- factor(integrated$beta)

beta.all_markers <- FindAllMarkers(
  integrated,
  assay = "RNA",
  group.by = "beta",
  verbose = FALSE
)

# ---- save output ----
saveRDS(integrated, file.path(outdir, "integrated.v4.RDS"))