#!/usr/bin/env Rscript
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

# Set libraries
library(Seurat)

seurat_integrate_24.2.2026_gh <- function(
    sobj_list,
    nobj,
    FindClusters.res,
    FindNeighbors.dims
  ) {
  
  cat(">>> Running Seurat integration\n")
  
  # Merge data for comparison (not used for analysis, just to check batch effects before integration)
  cat("  - Merging data for batch effect check\n")
  merged <- merge(
    x = sobj_list[[1]],
    y = sobj_list[-1],
    add.cell.ids = names(sobj_list),
    project = nobj
  )
  
  # extracting experiment from cell id
  cat("  - Creating experiment identity from cell id\n")
  merged@meta.data$experiment <- sapply(strsplit(rownames(merged@meta.data), "_"), "[", 1)
  
  # Check batch effects before integration
  cat("  - \n")
  merged <- SCTransform(merged, vars.to.regress = "percent.mt", verbose = FALSE)
  merged <- RunPCA(merged, verbose = FALSE)
  print(ElbowPlot(merged))
  merged <- FindNeighbors(merged, dims = FindNeighbors.dims, verbose = FALSE)
  merged <- FindClusters(merged, resolution = FindClusters.res, cluster.name = "unintegrated_clusters", verbose = FALSE)
  merged <- RunUMAP(merged, dims = FindNeighbors.dims, reduction.name = "umap.unintegrated", verbose = FALSE)
  
  # Visualization before integration
  cat("  - Plotting UMAPs before integration\n")
  p1 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = "experiment", raster = FALSE)
  p2 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", raster = FALSE)
  print(p1 + p2)
  
  # integrate layers with CCA
  merged <- IntegrateLayers(
    object = merged,
    method = CCAIntegration,
    normalization.method = "SCT",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )
  
  # Recluster
  merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = FindNeighbors.dims, verbose = FALSE)
  merged <- FindClusters(merged, resolution = FindClusters.res, verbose = FALSE)
  merged <- RunUMAP(merged, dims = FindNeighbors.dims, reduction = "integrated.cca", verbose = FALSE)
  
  # Plot integrated results
  cat("  - Plotting integrated UMAPs\n")
  p3 <- DimPlot(merged, reduction = "umap", group.by = "experiment", raster = FALSE)
  p4 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", raster = FALSE)
  print(p3 + p4)
  
  # re-join RNA layers after integration
  merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
  
  return(merged)
}
  