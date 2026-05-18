#!/usr/bin/env Rscript
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

# Set libraries
library(Seurat)

seurat_pipe_20.2.2026_gh <- function(
    sobj,
    FindClusters.res,
    FindNeighbors.dims
  ) {

  cat(">>> Running Seurat pipeline\n")
  
  # ---- SCTransform data ----
  cat("  - Running SCTransform\n")
  sobj <- SCTransform(sobj, vars.to.regress = "percent.mt", verbose = FALSE)
  
  # ---- linear dimensional reduction ----
  sobj <- RunPCA(sobj, verbose = FALSE)
  print(ElbowPlot(sobj)) 
  
  # ---- Cluster the cells ----
  cat("  - Clustering cells with resolution:", FindClusters.res, "\n")
  sobj <- FindNeighbors(sobj, dims = FindNeighbors.dims, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = FindClusters.res, verbose = FALSE)
  
  # ---- Run non-linear dimensional reduction ----
  cat("  - Plotting UMAPs with dimensions:", max(FindNeighbors.dims),"\n")
  sobj <- RunUMAP(sobj, dims = FindNeighbors.dims, verbose = FALSE)
  
  return(sobj)
}
  