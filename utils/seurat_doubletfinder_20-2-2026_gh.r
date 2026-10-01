#!/usr/bin/env Rscript
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

# Set libraries
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(SimDesign)

seurat_doubletfinder_20.2.2026_gh <- function(
  sobj,
  nobj,
  outdir,
  FindNeighbors.dims,
  seq_method
  ){
  
  # ---- Run DoubletFinder ----
  cat(">>> Running DoubletFinder\n")

  # pK Identification (no ground-truth) (hide unnecessary output)
  sweep.res <- quiet(paramSweep(sobj, PCs = FindNeighbors.dims, sct = FALSE))  # noisy
  sweep.stats <- quiet(summarizeSweep(sweep.res, GT = FALSE))
  pK_results <- quiet(find.pK(sweep.stats))
  best.pK <- as.numeric(as.character(pK_results$pK[which.max(pK_results$BCmetric)]))

  # Calculate doublet rate based on seq method
  if (seq_method == "plate") {
    doublet_rate <- 0.01    # https://www.biorxiv.org/content/10.1101/632216v2.full
    cat(sprintf("  - seq_method set to plate, running droplet rate at: %.3f\n", doublet_rate))
  } else if (seq_method == "droplet") {
    doublet_rate <- 0.0075 * (nrow(sobj@meta.data) / 1000)    # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
    cat(sprintf("  - seq_method set to droplet, running droplet rate at: %.3f\n", doublet_rate))
  } else{
    doublet_rate <- 0.075
    cat(sprintf("  WARNING: seq_method not set to plate or droplet, running droplet rate at: %.3f\n", doublet_rate))
  }

  # Homotypic Doublet Proportion Estimate
  annotations <- sobj@meta.data$seurat_clusters
  homotypic.prop <- quiet(modelHomotypic(annotations))
  nExp_poi <- round(doublet_rate * nrow(sobj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder with varying classification stringencies
  sobj <- quiet(doubletFinder(
      sobj, 
      PCs = FindNeighbors.dims, 
      pN = 0.25, 
      pK = best.pK, 
      nExp = nExp_poi.adj, 
      reuse.pANN = NULL, 
      sct = FALSE
  ))
  
  # Rename inconsistent doubletfinder metadata columns
  colnames(sobj@meta.data)[grep("pANN", colnames(sobj@meta.data))] <- "pANN"
  colnames(sobj@meta.data)[grep("DF.classifications", colnames(sobj@meta.data))] <- "doublet_call"
  
  # Plot doublets 
  print(DimPlot(sobj, reduction = "umap", label = TRUE) + NoLegend())
  print(DimPlot(sobj, group.by = "doublet_call"))

  # Drop doublets and invalid scale.data slot to save memory
  doublet_count <- sum(sobj$doublet_call == "Doublet")
  cat(sprintf("  - Removing doublets: %i cells (out of %i)\n", doublet_count, ncol(sobj)))
  sobj <- subset(sobj, subset = doublet_call == "Singlet")
  sobj[["RNA"]]$scale.data <- NULL
  
  return(list(
    sobj = sobj,
    params = list(
      seq_method = seq_method,
      doublet_rate = doublet_rate,
      best.pK = best.pK,
      nExp_poi = nExp_poi,
      nExp_poi.adj = nExp_poi.adj,
      doublet_count = doublet_count
    )
  ))

}
