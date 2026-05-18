#!/usr/bin/env Rscript
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

# Set libraries
library(Seurat)
library(DoubletFinder)
library(dplyr)

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
  sweep.res <- suppressMessages(paramSweep(sobj, PCs = FindNeighbors.dims, sct = TRUE))  # noisy
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  pK_results <- find.pK(sweep.stats)
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
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate * nrow(sobj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder with varying classification stringencies
  sobj <- doubletFinder(
    sobj, 
    PCs = FindNeighbors.dims, 
    pN = 0.25, 
    pK = best.pK, 
    nExp = nExp_poi.adj, 
    reuse.pANN = NULL, 
    sct = TRUE
  )
  
  # Rename inconsistent doubletfinder metadata columns
  colnames(sobj@meta.data)[grep("pANN", colnames(sobj@meta.data))] <- "pANN"
  colnames(sobj@meta.data)[grep("DF.classifications", colnames(sobj@meta.data))] <- "doublet_call"
  
  # Plot doublets 
  print(DimPlot(sobj, reduction = "umap", label = TRUE)+ NoLegend())
  print(DimPlot(sobj, group.by = "doublet_call"))
  doublets_per_cluster <- 
    table(sobj$seurat_clusters, sobj$doublet_call) %>% 
    as.data.frame() %>% 
    pivot_wider(names_from = Var2, values_from = Freq) %>% 
    rename(cluster = Var1)
  write.table(
    doublets_per_cluster,
    file = file.path(outdir, paste0(nobj,"_doubletfinder_stats.txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  # Drop doublets and now incorrect SCT layer
  doublet_count <- sum(sobj$doublet_call == "Doublet")
  cat(sprintf("  - Removing doublets: %i cells (out of %i)\n", doublet_count, ncol(sobj)))
  sobj <- subset(sobj, subset = doublet_call == "Singlet")
  DefaultAssay(sobj) <- "RNA"
  sobj[["SCT"]] <- NULL
  
  # Save doublet information
  params <- list(
    FindNeighbors.dims = max(FindNeighbors.dims),
    seq_method = seq_method,
    doublet_rate = doublet_rate,
    best.pK = best.pK,
    nExp_poi = nExp_poi,
    nExp_poi.adj = nExp_poi.adj,
    doublet_count = doublet_count
  )
  
  df <- data.frame(
    name  = names(params),
    value = unlist(params),
    row.names = NULL
  )
  
  write.table(
    df,
    file = file.path(outdir, paste0(nobj, "_doubletfinder_parameters.txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  return(sobj)

}
