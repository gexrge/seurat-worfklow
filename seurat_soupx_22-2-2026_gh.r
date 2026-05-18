#!/usr/bin/env Rscript
# https://github.com/constantAmateur/SoupX

library(Seurat)
library(SoupX)
library(ggplot2)

seurat_soupx_23.2.2026_gh <- function(
    sobj,
    nobj,
    FindNeighbors.dims,
    FindClusters.res,
    sobj.raw, 
    sobj.filt,
    outdir
  ){

  cat(">>> Running SoupX\n")
  
  # ---- SoupX ----
  sc <- SoupChannel(sobj.raw, sobj.filt) 
  sc <- setClusters(sc, sobj$seurat_clusters)
  sc <- tryCatch(
    autoEstCont(sc, verbose = FALSE, doPlot = FALSE, forceAccept = TRUE),
    error = function(e) {
      cat("  - WARNING: SoupX autoEstCont failed:", conditionMessage(e), "\n")
      cat("  - WARNING: Falling back to a manual contamination fraction of 0.2\n")
      setContaminationFraction(sc, contFrac = 0.2, forceAccept = TRUE)
    }
  )
  out <- tryCatch(
    adjustCounts(sc, verbose = FALSE),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("Contamination fractions must have already been calculated/set", msg)) {
        cat("  - WARNING: adjustCounts failed because contamination not set. Setting contamination to 0.2 and retrying.\n")
        sc <<- setContaminationFraction(sc, contFrac = 0.2, forceAccept = TRUE)
        return(adjustCounts(sc, verbose = FALSE))
      }
      stop(e)
    }
  )
  
  # Find marker genes vs soup genes
  cntSoggy = rowSums(sc$toc > 0)
  cntStrained = rowSums(out > 0)

  cat("  - Top 10 estimated marker genes (high expression, low correction):\n")
  est_markers <- sort(rowSums(out), decreasing = TRUE)
  print(names(head(est_markers, 10)))
  est_markers_df <- as.data.frame(est_markers)
  write.table(
    est_markers_df, 
    file = file.path(outdir, paste0(nobj, "_soupX_estimated_markers.txt")), 
    quote = FALSE, 
    sep = "\t", 
    col.names = NA
  )

  cat("  - Top 10 estimated soup genes (high contamination removed):\n")
  est_soup <- sort(rowSums(sc$toc - out), decreasing = TRUE)
  print(names(head(est_soup, 10)))
  est_soup_df <- as.data.frame(est_soup)
  write.table(
    est_soup_df,
    file = file.path(outdir, paste0(nobj, "_soupX_estimated_soup.txt")),
    quote = FALSE,
    sep = "\t",
    col.names = NA
  )

  # Plot genes of interest
  markers <- c(
    alpha = "GCG",
    beta = "INS",
    delta = "SST",
    pp = "PPY",
    epsilon = "GHRL",
    empdrop = "MALAT1"
  )
  
  # Filter islet markers by availability to avoid errors 
  markers_filt <- markers[markers %in% rownames(sobj)]
  
  # Report removed markers
  removed <- setdiff(markers, markers_filt)
  if (length(removed) > 0) {
    cat("  - Removed markers (not found in dataset):", paste(removed, collapse = ", "), "\n")
  }
  
  # Extract Seurat UMAP coords
  umap <- sobj@reductions$umap@cell.embeddings
  sc <- setDR(sc, umap)
  
  # Loop through available markers, plot difference
  cat("  - Plotting change in expression for islet markers:\n")
  print(FeaturePlot(sobj, features = markers_filt))
  for (marker in markers_filt) {
    print(plotChangeMap(sc, out, marker) + 
      labs(title = sprintf("Change in expression due to soup correction - %s", marker)))
  }

  # Recreate seurat object with corrected counts for downstream analysis
  cat("  - Creating new Seurat object with corrected counts\n")
  sobj <- CreateSeuratObject(
    counts = out,
    min.cells = 3,
    min.features = 200,
    project = nobj
  )

  return(sobj)

}