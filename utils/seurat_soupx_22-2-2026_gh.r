#!/usr/bin/env Rscript
# https://github.com/constantAmateur/SoupX

library(Seurat)
library(SoupX)
library(ggplot2)

seurat_soupx_23.2.2026_gh <- function(
    sobj,
    nobj,
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

  cat("  - Top 10 most zeroed genes:\n")
  mostZeroed = sort((cntSoggy - cntStrained) / cntSoggy)
  print(names(tail(mostZeroed, 10)))
  mostZeroed <- as.data.frame(mostZeroed)
  write.table(
    mostZeroed, 
    file = file.path(outdir, paste0(nobj, "_soupX_mostZeroed.txt")), 
    quote = FALSE, 
    sep = "\t", 
    col.names = NA
  )

  cat("  - Top 10 most corrected genes:\n")
  mostCorrected <- sort(rowSums(sc$toc > out) / rowSums(sc$toc > 0))
  print(names(head(mostCorrected, 10)))
  mostCorrected_df <- as.data.frame(mostCorrected)
  write.table(
    mostCorrected_df,
    file = file.path(outdir, paste0(nobj, "_soupX_mostCorrected.txt")),
    quote = FALSE,
    sep = "\t",
    col.names = NA
  )

  # Plot genes of interest
  markers <- c(
    alpha = "GCG",
    beta = "INS",
    delta = "SST",
    pp = "PPY"
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
