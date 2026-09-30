# Expression data -------------------------------------------------------------
# Pure functions that build and transform the cluster expression table.
# The expression table ("expr") has one row per cluster and the columns
#   <marker>_raw  median expression as imported
#   <marker>      the same values scaled to 0-1
#   freq          cluster frequency (clustering_prop)
#   cell          current annotation label of the cluster

# Column names of the raw (unscaled) marker values
raw_markers <- function(markers) paste0(markers, "_raw")

# Scale the expression values between 0 and 1 ---------------------------------
normalize01 <- function(hm) {

  eDR <- as.matrix(hm)
  rng <- colQuantiles(eDR, probs = c(0.01, 0.99))
  expr01 <- t((t(eDR) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1

  return(as.data.frame(expr01))
}

# Create the overall expression table ------------------------------------------
# `df_expr` holds one column per marker, in the order given by `markers`.
createExpressionDF <- function(df_expr, cell_freq, markers) {

  df_expr <- as.data.frame(df_expr)
  df01_expr <- normalize01(df_expr)
  colnames(df01_expr) <- markers
  colnames(df_expr) <- raw_markers(markers)

  df_expr <- cbind(df_expr, df01_expr)

  ## Add frequencies and annotation
  df_expr$freq <- cell_freq$clustering_prop
  df_expr$cell <- "Unassigned"

  return(df_expr)
}

# Build UMAP and return as DF -------------------------------------------------
buildUMAP <- function(df_expr, seed = 1234) {

  set.seed(seed)
  my_umap <- umap(df_expr)

  df_umap <- data.frame(
    u1 = my_umap$layout[, 1],
    u2 = my_umap$layout[, 2],
    my_umap$data,
    cluster_number = 1:length(my_umap$layout[, 1]),
    check.names = FALSE
  )

  return(df_umap)
}

# Filter the Heatmap ----------------------------------------------------------
# Keep the rows of DF that are above the threshold for every marker in posList
# and below the threshold for every marker in negList.
filterHM <- function(DF, posList, negList, th) {

  for (marker in posList) {
    DF <- DF[DF[marker] >= th$threshold[match(marker, th$cell)], ]
  }
  for (marker in negList) {
    DF <- DF[DF[marker] < th$threshold[match(marker, th$cell)], ]
  }
  return(as.data.frame(DF))
}

# Filter the color for UMAP ---------------------------------------------------
filterColor <- function(DF, hm) {
  ifelse(rownames(DF) %in% rownames(hm), "selected phenotype", "other clusters")
}

# Set a phenotype name from marker selection ----------------------------------
# Currently not used!
setPhenotypeName <- function(markers, s, ph_name) {
  if (s == "pos") {
    ph_name <- paste0(ph_name, paste0(markers, collapse = "+"), "+")
  } else if (s == "neg") {
    ph_name <- paste0(ph_name, paste0(markers, collapse = "-"), "-")
  }
  return(ph_name)
}
