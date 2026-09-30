# Data import -----------------------------------------------------------------
# Every reader returns list(expr, cell_freq, markers), the input that
# load_dataset() expects:
#   expr       one row per cluster, one column per marker (raw values)
#   cell_freq  data.frame with at least the column clustering_prop
#   markers    character vector of marker names, in column order of `expr`

# GigaSOM / FlowSOM CSV export -------------------------------------------------
read_gigasom_csv <- function(path_expr, path_freq) {
  expr <- read.csv(path_expr)
  list(expr = expr, cell_freq = read.csv(path_freq), markers = colnames(expr))
}

# RemoteSOM JSON export --------------------------------------------------------
read_remotesom_json <- function(path_features, path_counts, path_medians) {

  markers <- read_json(path_features, simplifyVector = TRUE)
  if (!length(markers)) {
    stop("feature-names.json could not be parsed.")
  }

  counts <- read_json(path_counts, simplifyVector = TRUE)
  cell_freq <- data.frame(cluster = seq_along(counts),
                          clustering_prop = counts / sum(counts))

  expr <- as.data.frame(read_json(path_medians, simplifyVector = TRUE))
  colnames(expr) <- markers

  list(expr = expr, cell_freq = cell_freq, markers = markers)
}

# Demo data shipped with the package ------------------------------------------
demo_dataset <- function() {
  expr <- getExportedValue("cycadas", "df_expr_demoData")
  list(expr = expr,
       cell_freq = getExportedValue("cycadas", "cluster_freq_demoData"),
       markers = colnames(expr))
}

# Drop the row-name column that write.csv() adds
drop_index_column <- function(df) {
  df$X <- NULL
  df
}

# CATALYST ---------------------------------------------------------------------
catalyst_available <- function() {
  requireNamespace("CATALYST", quietly = TRUE) &&
    requireNamespace("SingleCellExperiment", quietly = TRUE)
}

# Cluster frequencies of a CATALYST object at a clustering level
catalyst_cell_freq <- function(sce, level = NULL) {
  cell_freq <- as.data.frame(table(c_id = CATALYST::cluster_ids(sce, level)))
  cell_freq$clustering_prop <- cell_freq$Freq / sum(cell_freq$Freq)
  cell_freq$Freq <- NULL
  colnames(cell_freq) <- c("cluster", "clustering_prop")
  cell_freq
}

# Initial dataset from a CATALYST object: the SOM codes
read_catalyst <- function(sce) {
  rd <- SummarizedExperiment::rowData(sce)
  list(expr = as.data.frame(sce@metadata$SOM_codes),
       cell_freq = catalyst_cell_freq(sce),
       markers = rd$marker_name[rd$marker_class == "type"])
}

# Meta-clustering levels offered in the UI
catalyst_meta_levels <- function(sce) {
  mc_levels <- colnames(sce@metadata$cluster_codes)
  # remove the meta2 to meta5 levels
  mc_levels[-(2:5)]
}

# Dataset at a meta-clustering level: SOM codes averaged per meta-cluster
catalyst_meta_level_data <- function(sce, level, markers) {

  somExpr <- as.data.frame(sce@metadata$SOM_codes)
  somExpr$metaLevels <- sce@metadata$cluster_codes[[level]]

  somExpr <- somExpr %>%
    group_by(metaLevels) %>%
    summarise_if(is.numeric, mean)

  list(expr = as.data.frame(somExpr[, markers]),
       cell_freq = catalyst_cell_freq(sce, level),
       markers = markers)
}

# Write the annotation back into the CATALYST object
merge_catalyst_annotation <- function(sce, expr, level) {
  merge_table <- data.frame(old_cluster = as.numeric(rownames(expr)),
                            new_cluster = expr$cell)
  CATALYST::mergeClusters(sce, k = level, table = merge_table,
                          id = paste0("merging_", level))
}
