# Differential abundance ------------------------------------------------------
# counts:  one row per cluster (aligned with the expression table), one column
#          per sample
# md:      metadata with at least the columns sample_id and condition
# cells:   the `cell` annotation column of the expression table

# Sum the cluster counts per annotation label ---------------------------------
aggregate_counts_by_cell <- function(counts, cells) {

  if (is.null(cells)) {
    stop("Missing 'cell' labels in expression table.")
  }
  if (nrow(counts) != length(cells)) {
    stop("Row mismatch between counts table and expression table. Check inputs.")
  }
  counts$cell <- cells
  agg <- aggregate(. ~ cell, counts, sum, na.rm = TRUE)
  rownames(agg) <- agg$cell
  agg$cell <- NULL
  agg
}

# Column-wise percentages ------------------------------------------------------
to_percent <- function(counts) {
  sweep(counts, 2, colSums(counts, na.rm = TRUE), "/") * 100
}

# Append "_remaining" to labels of nodes that have children -------------------
label_remaining <- function(labels, graph) {
  if (is.null(graph)) return(labels)
  vapply(labels, function(x) {
    nid <- graph$nodes$id[graph$nodes$label == x]
    if (length(nid) > 0 && node_has_children(graph, nid)) paste0(x, "_remaining") else x
  }, character(1), USE.NAMES = FALSE)
}

# Proportion table of the merged phenotypes -----------------------------------
merged_prop_table <- function(counts, cells, graph) {

  agg <- aggregate_counts_by_cell(counts, cells)
  rownames(agg) <- label_remaining(rownames(agg), graph)

  as.matrix(to_percent(agg))
}

# Pairwise Wilcoxon test per phenotype ----------------------------------------
# Signals warning() for dropped samples and stop() for unusable input, so that
# callers can decide how to report them.
run_da <- function(counts, md, cells, graph, correction_method = "holm") {

  if (!all(c("sample_id", "condition") %in% colnames(md))) {
    stop("Metadata must contain 'sample_id' and 'condition'.")
  }
  md$sample_id <- as.character(md$sample_id)

  # --- 1) Counts aggregated by cell label ----
  countsTable <- aggregate_counts_by_cell(counts, cells)

  # --- 2) Drop samples with zero counts ----
  cs <- colSums(countsTable, na.rm = TRUE)
  if (any(cs == 0)) {
    warning(sprintf("Some samples have zero total counts and will be dropped: %s",
                    paste(names(cs)[cs == 0], collapse = ", ")), call. = FALSE)
  }
  if (!any(cs > 0)) {
    stop("All samples have zero counts; cannot compute proportions.")
  }
  props_table <- to_percent(countsTable[, cs > 0, drop = FALSE])

  # --- 3) Match samples to metadata ----
  mm <- match(colnames(props_table), md$sample_id)
  if (anyNA(mm)) {
    warning(sprintf("Samples missing in metadata and dropped: %s",
                    paste(colnames(props_table)[is.na(mm)], collapse = ", ")), call. = FALSE)
  }
  keep <- which(!is.na(mm))
  if (length(keep) == 0) {
    stop("No overlap between counts table columns and metadata sample_id.")
  }
  props_table <- props_table[, keep, drop = FALSE]
  tmp_cond <- droplevels(as.factor(md$condition[mm[keep]]))

  if (nlevels(tmp_cond) < 2) {
    stop("Need at least two conditions to run DA tests.")
  }

  # --- 4) Pairwise Wilcoxon per row (cell) ----
  DA_df <- lapply(seq_len(nrow(props_table)), function(i) {
    y <- as.numeric(props_table[i, ])
    # Skip if all NA or constant
    if (all(!is.finite(y)) || length(unique(y[is.finite(y)])) < 2) return(NULL)
    # Ensure each group has at least one finite value
    if (any(!tapply(y, tmp_cond, function(v) sum(is.finite(v)) > 0))) return(NULL)

    # ties only affect the exactness of the p-value, not worth reporting
    pw <- tryCatch(
      suppressWarnings(pairwise.wilcox.test(y, tmp_cond, p.adjust.method = correction_method)),
      error = function(e) NULL
    )
    if (is.null(pw) || is.null(pw$p.value)) return(NULL)

    m <- reshape2::melt(pw$p.value, varnames = c("Cond1", "Cond2"), value.name = "p.value")
    m <- subset(m, is.finite(p.value))
    if (nrow(m) == 0) return(NULL)

    m$Cell <- rownames(props_table)[i]
    m$Naming <- label_remaining(m$Cell, graph)
    m
  })
  DA_df <- do.call(rbind, DA_df)

  if (is.null(DA_df)) {
    warning("No pairwise results computed (check inputs/conditions).", call. = FALSE)
    DA_df <- data.frame(Cond1 = character(0), Cond2 = character(0), p.value = numeric(0),
                        Cell = character(0), Naming = character(0))
  }

  colnames(DA_df) <- c("Cond1", "Cond2", "p-value", "Cell", "Naming")
  DA_df
}

# Per-sample proportion of a node and all its descendants ---------------------
# Returns data.frame(value, cond) with one row per sample.
node_proportions <- function(counts, cells, graph, node_id, md) {

  selected_labels <- node_with_descendant_labels(graph, node_id)

  props_table <- to_percent(aggregate_counts_by_cell(counts, cells))
  props_table <- props_table[rownames(props_table) %in% selected_labels, , drop = FALSE]

  props <- data.frame(value = colSums(props_table))
  mm <- match(rownames(props), md$sample_id)
  props$cond <- as.factor(md$condition[mm])
  props
}

# Pairwise Wilcoxon test for a single node ------------------------------------
node_da <- function(props, label) {

  pw <- suppressWarnings(pairwise.wilcox.test(props$value, props$cond, p.adjust.method = "none"))
  df <- subset(reshape2::melt(pw$p.value), value != 0)
  df$cell <- rep(label, nrow(df))

  colnames(df) <- c("Var1", "Var2", "p-value", "Cell")
  df
}

# Function to get all pairs of combinations
get_pairs <- function(vec) {
    pairs <- combn(vec, 2)
    pairs_list <- split(pairs, col(pairs))
    return(pairs_list)
}
