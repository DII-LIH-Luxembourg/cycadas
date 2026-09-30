# Shared application state ----------------------------------------------------
# One state object is created per session in app_server() and passed to every
# module. Modules read from it and write to it; nothing is kept in global
# variables.

new_app_state <- function() {
  reactiveValues(
    expr = NULL,          # expression table, see fct_expression.R
    cell_freq = NULL,     # cluster frequencies
    markers = character(0),
    umap = NULL,          # UMAP coordinates + raw marker values
    th = NULL,            # thresholds, see fct_thresholds.R
    graph = NULL,         # annotation tree, see fct_tree.R
    md = NULL,            # sample metadata
    counts = NULL,        # cluster x sample counts
    sce = NULL,           # CATALYST SingleCellExperiment
    meta_level = NULL,    # CATALYST clustering level in use
    dataset_version = 0,  # bumped whenever a new dataset is loaded
    reset_version = 0     # bumped when the workspace is cleared
  )
}

# Load a new dataset ----------------------------------------------------------
# Replaces expression data, UMAP, thresholds and the tree. `data` is the
# list returned by the read_*() functions in fct_import.R. Thresholds are
# estimated unless `th` is given.
load_dataset <- function(state, data, th = NULL) {

  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Loading data...", value = 0.2)

  markers <- data$markers
  expr <- createExpressionDF(data$expr, data$cell_freq, markers)

  progress$set(message = "Building the UMAP...", value = 0.4)
  # UMAP needs more rows than neighbours
  umap <- if (nrow(expr) > 20) buildUMAP(expr[, raw_markers(markers)]) else NULL

  progress$set(message = "Estimating thresholds...", value = 0.7)
  th <- if (is.null(th)) kmeansTH(expr[, markers]) else add_estimates(th, expr, markers)

  state$markers <- markers
  state$cell_freq <- data$cell_freq
  state$expr <- expr
  state$umap <- umap
  state$th <- th
  state$graph <- initTree()
  state$dataset_version <- isolate(state$dataset_version) + 1

  invisible(state)
}

# Re-assign all clusters after thresholds or the tree changed -----------------
rebuild_annotation <- function(state) {
  req(state$expr, state$graph, state$th)
  state$expr$cell <- rebuiltTree(state$graph, state$expr, state$th, state$markers)
}

# Clear everything ------------------------------------------------------------
reset_app_state <- function(state) {
  for (field in c("expr", "cell_freq", "umap", "th", "graph", "md", "counts",
                  "sce", "meta_level")) {
    state[[field]] <- NULL
  }
  state$markers <- character(0)
  state$dataset_version <- isolate(state$dataset_version) + 1
  state$reset_version <- isolate(state$reset_version) + 1
}

# Workspace <-> state ---------------------------------------------------------
has_df <- function(x) is.data.frame(x) && nrow(x) > 0 && ncol(x) > 0
has_val <- function(x) !is.null(x) && length(x) > 0

workspace_from_state <- function(state) {
  df_or_null <- function(x) if (has_df(x)) x else NULL

  list(
    schema_version = WORKSPACE_SCHEMA_VERSION,
    app_version    = as.character(utils::packageVersion("cycadas")),
    saved_at       = Sys.time(),
    median_expr    = df_or_null(state$expr),
    lineage_marker = state$markers,
    lineage_marker_raw = raw_markers(state$markers),
    cl_freq        = df_or_null(state$cell_freq),
    metadata       = df_or_null(state$md),
    thresholds     = df_or_null(state$th),
    counts_table   = df_or_null(state$counts),
    nodes_df       = if (!is.null(state$graph)) serialize_nodes(state$graph$nodes) else NULL,
    edges_df       = df_or_null(state$graph$edges),
    annotation     = state$graph$nodes$label,
    umap_coords    = df_or_null(state$umap)
  )
}

restore_workspace <- function(state, ws) {

  reset_app_state(state)

  if (!is.null(ws$median_expr) && !is.null(ws$cl_freq)) {
    # workspaces saved before 1.5 may hold fractions
    ws$median_expr$freq <- freq_percent(ws$median_expr$freq)
    state$expr <- ws$median_expr
    state$cell_freq <- ws$cl_freq
    state$markers <- ws$lineage_marker
    state$graph <- initTree()
  }
  if (!is.null(ws$umap_coords)) state$umap <- ws$umap_coords
  if (!is.null(ws$counts_table)) state$counts <- ws$counts_table
  if (!is.null(ws$metadata)) state$md <- ws$metadata

  if (!is.null(ws$thresholds) && !is.null(ws$median_expr)) {
    state$th <- add_estimates(drop_index_column(ws$thresholds), ws$median_expr, ws$lineage_marker)
  } else if (!is.null(ws$thresholds)) {
    state$th <- drop_index_column(ws$thresholds)
  } else if (!is.null(ws$median_expr)) {
    state$th <- kmeansTH(ws$median_expr[, ws$lineage_marker])
  }

  if (!is.null(ws$nodes_df) && !is.null(ws$edges_df)) {
    state$graph <- getGraphFromLoad(ws$nodes_df, ws$edges_df)
  }

  state$dataset_version <- isolate(state$dataset_version) + 1
  invisible(state)
}
