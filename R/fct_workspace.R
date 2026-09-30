# Workspace files -------------------------------------------------------------
# A workspace is an RDS file holding a list of plain R objects:
#   schema_version, app_version, saved_at, median_expr, lineage_marker,
#   lineage_marker_raw, cl_freq, metadata, thresholds, counts_table,
#   nodes_df, edges_df, annotation, umap_coords

WORKSPACE_SCHEMA_VERSION <- "1.0"

save_workspace <- function(path, state) {
  stopifnot(is.list(state), !is.null(state$schema_version))
  saveRDS(state, file = path, compress = "xz")
  invisible(path)
}

load_workspace <- function(path) {
  ws <- if (grepl("\\.qs$", path)) qs::qread(path) else readRDS(path)

  # minimal validation
  stopifnot(is.list(ws), !is.null(ws$schema_version))

  ws
}

# Export the annotation (expression table, nodes, edges) as a zip of CSVs -----
export_annotation_zip <- function(file, expr, graph) {

  temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
  dir.create(temp_directory)

  download_list <- list(annTable = expr,
                        nodesTable = serialize_nodes(graph$nodes),
                        edgesTable = graph$edges)

  download_list %>%
    imap(function(x,y){
      if(!is.null(x)){
        readr::write_csv(x, file.path(temp_directory, glue("{y}_data.csv")))
      }
    })

  zip::zip(
    zipfile = file,
    files = dir(temp_directory),
    root = temp_directory
  )
}
