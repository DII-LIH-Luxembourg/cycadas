# Folder sessions -------------------------------------------------------------
# cycadas(dir) opens the cluster data in a folder and saves the annotation
# back into it. A host application (CytoSuite) uses this to start CyCadas on
# its own clustering and to read the result. Files in the folder:
#
#   median_expr.csv     input: one row per cluster, one column per marker
#   cluster_freq.csv    input: column clustering_prop, one row per cluster
#   cluster_counts.csv  optional input: one row per cluster, one column per sample
#   metadata.csv        optional input: sample_id, condition, ... (with counts)
#   workspace.rds       written by "Save"; reopened in place of the inputs
#   labels.csv          written by "Save": cluster (row number of the input),
#                       label, path (the node's ancestors and itself, " > ")
#
# Both outputs are written to a temporary file and renamed, so a reader never
# sees a half-written file.

FOLDER_WORKSPACE <- "workspace.rds"
FOLDER_LABELS <- "labels.csv"

# Fill the state from the folder: the saved workspace if there is one,
# else the input files.
load_folder <- function(state, dir) {
  ws <- file.path(dir, FOLDER_WORKSPACE)
  if (file.exists(ws)) {
    restore_workspace(state, load_workspace(ws))
    return(invisible(state))
  }
  load_dataset(state, read_gigasom_csv(file.path(dir, "median_expr.csv"),
                                       file.path(dir, "cluster_freq.csv")))
  counts <- file.path(dir, "cluster_counts.csv")
  md <- file.path(dir, "metadata.csv")
  if (file.exists(counts) && file.exists(md)) {
    state$counts <- drop_index_column(read.csv(counts, check.names = FALSE))
    state$md <- drop_index_column(read.csv(md, check.names = FALSE))
  }
  invisible(state)
}

# Path from the first node below the root to the node with `label`
node_path <- function(graph, label) {
  nodes <- graph$nodes
  edges <- graph$edges
  id <- nodes$id[nodes$label == label][1]
  if (is.na(id)) return(label)
  path <- character(0)
  while (!is.na(id) && id != 1) {
    path <- c(nodes$label[nodes$id == id], path)
    id <- edges$to[edges$from == id][1]
  }
  if (length(path) == 0) label else paste(path, collapse = " > ")
}

# Cluster labels as written to labels.csv
folder_labels <- function(expr, graph) {
  labels <- unique(expr$cell)
  paths <- vapply(labels, function(l) node_path(graph, l), character(1))
  data.frame(cluster = seq_len(nrow(expr)), label = expr$cell,
             path = unname(paths[expr$cell]))
}

write_replacing <- function(dir, name, write) {
  tmp <- file.path(dir, paste0(".", name, ".tmp"))
  write(tmp)
  if (!file.rename(tmp, file.path(dir, name))) {
    unlink(tmp)
    stop("Could not write ", name, " in ", dir)
  }
}

# Write the workspace and the cluster labels into the folder
save_folder <- function(state, dir) {
  ws <- workspace_from_state(state)
  labels <- folder_labels(state$expr, state$graph)
  write_replacing(dir, FOLDER_WORKSPACE, function(p) save_workspace(p, ws))
  write_replacing(dir, FOLDER_LABELS, function(p) write.csv(labels, p, row.names = FALSE))
  invisible(dir)
}

# Module: load the folder at start, "Save" button in the navbar ----------------

mod_folder_ui <- function(id) {
  ns <- NS(id)
  actionButton(ns("save"), "Save", icon = icon("floppy-disk"), class = "btn-primary btn-sm",
               title = "Save the thresholds and the annotation tree into the folder this session was opened on")
}

mod_folder_server <- function(id, state, dir) {
  moduleServer(id, function(input, output, session) {
    # a failure here would end the session and leave an empty, greyed-out page
    tryCatch(isolate(load_folder(state, dir)), error = function(e) {
      showModal(modalDialog(title = "Could not open the data",
                            tags$p(conditionMessage(e)), tags$p(class = "small text-body-secondary", dir),
                            easyClose = TRUE))
    })

    observeEvent(input$save, {
      ok <- tryCatch({ save_folder(state, dir); TRUE }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
        FALSE
      })
      if (ok) showNotification(sprintf("Saved at %s.", format(Sys.time(), "%H:%M")), type = "message")
    })
  })
}
