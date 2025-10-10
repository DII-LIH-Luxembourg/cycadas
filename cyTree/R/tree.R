#' Create an empty annotation tree
#'
#' @return A list with `nodes` and `edges` data frames.
#' @export
initTree <- function() {
  nodes <- data.frame(
    id = 1L,
    label = "Unassigned",
    pm = I(list(character())),
    nm = I(list(character())),
    color = "blue",
    stringsAsFactors = FALSE
  )

  edges <- data.frame(
    from = integer(),
    to = integer(),
    stringsAsFactors = FALSE
  )

  list(nodes = nodes, edges = edges)
}

#' Add a node to an annotation tree
#'
#' @param graph A tree object as returned by [initTree()].
#' @param parent The label of the parent node.
#' @param name The label for the new node.
#' @param pos_m A list or vector of positive markers.
#' @param neg_m A list or vector of negative markers.
#' @param color The display color for the node.
#'
#' @return The updated tree.
#' @export
add_node <- function(graph, parent, name, pos_m, neg_m, color = "blue") {
  if (is.null(graph) || !is.list(graph)) {
    stop("`graph` must be a tree created by initTree().", call. = FALSE)
  }

  nodes <- graph$nodes
  edges <- graph$edges

  if (!parent %in% nodes$label) {
    stop("Parent node not found in tree.", call. = FALSE)
  }

  next_id <- if (nrow(nodes) == 0) 1L else max(nodes$id) + 1L
  parent_id <- nodes$id[nodes$label == parent][1]

  new_row <- data.frame(
    id = next_id,
    label = name,
    pm = I(list(as.character(unlist(pos_m)))),
    nm = I(list(as.character(unlist(neg_m)))),
    color = color,
    stringsAsFactors = FALSE
  )

  graph$nodes <- rbind(nodes, new_row)
  new_edge <- data.frame(from = next_id, to = parent_id, stringsAsFactors = FALSE)
  graph$edges <- rbind(edges, new_edge)
  graph
}

#' Recursively delete a node and its children
#'
#' @param graph_data Tree data structure.
#' @param node_id Node identifier.
#'
#' @return The modified tree.
#' @export
delete_child_nodes <- function(graph_data, node_id) {
  children <- graph_data$edges$from[graph_data$edges$to == node_id]

  if (length(children) > 0) {
    for (child in children) {
      graph_data <- delete_child_nodes(graph_data, child)
    }
  }

  graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, , drop = FALSE]
  keep_edges <- !(graph_data$edges$from == node_id | graph_data$edges$to == node_id)
  graph_data$edges <- graph_data$edges[keep_edges, , drop = FALSE]
  graph_data
}

#' Return all children for a node
#'
#' @param graph_data Tree data structure.
#' @param node_id Node identifier.
#'
#' @return A vector of node ids representing descendants.
#' @export
all_my_children <- function(graph_data, node_id) {
  if (node_id <= 1) {
    return(integer())
  }

  children <- graph_data$edges$from[graph_data$edges$to == node_id]
  if (!length(children)) {
    return(children)
  }

  idx <- 1L
  while (idx <= length(children)) {
    current <- children[idx]
    new_children <- graph_data$edges$from[graph_data$edges$to == current]
    if (length(new_children)) {
      children <- unique(c(children, new_children))
    }
    idx <- idx + 1L
  }
  children
}

#' Reconstruct tree data from stored tables
#'
#' @param df_nodes Data frame describing nodes.
#' @param df_edges Data frame describing edges.
#'
#' @return A tree compatible with [initTree()].
#' @export
getGraphFromLoad <- function(df_nodes, df_edges) {
  if (!nrow(df_nodes)) {
    return(initTree())
  }

  df_nodes$pm[is.na(df_nodes$pm)] <- ""
  df_nodes$nm[is.na(df_nodes$nm)] <- ""

  df_nodes$pm <- lapply(strsplit(df_nodes$pm, "\\|", fixed = FALSE), function(x) x[nzchar(x)])
  df_nodes$nm <- lapply(strsplit(df_nodes$nm, "\\|", fixed = FALSE), function(x) x[nzchar(x)])

  df_nodes$pm <- I(df_nodes$pm)
  df_nodes$nm <- I(df_nodes$nm)

  df_edges <- df_edges[, c("from", "to")]
  df_edges$from <- as.integer(df_edges$from)
  df_edges$to <- as.integer(df_edges$to)

  list(nodes = df_nodes, edges = df_edges)
}

#' Delete a leaf node from the tree
#'
#' @param graph_data Tree data structure.
#' @param node_id Identifier of the node to delete.
#' @param notify Optional callback invoked when the node is not a leaf. The
#'   callback receives no arguments.
#'
#' @return The updated tree.
#' @export
delete_leaf_node <- function(graph_data, node_id, notify = NULL) {
  children <- graph_data$edges$from[graph_data$edges$to == node_id]

  if (length(children) > 0) {
    if (!is.null(notify) && is.function(notify)) {
      notify()
    }
    return(graph_data)
  }

  graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, , drop = FALSE]
  keep_edges <- !(graph_data$edges$from == node_id | graph_data$edges$to == node_id)
  graph_data$edges <- graph_data$edges[keep_edges, , drop = FALSE]
  graph_data
}

#' Rebuild the annotation for all clusters based on the stored tree
#'
#' @param graph Tree returned by [initTree()].
#' @param df_expr Expression data frame that includes a `cell` column.
#' @param th Threshold table.
#' @param markers Character vector of marker column names.
#' @param filter_fun Function used to filter clusters based on markers. The
#'   function must accept `(df, pos, neg, th)` and return a data frame subset of
#'   `df`.
#'
#' @return A character vector with the rebuilt annotations.
#' @export
rebuiltTree <- function(graph, df_expr, th, markers, filter_fun) {
  if (is.null(graph) || !is.list(graph) || !length(graph$nodes)) {
    return(df_expr$cell)
  }

  df_expr$cell <- "Unassigned"
  nodes <- graph$nodes
  edges <- graph$edges

  if (nrow(nodes) <= 1) {
    return(df_expr$cell)
  }

  parent_lookup <- stats::setNames(edges$to, edges$from)

  for (i in seq_len(nrow(nodes))) {
    node_id <- nodes$id[i]
    if (node_id == 1L) {
      next
    }

    parent_id <- parent_lookup[as.character(node_id)]
    if (is.na(parent_id)) {
      next
    }

    parent_label <- nodes$label[nodes$id == parent_id]
    if (!length(parent_label)) {
      next
    }

    posMarker <- unique(as.character(unlist(nodes$pm[[i]])))
    negMarker <- unique(as.character(unlist(nodes$nm[[i]])))
    posMarker <- posMarker[nzchar(posMarker)]
    negMarker <- negMarker[nzchar(negMarker)]

    tmp_parent <- df_expr[df_expr$cell == parent_label, , drop = FALSE]
    if (!nrow(tmp_parent)) {
      next
    }

    marker_df <- tmp_parent[, markers, drop = FALSE]
    filtered <- filter_fun(marker_df, posMarker, negMarker, th)

    if (nrow(filtered)) {
      df_expr[rownames(filtered), "cell"] <- nodes$label[i]
    }
  }

  df_expr$cell
}
