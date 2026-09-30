# Annotation tree -------------------------------------------------------------
# The tree ("graph") is a list with
#   nodes  tibble(id, label, pm, nm, color); pm / nm are list columns holding
#          the positive / negative markers that define the node
#   edges  data.frame(from, to) pointing from a child to its parent.
#          The root node (id 1, "Unassigned") has a self-loop 1 -> 1.

# Initialize the Tree at start ------------------------------------------------
initTree <- function() {

  # create initial master node of all Unassigned clusters
  nodes <- tibble(id = 1,
                  label = "Unassigned",
                  pm = list(""),
                  nm = list(""),
                  color = "blue"
  )

  edges <- data.frame(from = c(1), to = c(1))

  return (list(nodes = nodes, edges = edges))
}

# Add a new node below `parent` -----------------------------------------------
add_node <- function(graph,parent,name,pos_m,neg_m,color) {

  next_id <- max(graph$nodes$id) + 1
  parent_row <- graph$nodes %>% dplyr::filter(label == parent)
  graph$nodes <- graph$nodes %>% add_row(id = next_id,
                                         label = name,
                                         pm = pos_m,
                                         nm = neg_m,
                                         color = color)

  graph$edges <- graph$edges %>% add_row(from = next_id, to = parent_row$id)

  return(graph)
}

# Direct children of a node ----------------------------------------------------
node_children <- function(graph_data, node_id) {
  graph_data$edges$from[graph_data$edges$to == node_id & graph_data$edges$from != node_id]
}

# TRUE if the node has an incoming edge. The root counts as having children
# because of its self-loop.
node_has_children <- function(graph_data, node_id) {
  any(graph_data$edges$to %in% node_id)
}

# Get all Children from a node selection --------------------------------------
all_my_children <- function(graph_data, node_id) {

  # we don't want the root node
  if (node_id > 1) {
    children <- graph_data$edges$from[graph_data$edges$to == node_id]
    i <- 1
    while (i <= length(children)) {
      children <- append(children, graph_data$edges$from[graph_data$edges$to == children[i]])
      i <- i + 1
    }
    return (children)
  } else {
    return (c())
  }
}

# Labels of a node and all its descendants ------------------------------------
node_with_descendant_labels <- function(graph_data, node_id) {
  ids <- c(all_my_children(graph_data, node_id), node_id)
  graph_data$nodes$label[graph_data$nodes$id %in% ids]
}

# define recursive function to delete a node and all its children -------------
delete_child_nodes <- function(graph_data, node_id) {

  for (child in node_children(graph_data, node_id)) {
    graph_data <- delete_child_nodes(graph_data, child)
  }

  # delete the node and its edges from the graph data
  graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, ]
  graph_data$edges <- graph_data$edges[!(graph_data$edges$from == node_id | graph_data$edges$to == node_id), ]

  return(graph_data)
}

# Delete a leaf node ----------------------------------------------------------
# Callers must check `node_has_children()` first; non-leaf nodes are left as is.
delete_leaf_node <- function(graph_data, node_id) {

  if (node_has_children(graph_data, node_id)) {
    return(graph_data)
  }
  graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, ]
  graph_data$edges <- graph_data$edges[!(graph_data$edges$from == node_id | graph_data$edges$to == node_id), ]

  return(graph_data)
}

# Build Graph from nodes and edges --------------------------------------------
getGraphFromLoad <- function(df_nodes, df_edges) {

  df_nodes$pm[is.na(df_nodes$pm)] <- ""
  df_nodes$pm <- strsplit(df_nodes$pm, "\\|")

  df_nodes$nm[is.na(df_nodes$nm)] <- ""
  df_nodes$nm <- strsplit(df_nodes$nm, "\\|")

  return (list(nodes = df_nodes, edges = df_edges))
}

# Collapse the pm / nm list columns to "|" separated strings for export -------
serialize_nodes <- function(nodes) {
  collapse <- function(x) vapply(x, function(v) paste(unlist(v), collapse = "|"), character(1))
  nodes$pm <- collapse(nodes$pm)
  nodes$nm <- collapse(nodes$nm)
  nodes
}

# Rebuild Tree ----------------------------------------------------------------
# Re-assign every cluster to a node, following the tree from the root down.
# Used when thresholds change or a tree is loaded. Returns the new `cell` column.
rebuiltTree <- function(graph, df_expr, th, markers) {

  df_expr$cell <- "Unassigned"

  if (nrow(graph$nodes) > 1) {

    # nodes are stored in creation order, so a parent is always handled
    # before its children
    for (i in 1:nrow(graph$nodes)) {

      nodeID <- graph$nodes$id[i]
      parentID <- graph$edges$to[match(nodeID, graph$edges$from)]

      # remove the empty strings in the markers
      posMarker <- unlist(graph$nodes$pm[i])
      negMarker <- unlist(graph$nodes$nm[i])
      posMarker <- unique(posMarker[nzchar(posMarker)])
      negMarker <- unique(negMarker[nzchar(negMarker)])

      # filter the clusters of the parent by the node's markers
      parent_label <- graph$nodes$label[graph$nodes$id == parentID]

      tmp_parent <- df_expr[df_expr$cell == parent_label, ]
      tmp <- filterHM(tmp_parent[, markers], posMarker, negMarker, th)

      df_expr[rownames(tmp), "cell"] <- graph$nodes$label[i]
    }
  }

  return(df_expr$cell)
}

# Node colors: grey if no cluster carries the label, blue otherwise -----------
color_nodes_by_usage <- function(nodes, cells) {
  nodes$color <- ifelse(nodes$label %in% cells, "blue", "grey")
  nodes
}

# visNetwork widget of the tree -----------------------------------------------
# Clicking a node sets the Shiny input `input_id` to the node id.
tree_network <- function(graph, input_id, hierarchical = FALSE) {

  net <- visNetwork(graph$nodes, graph$edges, width = "100%")
  if (hierarchical) {
    net <- net %>%
      visEdges(arrows = "from") %>%
      visHierarchicalLayout()
  }
  net %>%
    visEvents(select = sprintf(
      "function(nodes) { Shiny.setInputValue('%s', nodes.nodes, {priority: 'event'}); }",
      input_id
    ))
}
