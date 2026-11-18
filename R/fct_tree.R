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

# Annotation order ------------------------------------------------------------
# rebuiltTree() handles the nodes in row order. A parent therefore has to come
# before its children, and among siblings the earlier row takes the clusters
# that match several siblings. The functions below change that order or the
# parent of a node; the caller re-assigns the clusters with rebuiltTree().

node_parent <- function(graph, node_id) {
  graph$edges$to[match(node_id, graph$edges$from)]
}

# Children of a node in annotation order
ordered_children <- function(graph, node_id) {
  kids <- node_children(graph, node_id)
  graph$nodes$id[graph$nodes$id %in% kids]
}

# Nodes and edges in depth-first order from the root, siblings in their
# current row order
sort_tree <- function(graph) {
  visit <- function(id) c(id, unlist(lapply(ordered_children(graph, id), visit)))
  ids <- visit(1)
  graph$nodes <- graph$nodes[match(ids, graph$nodes$id), ]
  graph$edges <- graph$edges[match(ids, graph$edges$from), ]
  rownames(graph$edges) <- NULL
  graph
}

# Swap a node with its previous (by = -1) or next (by = 1) sibling. Returns
# the graph unchanged if there is no sibling in that direction.
move_sibling <- function(graph, node_id, by) {
  if (node_id == 1) return(graph)
  sibs <- ordered_children(graph, node_parent(graph, node_id))
  i <- match(node_id, sibs)
  j <- i + by
  if (is.na(i) || j < 1 || j > length(sibs)) return(graph)
  sibs[c(i, j)] <- sibs[c(j, i)]
  rows <- seq_len(nrow(graph$nodes))
  rows[graph$nodes$id %in% sibs] <- match(sibs, graph$nodes$id)
  graph$nodes <- graph$nodes[rows, ]
  sort_tree(graph)
}

# Nodes that can become the parent of `node_id`: all except the node itself,
# its descendants and its current parent
possible_parents <- function(graph, node_id) {
  excluded <- c(node_id, all_my_children(graph, node_id), node_parent(graph, node_id))
  graph$nodes$id[!graph$nodes$id %in% excluded]
}

# Attach a node, with its subtree, to another parent. It becomes the last
# child there. The node keeps its own markers; its lineage changes.
move_node <- function(graph, node_id, new_parent_id) {
  if (node_id == 1 || !new_parent_id %in% possible_parents(graph, node_id)) {
    stop("Cannot move this node below the selected parent.")
  }
  graph$edges$to[graph$edges$from == node_id] <- new_parent_id
  row <- match(node_id, graph$nodes$id)
  graph$nodes <- graph$nodes[c(setdiff(seq_len(nrow(graph$nodes)), row), row), ]
  sort_tree(graph)
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

# Ancestors of a node, from the root down to the node itself ----------------
node_lineage <- function(graph_data, node_id) {
  ids <- node_id
  repeat {
    parent <- graph_data$edges$to[match(ids[1], graph_data$edges$from)]
    if (is.na(parent) || parent %in% ids) break
    ids <- c(parent, ids)
  }
  ids
}

# Markers that define a node, collected along its lineage ---------------------
lineage_markers <- function(graph_data, node_id) {
  rows <- match(node_lineage(graph_data, node_id), graph_data$nodes$id)
  collect <- function(col) {
    m <- unlist(graph_data$nodes[[col]][rows])
    unique(m[!is.na(m) & nzchar(m)])
  }
  list(pos = collect("pm"), neg = collect("nm"))
}

# Markers required both positive and negative in the lineage of a node or of
# one of its descendants. No cluster can match such a node.
lineage_conflicts <- function(graph_data, node_id) {
  ids <- c(node_id, all_my_children(graph_data, node_id))
  unique(unlist(lapply(ids, function(id) {
    m <- lineage_markers(graph_data, id)
    intersect(m$pos, m$neg)
  })))
}

# Node colors: grey if no cluster carries the label, blue otherwise -----------
TREE_COLOR_USED <- "#2c7be5"
TREE_COLOR_EMPTY <- "#adb5bd"

color_nodes_by_usage <- function(nodes, cells) {
  nodes$color <- ifelse(nodes$label %in% cells, TREE_COLOR_USED, TREE_COLOR_EMPTY)
  nodes
}

# Hover text: own marker definition, cluster count and cell share -----------
node_tooltips <- function(graph, expr) {
  comp <- phenotype_composition(graph, expr)
  comp <- comp[match(graph$nodes$id, comp$id), ]
  nodes <- graph$nodes
  vapply(seq_len(nrow(nodes)), function(i) {
    pm <- unlist(nodes$pm[i]); nm <- unlist(nodes$nm[i])
    def <- c(if (length(pm)) paste0(pm[nzchar(pm)], "+"), if (length(nm)) paste0(nm[nzchar(nm)], "\u2212"))
    paste0("<b>", htmltools::htmlEscape(nodes$label[i]), "</b><br>",
           if (length(def)) paste0(paste(def, collapse = " "), "<br>"),
           sprintf("%.2f%% of cells (%d clusters)<br>remaining %.2f%% (%d clusters)",
                   comp$total[i], comp$clusters_total[i], comp$remaining[i], comp$clusters_remaining[i]),
           if (!is.na(comp$of_parent[i])) sprintf("<br>%.1f%% of %s", comp$of_parent[i],
                                                  htmltools::htmlEscape(comp$parent[i])))
  }, character(1))
}

# Node styles: fill by usage, a dark ring around the selected node -----------
tree_selection_style <- function(nodes, expr, selected = NULL) {
  fill <- if (is.null(expr)) rep(TREE_COLOR_USED, nrow(nodes)) else
    color_nodes_by_usage(nodes, expr$cell)$color
  is_sel <- nodes$id %in% selected
  data.frame(id = nodes$id,
             color.background = fill,
             color.border = ifelse(is_sel, "#0b1f3a", fill),
             color.highlight.background = fill,
             color.highlight.border = "#0b1f3a",
             font.color = ifelse(fill == TREE_COLOR_EMPTY, "#343a40", "#ffffff"),
             borderWidth = ifelse(is_sel, 3, 1))
}

# visNetwork widget of the tree -----------------------------------------------
# Clicking a node sets the Shiny input `input_id` to the node id. With `expr`
# the nodes are colored by usage and get hover details.
tree_network <- function(graph, input_id, expr = NULL, selected = NULL, height = NULL) {

  nodes <- data.frame(tree_selection_style(graph$nodes, expr, selected),
                      label = graph$nodes$label, check.names = FALSE)
  if (!is.null(expr)) nodes$title <- node_tooltips(graph, expr)

  # edges point from parent to child for the layout; the root's self-loop is
  # not drawn
  edges <- graph$edges[graph$edges$from != graph$edges$to, ]
  edges <- data.frame(from = edges$to, to = edges$from)

  visNetwork(nodes, edges, width = "100%", height = height) %>%
    visNodes(shape = "box", margin = list(top = 6, bottom = 6, left = 10, right = 10),
             shapeProperties = list(borderRadius = 6),
             font = list(size = 15, face = "system-ui, sans-serif")) %>%
    visEdges(color = list(color = "#ced4da", highlight = TREE_COLOR_USED), width = 1.5,
             smooth = list(type = "cubicBezier", forceDirection = "horizontal", roundness = 0.5)) %>%
    visHierarchicalLayout(direction = "LR", sortMethod = "directed", shakeTowards = "roots",
                          levelSeparation = 230, nodeSpacing = 42) %>%
    visPhysics(enabled = FALSE) %>%
    visInteraction(hover = TRUE, tooltipDelay = 150,
                   keyboard = list(enabled = TRUE, bindToWindow = FALSE)) %>%
    visEvents(
      select = sprintf(
        "function(nodes) { Shiny.setInputValue('%s', nodes.nodes, {priority: 'event'}); }",
        input_id),
      # first view: see cycadasTree.readable() in inst/app/www/cycadas.js
      afterDrawing = "function() {
        if (this.cycadasFitted) return;
        this.cycadasFitted = true;
        if (window.cycadasTree) window.cycadasTree.readable(this);
      }"
    )
}
