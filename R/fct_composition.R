# Phenotype composition -------------------------------------------------------
# Approximate share of cells per phenotype, from the cluster frequencies of
# the expression table (expr$freq, in percent of all cells).

# Cluster frequencies as percent of all cells, whatever the input scale
# (the CSV import gives percent, CATALYST and RemoteSOM give fractions)
freq_percent <- function(freq) {
  if (is.null(freq)) return(NULL)
  total <- sum(freq, na.rm = TRUE)
  if (total > 0) freq / total * 100 else freq
}

# Node ids in depth-first order, children in creation order
tree_order <- function(graph, node_id = 1, depth = 0) {
  out <- data.frame(id = node_id, depth = depth)
  for (child in node_children(graph, node_id)) {
    out <- rbind(out, tree_order(graph, child, depth + 1))
  }
  out
}

# One row per phenotype in tree order:
#   total      % of all cells in the phenotype including its subpopulations
#   remaining  % of all cells assigned to the phenotype itself
#   of_parent  total as % of the parent's total (NA for the root)
phenotype_composition <- function(graph, expr) {

  ord <- tree_order(graph)
  nodes <- graph$nodes[match(ord$id, graph$nodes$id), ]
  own_freq <- vapply(nodes$label, function(l) sum(expr$freq[expr$cell == l]), numeric(1))
  own_n <- vapply(nodes$label, function(l) sum(expr$cell == l), integer(1))

  subtree <- lapply(ord$id, function(id) node_with_descendant_labels(graph, id))
  total <- vapply(seq_along(subtree), function(i) {
    if (ord$id[i] == 1) sum(expr$freq) else sum(own_freq[nodes$label %in% subtree[[i]]])
  }, numeric(1))
  total_n <- vapply(seq_along(subtree), function(i) {
    if (ord$id[i] == 1) nrow(expr) else sum(own_n[nodes$label %in% subtree[[i]]])
  }, numeric(1))

  parent_id <- graph$edges$to[match(ord$id, graph$edges$from)]
  parent_total <- total[match(parent_id, ord$id)]
  of_parent <- ifelse(ord$id == 1 | is.na(parent_total) | parent_total == 0, NA,
                      total / parent_total * 100)

  data.frame(
    id = ord$id,
    phenotype = nodes$label,
    parent = ifelse(ord$id == 1, NA, graph$nodes$label[match(parent_id, graph$nodes$id)]),
    depth = ord$depth,
    clusters_total = total_n,
    clusters_remaining = unname(own_n),
    total = total,
    remaining = unname(own_freq),
    of_parent = of_parent,
    stringsAsFactors = FALSE
  )
}

# Composition table as written by the download
write_composition <- function(comp, file) {
  out <- comp[, c("phenotype", "parent", "depth", "clusters_total", "clusters_remaining",
                  "total", "remaining", "of_parent")]
  names(out)[6:8] <- c("pct_total", "pct_remaining", "pct_of_parent")
  out[6:8] <- lapply(out[6:8], round, 3)
  write.csv(out, file, row.names = FALSE)
}
