demo_graph <- function() getGraphFromLoad(nodes_demo_data, edges_demo_data)

test_that("add_node() appends a node and an edge to its parent", {

  g <- initTree()
  g <- add_node(g, "Unassigned", "T cells", list("CD3"), list(""), "blue")
  g <- add_node(g, "T cells", "CD4 T", list("CD4"), list("CD8"), "blue")

  expect_equal(g$nodes$label, c("Unassigned", "T cells", "CD4 T"))
  expect_equal(g$edges$to[g$edges$from == 3], 2)
  expect_equal(unlist(g$nodes$nm[3]), "CD8")
})

test_that("delete_leaf_node() only deletes leaves", {

  g <- initTree()
  g <- add_node(g, "Unassigned", "A", list("m1"), list(""), "blue")
  g <- add_node(g, "A", "B", list("m2"), list(""), "blue")

  expect_true(node_has_children(g, 2))
  expect_identical(delete_leaf_node(g, 2), g)

  g2 <- delete_leaf_node(g, 3)
  expect_equal(g2$nodes$label, c("Unassigned", "A"))
  expect_false(3 %in% g2$edges$from)
})

test_that("the root always counts as having children", {
  expect_true(node_has_children(initTree(), 1))
})

test_that("node_with_descendant_labels() collects the whole subtree", {

  g <- initTree()
  g <- add_node(g, "Unassigned", "A", list("m1"), list(""), "blue")
  g <- add_node(g, "A", "B", list("m2"), list(""), "blue")
  g <- add_node(g, "B", "C", list("m3"), list(""), "blue")
  g <- add_node(g, "Unassigned", "D", list(""), list("m1"), "blue")

  expect_setequal(node_with_descendant_labels(g, 2), c("A", "B", "C"))
  expect_setequal(node_with_descendant_labels(g, 5), "D")
})

test_that("serialize_nodes() and getGraphFromLoad() round-trip", {

  g <- demo_graph()
  g2 <- getGraphFromLoad(serialize_nodes(g$nodes), g$edges)

  expect_equal(g2$nodes$pm, g$nodes$pm)
  expect_equal(g2$nodes$nm, g$nodes$nm)
})

test_that("rebuiltTree() assigns every demo cluster to a tree label", {

  markers <- colnames(df_expr_demoData)
  expr <- createExpressionDF(df_expr_demoData, cluster_freq_demoData, markers)
  th <- prepare_thresholds(marker_th_demo_data)
  g <- demo_graph()

  cells <- rebuiltTree(g, expr, th, markers)

  expect_length(cells, nrow(expr))
  expect_true(all(cells %in% g$nodes$label))
  expect_gt(length(unique(cells)), 10)
})

test_that("rebuiltTree() matches a node created by hand", {

  markers <- c("m1", "m2")
  expr <- createExpressionDF(data.frame(m1 = c(1, 5, 9, 2), m2 = c(9, 1, 8, 2)),
                             data.frame(clustering_prop = rep(0.25, 4)), markers)
  th <- data.frame(cell = markers, threshold = c(0.5, 0.5))

  g <- add_node(initTree(), "Unassigned", "m1+", list("m1"), list(""), "blue")
  g <- add_node(g, "m1+", "m1+m2-", list(""), list("m2"), "blue")

  expect_equal(rebuiltTree(g, expr, th, markers),
               c("Unassigned", "m1+m2-", "m1+", "Unassigned"))
})

test_that("node_lineage() and lineage_markers() follow the path to the root", {

  g <- initTree()
  g <- add_node(g, "Unassigned", "A", list("m1"), list(""), "blue")
  g <- add_node(g, "A", "B", list(""), list("m2"), "blue")
  g <- add_node(g, "B", "C", list("m3"), list("m4"), "blue")

  expect_equal(node_lineage(g, 4), c(1, 2, 3, 4))
  expect_equal(node_lineage(g, 1), 1)
  expect_equal(lineage_markers(g, 4), list(pos = c("m1", "m3"), neg = c("m2", "m4")))
  expect_equal(lineage_markers(g, 1), list(pos = character(0), neg = character(0)))
})

test_that("tree_network() draws edges from parent to child without the root loop", {

  g <- add_node(initTree(), "Unassigned", "A", list("m1"), list(""), "blue")
  expr <- data.frame(cell = c("A", "A", "Unassigned"), freq = c(10, 20, 70))
  net <- tree_network(g, "tree_click", expr = expr, selected = 2)

  expect_equal(net$x$edges$from, 1)
  expect_equal(net$x$edges$to, 2)
  expect_equal(net$x$nodes$borderWidth, c(1, 3))
  expect_match(net$x$nodes$title[2], "30.00% of cells (2 clusters)", fixed = TRUE)
})
