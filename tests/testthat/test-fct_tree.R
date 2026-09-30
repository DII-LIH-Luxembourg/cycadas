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
