annotated_demo <- function() {
  markers <- colnames(df_expr_demoData)
  expr <- createExpressionDF(df_expr_demoData, cluster_freq_demoData, markers)
  graph <- getGraphFromLoad(nodes_demo_data, edges_demo_data)
  expr$cell <- rebuiltTree(graph, expr, prepare_thresholds(marker_th_demo_data), markers)
  list(expr = expr, graph = graph,
       counts = drop_index_column(cluster_counts_demoData),
       md = drop_index_column(meta_demo_data))
}

test_that("merged_prop_table() gives percentages per sample", {

  d <- annotated_demo()
  props <- merged_prop_table(d$counts, d$expr$cell, d$graph)

  expect_equal(unname(colSums(props)), rep(100, ncol(d$counts)))

  # labels of inner nodes get the "_remaining" suffix
  inner <- d$graph$nodes$label[vapply(d$graph$nodes$id, node_has_children, logical(1),
                                      graph_data = d$graph)]
  used_inner <- intersect(inner, d$expr$cell)
  expect_gt(length(used_inner), 0)
  expect_true(all(paste0(used_inner, "_remaining") %in% rownames(props)))
})

test_that("aggregate_counts_by_cell() rejects misaligned input", {
  expect_error(aggregate_counts_by_cell(data.frame(s1 = 1:3), c("a", "b")), "Row mismatch")
})

test_that("run_da() returns one row per phenotype and condition pair", {

  d <- annotated_demo()
  res <- run_da(d$counts, d$md, d$expr$cell, d$graph, "BH")

  expect_named(res, c("Cond1", "Cond2", "p-value", "Cell", "Naming"))
  expect_gt(nrow(res), 0)
  expect_true(all(res$`p-value` >= 0 & res$`p-value` <= 1))
  expect_true(all(res$Cell %in% d$graph$nodes$label))
})

test_that("run_da() stops on unusable metadata", {
  d <- annotated_demo()
  expect_error(run_da(d$counts, data.frame(x = 1), d$expr$cell, d$graph), "sample_id")
})

test_that("node_proportions() of the root only counts unassigned clusters", {

  d <- annotated_demo()
  props <- node_proportions(d$counts, d$expr$cell, d$graph, 1, d$md)

  pct <- to_percent(aggregate_counts_by_cell(d$counts, d$expr$cell))
  unassigned <- if ("Unassigned" %in% rownames(pct)) unname(unlist(pct["Unassigned", ])) else 0
  expect_equal(props$value, rep_len(unassigned, ncol(d$counts)))
  expect_s3_class(props$cond, "factor")
})

test_that("node_proportions() includes descendants and node_da() tests them", {

  d <- annotated_demo()
  node_id <- d$graph$nodes$id[d$graph$nodes$label == "Granulocytes"]
  labels <- node_with_descendant_labels(d$graph, node_id)
  props <- node_proportions(d$counts, d$expr$cell, d$graph, node_id, d$md)

  pct <- to_percent(aggregate_counts_by_cell(d$counts, d$expr$cell))
  expected <- colSums(pct[rownames(pct) %in% labels, , drop = FALSE])
  expect_equal(props$value, unname(expected))

  res <- node_da(props, "Granulocytes")
  expect_named(res, c("Var1", "Var2", "p-value", "Cell"))
})
