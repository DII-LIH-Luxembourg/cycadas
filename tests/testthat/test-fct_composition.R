test_that("freq_percent() scales fractions and percentages to percent", {
  expect_equal(freq_percent(c(0.25, 0.75)), c(25, 75))
  expect_equal(freq_percent(c(25, 75)), c(25, 75))
})

small_tree <- function() {
  g <- initTree()
  g <- add_node(g, "Unassigned", "T", list("CD3"), list(""), "blue")
  g <- add_node(g, "T", "CD4 T", list("CD4"), list(""), "blue")
  g <- add_node(g, "T", "CD8 T", list("CD8"), list(""), "blue")
  g <- add_node(g, "Unassigned", "B", list("CD19"), list(""), "blue")
  g
}

test_that("phenotype_composition() gives total, remaining and share of parent", {

  g <- small_tree()
  expr <- data.frame(cell = c("Unassigned", "T", "CD4 T", "CD4 T", "CD8 T", "B"),
                     freq = c(10, 5, 30, 20, 25, 10))
  comp <- phenotype_composition(g, expr)

  # depth-first order, children in creation order
  expect_equal(comp$phenotype, c("Unassigned", "T", "CD4 T", "CD8 T", "B"))
  expect_equal(comp$depth, c(0, 1, 2, 2, 1))
  expect_equal(comp$total, c(100, 80, 50, 25, 10))
  expect_equal(comp$remaining, c(10, 5, 50, 25, 10))
  expect_equal(comp$of_parent, c(NA, 80, 62.5, 31.25, 10))
  expect_equal(comp$clusters_total, c(6, 4, 2, 1, 1))
  expect_equal(comp$parent, c(NA, "Unassigned", "T", "T", "Unassigned"))
})

test_that("phenotype_composition() handles empty phenotypes", {
  g <- small_tree()
  expr <- data.frame(cell = c("Unassigned", "T"), freq = c(40, 60))
  comp <- phenotype_composition(g, expr)
  expect_equal(comp$total, c(100, 60, 0, 0, 0))
  expect_equal(comp$of_parent[3], 0)
})

test_that("phenotype_composition() adds up on the annotated demo", {

  markers <- colnames(df_expr_demoData)
  expr <- createExpressionDF(df_expr_demoData, cluster_freq_demoData, markers)
  g <- getGraphFromLoad(nodes_demo_data, edges_demo_data)
  expr$cell <- rebuiltTree(g, expr, prepare_thresholds(marker_th_demo_data), markers)
  comp <- phenotype_composition(g, expr)

  expect_equal(nrow(comp), nrow(g$nodes))
  expect_equal(sum(comp$remaining), 100)
  # a parent's total is its remaining share plus its children's totals
  for (i in seq_len(nrow(comp))) {
    kids <- comp$parent %in% comp$phenotype[i]
    expect_equal(comp$total[i], comp$remaining[i] + sum(comp$total[kids]))
  }
})

test_that("write_composition() writes percentages rounded to 3 digits", {
  g <- small_tree()
  expr <- data.frame(cell = c("Unassigned", "T"), freq = c(100 / 3, 200 / 3))
  f <- tempfile(fileext = ".csv")
  write_composition(phenotype_composition(g, expr), f)
  back <- read.csv(f)
  expect_named(back, c("phenotype", "parent", "depth", "clusters_total", "clusters_remaining",
                       "pct_total", "pct_remaining", "pct_of_parent"))
  expect_equal(back$pct_total[2], 66.667)
})

test_that("CATALYST-style fractions are shown as percent", {
  expr <- createExpressionDF(data.frame(a = 1:4), data.frame(clustering_prop = rep(0.25, 4)), "m")
  expect_equal(expr$freq, rep(25, 4))
})
