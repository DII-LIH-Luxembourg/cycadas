folder_with_inputs <- function() {
  dir <- tempfile("cycadas-folder")
  dir.create(dir)
  write.csv(data.frame(m1 = c(1, 6, 9, 2, 8), m2 = c(9, 1, 8, 2, 1)),
            file.path(dir, "median_expr.csv"), row.names = FALSE)
  write.csv(data.frame(cluster = 1:5, clustering_prop = rep(0.2, 5)),
            file.path(dir, "cluster_freq.csv"), row.names = FALSE)
  dir
}

test_that("node_path lists the ancestors below the root", {
  g <- initTree()
  g <- add_node(g, "Unassigned", "T cells", list("m1"), list(""), "red")
  g <- add_node(g, "T cells", "CD4 T", list("m2"), list(""), "red")
  expect_equal(node_path(g, "CD4 T"), "T cells > CD4 T")
  expect_equal(node_path(g, "T cells"), "T cells")
  expect_equal(node_path(g, "Unassigned"), "Unassigned")
})

test_that("a folder session loads the inputs, saves and reopens the workspace", {
  dir <- folder_with_inputs()
  state <- new_app_state()
  testServer(mod_folder_server, args = list(state = state, dir = dir), {
    expect_equal(nrow(state$expr), 5)
    state$graph <- add_node(state$graph, "Unassigned", "high m1", list("m1"), list(""), "red")
    state$expr$cell[c(2, 3, 5)] <- "high m1"
    session$setInputs(save = 1)
  })
  labels <- read.csv(file.path(dir, FOLDER_LABELS))
  expect_equal(labels$cluster, 1:5)
  expect_equal(labels$label, c("Unassigned", "high m1", "high m1", "Unassigned", "high m1"))
  expect_true(file.exists(file.path(dir, FOLDER_WORKSPACE)))

  state2 <- new_app_state()
  testServer(mod_folder_server, args = list(state = state2, dir = dir), {
    expect_equal(state$expr$cell[2], "high m1")
    expect_true("high m1" %in% state$graph$nodes$label)
  })
})
