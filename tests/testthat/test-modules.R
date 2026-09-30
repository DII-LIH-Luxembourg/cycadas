# Module tests drive the servers with shiny::testServer() and check the
# shared state. The state is filled directly to keep the tests fast.

small_state <- function() {
  markers <- c("m1", "m2")
  # scaled to 0-1, m1 is above 0.5 in rows 2, 3, 5 and m2 below 0.5 in rows 2, 4, 5
  cell_freq <- data.frame(cluster = 1:5, clustering_prop = rep(0.2, 5))
  expr <- createExpressionDF(data.frame(m1 = c(1, 6, 9, 2, 8), m2 = c(9, 1, 8, 2, 1)),
                             cell_freq, markers)
  state <- new_app_state()
  isolate({
    state$markers <- markers
    state$cell_freq <- cell_freq
    state$expr <- expr
    state$th <- data.frame(cell = markers, threshold = c(0.5, 0.5), color = "blue",
                           bi_mod = 0.6, row.names = markers)
    state$graph <- initTree()
  })
  state
}

test_that("annotation: creating a node labels the selected clusters", {

  state <- small_state()
  testServer(mod_annotation_server, args = list(state = state), {
    session$setInputs(parentPicker = "Unassigned", treePickerPos = "m1", treePickerNeg = "m2")
    expect_equal(rownames(preview()), c("2", "5"))

    session$setInputs(newNode = "m1+m2-", createNodeBtn = 1)

    expect_equal(state$graph$nodes$label, c("Unassigned", "m1+m2-"))
    expect_equal(state$expr$cell,
                 c("Unassigned", "m1+m2-", "Unassigned", "Unassigned", "m1+m2-"))
    expect_equal(selected_label(), "m1+m2-")
  })
})

test_that("annotation: duplicate names and empty selections are rejected", {

  state <- small_state()
  testServer(mod_annotation_server, args = list(state = state), {
    session$setInputs(parentPicker = "Unassigned", treePickerPos = "m1",
                      newNode = "Unassigned", createNodeBtn = 1)
    expect_equal(nrow(state$graph$nodes), 1)

    session$setInputs(treePickerPos = NULL, treePickerNeg = NULL,
                      newNode = "x", createNodeBtn = 2)
    expect_equal(nrow(state$graph$nodes), 1)
  })
})

test_that("annotation: deleting a leaf gives its clusters back to the parent", {

  state <- small_state()
  testServer(mod_annotation_server, args = list(state = state), {
    session$setInputs(parentPicker = "Unassigned", treePickerPos = "m1",
                      newNode = "m1+", createNodeBtn = 1)
    expect_equal(sum(state$expr$cell == "m1+"), 3)

    session$setInputs(parentPicker = "m1+", deleteNodeBtn = 1)
    expect_equal(state$graph$nodes$label, "Unassigned")
    expect_true(all(state$expr$cell == "Unassigned"))
    expect_equal(selected_label(), "Unassigned")
  })
})

test_that("annotation: a tree click selects the node", {

  state <- small_state()
  isolate(state$graph <- add_node(state$graph, "Unassigned", "A", list("m1"), list(""), "blue"))
  testServer(mod_annotation_server, args = list(state = state), {
    session$setInputs(tree_click = 2)
    expect_equal(selected_label(), "A")
    expect_equal(current_node()$id, 2)
  })
})

test_that("thresholds: clicking the plot moves the threshold and re-assigns clusters", {

  state <- small_state()
  isolate({
    state$graph <- add_node(state$graph, "Unassigned", "m1+", list("m1"), list(""), "blue")
    rebuild_annotation(state)
  })
  expect_equal(sum(isolate(state$expr$cell) == "m1+"), 3)

  testServer(mod_thresholds_server, args = list(state = state), {
    session$setInputs(table_rows_selected = 1)
    expect_equal(selected()$cell, "m1")

    session$setInputs(scatter_click = list(x = 0.95))
    expect_equal(state$th["m1", "threshold"], 0.95)
    expect_equal(sum(state$expr$cell == "m1+"), 1)
  })
})

test_that("workspace: clearing resets the shared state", {

  state <- small_state()
  testServer(mod_workspace_server, args = list(state = state), {
    session$setInputs(btnClearWorkspace = 1)
    expect_null(state$expr)
    expect_null(state$graph)
    expect_equal(state$reset_version, 1)
  })
})

test_that("workspace: a saved workspace restores the annotation", {

  state <- small_state()
  isolate({
    state$graph <- add_node(state$graph, "Unassigned", "m1+", list("m1"), list(""), "blue")
    rebuild_annotation(state)
  })
  path <- tempfile(fileext = ".rds")
  isolate(save_workspace(path, workspace_from_state(state)))

  restored <- new_app_state()
  isolate({
    restore_workspace(restored, load_workspace(path))
    expect_equal(restored$expr, state$expr)
    expect_equal(restored$markers, state$markers)
    expect_equal(restored$graph$nodes$label, c("Unassigned", "m1+"))
    expect_equal(unlist(restored$graph$nodes$pm[2]), "m1")
  })
})

test_that("DA: results are computed from the shared state", {

  markers <- colnames(df_expr_demoData)
  state <- new_app_state()
  isolate({
    state$markers <- markers
    state$expr <- createExpressionDF(df_expr_demoData, cluster_freq_demoData, markers)
    state$th <- prepare_thresholds(marker_th_demo_data)
    state$graph <- getGraphFromLoad(nodes_demo_data, edges_demo_data)
    state$counts <- drop_index_column(cluster_counts_demoData)
    state$md <- drop_index_column(meta_demo_data)
    rebuild_annotation(state)
  })

  testServer(mod_da_server, args = list(state = state), {
    session$setInputs(correction_method = "BH", doDA = 1)
    expect_gt(nrow(result()), 0)
  })

  testServer(mod_da_tree_server, args = list(state = state), {
    session$setInputs(tree_click = 2)
    expect_equal(node_label(), "Granulocytes")
    expect_equal(nrow(props()), ncol(state$counts))
  })
})
