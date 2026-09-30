test_that("createExpressionDF() adds scaled, raw, freq and cell columns", {

  markers <- c("m1", "m2")
  expr <- createExpressionDF(data.frame(a = 1:10, b = 10:1),
                             data.frame(clustering_prop = rep(0.1, 10)), markers)

  expect_named(expr, c("m1_raw", "m2_raw", "m1", "m2", "freq", "cell"))
  expect_true(all(expr$m1 >= 0 & expr$m1 <= 1))
  expect_equal(expr$m1_raw, 1:10)
  expect_true(all(expr$cell == "Unassigned"))
})

test_that("filterHM() does not depend on the order of the threshold table", {

  df <- data.frame(m1 = c(0.1, 0.9, 0.9), m2 = c(0.9, 0.1, 0.9))
  th <- data.frame(cell = c("m2", "m1"), threshold = c(0.5, 0.2))

  res <- filterHM(df, c("m1"), c("m2"), th)
  expect_equal(rownames(res), "2")

  res <- filterHM(df, c("m2", "m1"), c(), th)
  expect_equal(rownames(res), "3")
})

test_that("filterColor() marks the selected rows", {
  df <- data.frame(x = 1:4)
  expect_equal(filterColor(df, df[c(2, 4), , drop = FALSE]),
               c("other clusters", "selected phenotype", "other clusters", "selected phenotype"))
})

test_that("prepare_thresholds() colors markers by bimodality", {
  th <- prepare_thresholds(marker_th_demo_data)
  expect_null(th$X)
  expect_equal(rownames(th), th$cell)
  expect_equal(th$color == "red", th$bi_mod < BIMODALITY_CUTOFF)
})
