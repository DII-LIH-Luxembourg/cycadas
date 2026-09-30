test_that("kmeansTH() stores the estimate next to the threshold", {
  df <- data.frame(m1 = c(rnorm(50, 0.2, 0.05), rnorm(50, 0.8, 0.05)))
  th <- kmeansTH(df)
  expect_equal(th$estimated, th$threshold)
  expect_false(any(threshold_modified(th)))
})

test_that("add_estimates() fills in missing estimates only", {

  markers <- colnames(df_expr_demoData)[1:3]
  expr <- createExpressionDF(df_expr_demoData[, markers], cluster_freq_demoData, markers)
  th <- prepare_thresholds(marker_th_demo_data[marker_th_demo_data$cell %in% markers, ])
  expect_true(all(is.na(th$estimated)))

  th$estimated[1] <- 0.123
  th <- add_estimates(th, expr, markers)

  expect_equal(th$estimated[1], 0.123)
  expected <- kmeansTH(expr[, markers[2:3]])
  expect_equal(th$estimated[2:3], expected$threshold[match(th$cell[2:3], expected$cell)])
})

test_that("reset_thresholds() restores estimates of the given markers", {

  th <- data.frame(cell = c("a", "b", "c"), threshold = c(0.5, 0.7, 0.2),
                   estimated = c(0.5, 0.4, NA))
  expect_equal(threshold_modified(th), c(FALSE, TRUE, FALSE))

  expect_equal(reset_thresholds(th, "b")$threshold, c(0.5, 0.4, 0.2))
  # no estimate: the threshold is kept
  expect_equal(reset_thresholds(th)$threshold, c(0.5, 0.4, 0.2))
})

test_that("write_thresholds() output can be uploaded again", {

  th <- data.frame(cell = c("a", "b"), threshold = c(0.5, 0.7), estimated = c(0.5, 0.4),
                   bi_mod = c(0.7, 0.4), color = c("blue", "red"))
  f <- tempfile(fileext = ".csv")
  write_thresholds(th, f)
  back <- prepare_thresholds(read.csv(f))

  expect_equal(back$threshold, th$threshold)
  expect_equal(back$estimated, th$estimated)
  expect_equal(back$color, th$color)
})
