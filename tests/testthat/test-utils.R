test_that("availMarkers() is empty if input is equla to colnames", {

  df <- data.frame (m1  = c("v1", "v2"),
                    m2 = c("v3", "v4")
  )
  expect_true( is_empty(availMarkers(c("m1","m2"))) )
})

gen_DF <- function() {

  y <- rnorm(100, mean = 2.0, sd = 1)
  z <- rnorm(100, mean = 8.0, sd = 1)

  a <- rnorm(100, mean = 2.5, sd = 0.5)
  b <- rnorm(100, mean = 6.0, sd = 0.5)

  df <- data.frame(m1 = c(z,y), m2 = c(a,b))

  return(df)

}

gen_TH <- function() {

  df <- data.frame(cell = c("m1", "m2") ,value = c(5.0, 2.0))

  return(df)
}

test_that("kmean_TH() selcted threshold is within middle range", {

  mydf <- gen_DF()
  th_tmp <- kmeansTH(mydf)

  expect_equal(th_tmp$value[1], (mean(mydf$m1)+mean(mydf$m2))/2 , tolerance = 0.5)
})

test_that("normalize01() works as expected", {

  mydf <- gen_DF()

  df01 <- normalize01(mydf)

  expect_true(min(df01$m1) >= 0 & max(df01$m1) <= 1)

})

test_that("filterHM() returns full ddf if pos and neg are empty", {

  mydf <- gen_DF()
  posMarker <- c()
  negMarker <- c()

  th <- c(1.2, 3.4)

  df_filtered <- filterHM(mydf, posMarker, negMarker, th)
  expect_setequal(mydf$m1, df_filtered$m1)
})

test_that("filterHM() returns correct HM of neg marker", {

  mydf <- gen_DF()
  posMarker <- c()
  negMarker <- c("m2")

  m2_below_2 <- mydf[mydf$m2 < 2.0, ]

  th <- gen_TH()

  df_filtered <- filterHM(mydf, posMarker, negMarker, th)

  # test for marker m2 and th < 2.0
  expect_setequal(m2_below_2$m2, df_filtered$m2)

})

test_that("filterHM() returns correct HM of pos marker", {

  mydf <- gen_DF()
  posMarker <- c("m1")
  negMarker <- c()

  m1_above_2 <- mydf[mydf$m1 > 5.0, ]

  th <- gen_TH()

  df_filtered <- filterHM(mydf, posMarker, negMarker, th)

  # test for marker m2 and th < 2.0
  expect_setequal(m1_above_2$m1, df_filtered$m1)

})

test_that("filterHM() returns correct HM of pos and neg marker", {

  mydf <- gen_DF()
  posMarker <- c("m1")
  negMarker <- c("m2")

  df_intersect <- mydf[mydf$m1 > 5.0 & mydf$m2 < 2.0, ]

  th <- gen_TH()

  df_filtered <- filterHM(mydf, posMarker, negMarker, th)

  # test for marker m2 and th < 2.0
  # expect_setequal(df_intersect$m1, df_filtered$m1)
  # expect_setequal(df_intersect$m2, df_filtered$m2)
  expect_equal(df_intersect$m1, df_filtered$m1)

})

test_that("setPhenotypeName set the name for pos and neg markers", {

  m <- c("CD3","CD4")
  ph_name <- ""

  ph_name <- setPhenotypeName(m, "pos", ph_name)
  # split by pos sign
  tmp <- strsplit(ph_name, "\\+")
  # check if both lists a equal
  expect_equal(unlist(tmp), m)
})










