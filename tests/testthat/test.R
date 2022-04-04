test_that("RAMClustR with xcms works", {
  skip_if_not_installed("xcms")
  wd <- getwd()
  tmp <- tempdir()
  load("testdata/test.fillpeaks")
  expected <- readLines("testdata/output.msp")

  setwd(tmp)
  ramclustr_obj <- ramclustR(xcmsObj = xdata)
  write.msp(ramclustr_obj, one.file = TRUE)
  diff <- setdiff(expected, readLines("spectra/fill.msp"))
  expect_true(length(diff) < 10)
  setwd(wd)
})

test_that("RAMClustR with csv works", {
  wd <- getwd()
  tmp <- tempdir()
  filename <- file.path(wd, "testdata/peaks.csv")
  expected <- readRDS("testdata/ramclustObj.rds")

  setwd(tmp)
  actual <- ramclustR(ms = filename, st = 5, maxt = 1, blocksize = 1000)
  actual$history <- NA
  expected$history <- NA

  expect_equal_labels(actual$labels, expected$labels)
  expect_equal_MSdata(actual$MSdata, expected$MSdata)

  actual$labels <- expected$labels <- actual$MSdata <- expected$MSdata <- NA

  expect_equal(actual, expected)
  setwd(wd)
})
