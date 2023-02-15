test_that("RAMClustR do.findmain", {
  wd <- getwd()
  tmp <- tempdir()
  ramclustObj <- readRDS("testdata/rc.qc.rds")
  expected <- readRDS("testdata/do.findmain.rds")

  setwd(tmp)
  actual <- do.findmain(ramclustObj = ramclustObj)

  expect_equal(actual, expected, tolerance = 0.02)
  setwd(wd)
})