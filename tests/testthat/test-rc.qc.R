test_that("RAMClustR rc.qc", {
  wd <- getwd()
  tmp <- tempdir()
  ramclustObj <- readRDS("testdata/rc.ramclustr.rds")
  expected <- readRDS("testdata/rc.qc.rds")

  setwd(tmp)

  actual <- rc.qc(ramclustObj = ramclustObj)

  expect_equal(actual, expected)

  setwd(wd)
})