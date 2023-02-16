test_that("RAMClustR rc.feature.normalize.qc", {
  ramclustObj <- readRDS("testdata/rc.feature.filter.blanks.rds")
  expected <- readRDS("testdata/rc.feature.normalize.qc.rds")

  actual <- rc.feature.normalize.qc(ramclustObj = ramclustObj, qc.tag = "QC")

  expect_equal(actual, expected)
})