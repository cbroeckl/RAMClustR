test_that("RAMClustR rc.restore.qc.samples", {
  ramclustObj <- readRDS("testdata/rc.remove.qc.rds")
  expected <- readRDS("testdata/rc.restore.qc.samples.rds")

  actual <- rc.restore.qc.samples(ramclustObj = ramclustObj)
  
  expect_equal(actual, expected)
})