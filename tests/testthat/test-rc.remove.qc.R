test_that("RAMClustR rc.remove.qc", {
  ramclustObj <- readRDS("testdata/rc.get.xcms.data.rds")
  expected <- readRDS("testdata/rc.remove.qc.rds")

  actual <- rc.remove.qc(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})