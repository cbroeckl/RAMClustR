test_that("RAMClustR rc.expand.sample.names", {
  skip_if_not_installed("xcms")
  
  load("testdata/test.rc.ramclustr.fillpeaks")
  expected <- readRDS("testdata/rc.expand.sample.names.rds")

  ramclustObj <- rc.get.xcms.data(xcmsObj = xdata)
  actual <- rc.expand.sample.names(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})