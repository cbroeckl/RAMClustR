test_that("RAMClustR rc.get.xcms.data", {
  skip_if_not_installed("xcms")
  
  load("testdata/test.rc.ramclustr.fillpeaks")
  expected <- readRDS("testdata/rc.get.xcms.data.rds")

  actual <- rc.get.xcms.data(xcmsObj = xdata)

  expect_equal(actual, expected)
})