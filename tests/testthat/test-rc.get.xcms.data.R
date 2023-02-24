test_that("RAMClustR rc.get.xcms.data", {
  skip_if_not_installed("xcms")
  
  load(file.path("testdata", "test.rc.ramclustr.fillpeaks"))
  expected <- readRDS(file.path("testdata", "rc.get.xcms.data.rds"))
  
  actual <- rc.get.xcms.data(xcmsObj = xdata)

  expect_equal(actual, expected)
})
