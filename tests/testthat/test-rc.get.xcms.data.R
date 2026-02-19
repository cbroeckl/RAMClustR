test_that("RAMClustR rc.get.xcms.data", {
  skip_if_not_installed("xcms")
  
  load(file.path("testdata", "test.fillpeaks"))
  expected <- readRDS(file.path("testdata", "rc.get.xcms.data.rds"))
  
  actual <- rc.get.xcms.data(xcmsObj = xdata)

  actual$history <- NA
  expected$history <- NA

  expect_equal(actual, expected)
})
