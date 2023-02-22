test_that("RAMClustR rc.expand.sample.names", {
  expected <- readRDS("testdata/rc.expand.sample.names.rds")
  ramclustObj <- readRDS("testdata/rc.get.xcms.data.rds")

  actual <- rc.expand.sample.names(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})