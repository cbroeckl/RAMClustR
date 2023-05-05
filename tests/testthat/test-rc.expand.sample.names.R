test_that("RAMClustR rc.expand.sample.names", {
  skip("disabled")
  expected <- readRDS(file.path("testdata", "rc.expand.sample.names.rds"))
  ramclustObj <- readRDS(file.path("testdata", "rc.get.xcms.data.rds"))
  
  actual <- rc.expand.sample.names(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})
