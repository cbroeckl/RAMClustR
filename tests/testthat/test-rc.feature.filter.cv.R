test_that("RAMClustR rc.feature.filter.cv", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.normalize.qc.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.filter.cv.rds"))
  
  actual <- rc.feature.filter.cv(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})
