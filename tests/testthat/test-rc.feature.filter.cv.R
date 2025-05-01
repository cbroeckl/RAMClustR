test_that("RAMClustR rc.feature.filter.cv", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.normalize.qc.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.filter.cv.rds"))
  
  actual <- rc.feature.filter.cv(ramclustObj = ramclustObj, qc.tag = c("QC", "sample.names.sample_name"))
  
  ## remove history and params to keep them from triggering failure
  expected$history <- NULL
  expected$params <- NULL
  actual$history <- NULL
  actual$params <- NULL

  expect_equal(actual, expected)
})
