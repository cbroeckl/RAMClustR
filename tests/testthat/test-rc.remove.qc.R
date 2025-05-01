test_that("RAMClustR rc.remove.qc", {
  ramclustObj <- readRDS(file.path("testdata", "rc.get.xcms.data.rds"))
  expected <- readRDS(file.path("testdata", "rc.remove.qc.rds"))
  
  actual <- rc.remove.qc(ramclustObj = ramclustObj, qc.tag = c("QC", "sample.names.sample_name"))

  ## remove history and params to keep them from triggering failure
  expected$history <- NULL
  expected$params <- NULL
  actual$history <- NULL
  actual$params <- NULL
  
  expect_equal(actual, expected)
})
