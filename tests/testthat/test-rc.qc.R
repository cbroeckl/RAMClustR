test_that("RAMClustR rc.qc", {

  tmp <- tempdir()
  ramclustObj <- readRDS(file.path("testdata", "clustering", "only_rt_rc.ramclustr.rds"))
  expected <- readRDS(file.path("testdata", "rc.qc.rds"))
  
  setwd(tmp)

  actual <- rc.qc(ramclustObj = ramclustObj, qc.tag = c("QC", "sample.names.sample_name"), out.dir = tmp)

  expect_equal(actual, expected)


})
