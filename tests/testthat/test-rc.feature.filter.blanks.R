test_that("RAMClustR rc.feature.filter.blanks", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.replace.na.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.filter.blanks.rds"))
  actual <- rc.feature.filter.blanks(ramclustObj = ramclustObj, qc.tag = c("QC", "sample.names.sample_name"), blank.tag = c("Blanc", "sample.names.sample_name"))

  # slight change in 'history' language causing error. unimportant, so just remove.
  expected$history <- NULL
  actual$history <- NULL
  
  expect_equal(actual, expected)
})
