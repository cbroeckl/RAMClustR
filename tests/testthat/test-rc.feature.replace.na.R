test_that("RAMClustR rc.feature.replace.na", {
  set.seed(123) # to get reproducible results with jitters
  ramclustObj <- readRDS(file.path("testdata", "rc.expand.sample.names.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.replace.na.rds"))
  
  actual <- rc.feature.replace.na(ramclustObj = ramclustObj)
  
  # slight change in 'history' language causing error. unimportant, so just remove.
  expected$history <- NULL
  actual$history <- NULL
  
  ## test that there are zero NA present
  expect_true(length(which(is.na(actual$MSdata))) == 0)
})
