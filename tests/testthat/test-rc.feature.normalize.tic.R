test_that("RAMClustR rc.feature.normalize.quantile", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.filter.blanks.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.normalize.tic.rds"))
  
  actual <- rc.feature.normalize.tic(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})
