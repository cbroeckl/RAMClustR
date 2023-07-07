test_that("RAMClustR rc.feature.normalize.quantile", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.filter.blanks.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.normalize.quantile.rds"))
  
  actual <- rc.feature.normalize.quantile(ramclustObj = ramclustObj)
  expect_equal(actual, expected)
})
