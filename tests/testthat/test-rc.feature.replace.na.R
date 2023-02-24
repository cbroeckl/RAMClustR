test_that("RAMClustR rc.feature.replace.na", {
  ramclustObj <- readRDS(file.path("testdata", "rc.expand.sample.names.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.replace.na.rds"))
  
  actual <- rc.feature.replace.na(ramclustObj = ramclustObj)

  expect_equal(actual, expected, tolerance = 0.02)
})
