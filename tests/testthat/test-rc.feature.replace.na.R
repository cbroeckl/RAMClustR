test_that("RAMClustR rc.feature.replace.na", {
  ramclustObj <- readRDS("testdata/rc.expand.sample.names.rds")
  expected <- readRDS("testdata/rc.feature.replace.na.rds")

  actual <- rc.feature.replace.na(ramclustObj = ramclustObj)

  expect_equal(actual, expected, tolerance = 0.02)
})