test_that("RAMClustR rc.ramclustr", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.filter.cv.rds"))
  expected <- readRDS(file.path("testdata", "rc.ramclustr.rds"))
  
  actual <- rc.ramclustr(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})
