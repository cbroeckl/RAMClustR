test_that("RAMClustR rc.ramclustr", {
  ramclustObj <- readRDS("testdata/rc.feature.filter.cv.rds")
  expected <- readRDS("testdata/rc.ramclustr.rds")

  actual <- rc.ramclustr(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})