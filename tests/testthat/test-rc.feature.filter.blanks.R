test_that("RAMClustR rc.feature.filter.blanks", {
  ramclustObj <- readRDS("testdata/rc.feature.replace.na.rds")
  expected <- readRDS("testdata/rc.feature.filter.blanks.rds")

  actual <- rc.feature.filter.blanks(ramclustObj = ramclustObj, blank.tag = "Blanc")

  expect_equal(actual, expected)
})