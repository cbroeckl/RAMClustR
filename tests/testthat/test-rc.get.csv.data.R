test_that("RAMClustR rc.get.csv.data", {
  wd <- getwd()
  filename <- file.path(wd, "testdata/peaks.csv")
  phenoData <- file.path(wd, "testdata/phenoData.csv")
  expected <- readRDS("testdata/rc.get.csv.data.rds")

  actual <- rc.get.csv.data(csv = filename, phenoData = phenoData, st = 5)
  
  expect_equal(actual, expected)
})