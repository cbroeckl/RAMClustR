test_that("RAMClustR rc.get.csv.data", {
  filename <- file.path("testdata", "peaks.csv")
  phenoData <- file.path("testdata", "phenoData.csv")
  expected <- readRDS(file.path("testdata", "rc.get.csv.data.rds"))
  
  actual <- rc.get.csv.data(csv = filename, phenoData = phenoData, st = 5)
  
  expect_equal(actual, expected)
})
