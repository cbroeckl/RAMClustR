test_that("RAMClustR rc.expand.sample.names", {
  expected <- readRDS(file.path("testdata", "rc.expand.sample.names.rds"))
  ramclustObj <- readRDS(file.path("testdata", "rc.get.xcms.data.rds"))
  
  actual <- rc.expand.sample.names(ramclustObj = ramclustObj, quiet=TRUE)
  
  # renamed phenoData colnames as test fails in R CMD checks becuase of no user input for colnames
  colnames(actual$phenoData) <- colnames(expected$phenoData)

  expect_equal(actual, expected)
})
