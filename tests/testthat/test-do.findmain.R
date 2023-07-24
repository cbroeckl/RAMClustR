test_that("RAMClustR do.findmain", {
  wd <- getwd()
  tmp <- tempdir()
  ramclustObj <- readRDS(file.path("testdata", "rc.qc.rds"))
  expected <- readRDS(file.path("testdata", "do.findmain.rds"))

  setwd(tmp)
  actual <- do.findmain(ramclustObj = ramclustObj)

  actual$history <- NA
  expected$history <- NA
  actual$params <- NA
  expected$params <- NA
  actual$phenoData <- NA
  expected$phenoData <- NA

  expect_equal(actual, expected)
  setwd(wd)
})
