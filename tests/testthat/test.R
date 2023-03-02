test_that("RAMClustR with xcms works", {
  skip_if_not_installed("xcms")
  set.seed(123) # to get reproducible results with jitters
  wd <- getwd()
  tmp <- tempdir()
  load("testdata/test.fillpeaks")
  expected <- readLines("testdata/output.msp")
  setwd(tmp)

  ramclustr_obj <- ramclustR(xcmsObj = xdata, maxt = 20, sr = 0.5, mzdec = 4)
  write.msp(ramclustr_obj, one.file = TRUE)
  mismatches <- setdiff(expected, readLines("spectra/fill.msp"))

  expect_true(length(mismatches) < 10)
  setwd(wd)
})

test_that("RAMClustR with csv works", {
  skip("disabled")
  wd <- getwd()
  tmp <- tempdir()
  filename <- file.path(wd, "testdata/peaks.csv")
  expected <- readRDS("testdata/test_csv.rds")
  pheno <- file.path(wd, "testdata/phenoData.csv")

  setwd(tmp)
  actual <- ramclustR(
    ms = filename,
    pheno_csv = pheno,
    st = 5,
    maxt = 1,
    blocksize = 1000
  )

  actual$history <- NA
  expected$history <- NA

  expect_equal_labels(actual$labels, expected$labels)
  expect_equal_MSdata(actual$MSdata, expected$MSdata)

  actual$labels <- expected$labels <- actual$MSdata <- expected$MSdata <- NA

  expect_equal(actual, expected)
  setwd(wd)
})
