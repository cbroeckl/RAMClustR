test_that("RAMClustR works", {
  wd <- getwd()
  tmp <- tempdir()
  load("testdata/test.fillpeaks")
  expected <- readLines("testdata/output.msp")

  setwd(tmp)
  ramclustr_obj <- ramclustR(xcmsObj = xdata)
  write.msp(ramclustr_obj, one.file = TRUE)
  diff <- setdiff(expected, readLines("spectra/fill.msp"))
  expect_true(length(diff) < 10)
  setwd(wd)
})