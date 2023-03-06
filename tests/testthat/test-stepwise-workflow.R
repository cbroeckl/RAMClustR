test_that("RAMClustR workflow with xcms works", {
  skip_if_not_installed("xcms")
  set.seed(123) # to get reproducible results with jitters
  wd <- getwd()
  tmp <- tempdir()
  load(file.path("testdata", "test.rc.ramclustr.fillpeaks"))
  expected <- readRDS(file.path("testdata", "do.findmain.rds"))
  setwd(tmp)

  ramclustObj <- rc.get.xcms.data(xcmsObj = xdata)
  ramclustObj <- rc.expand.sample.names(ramclustObj = ramclustObj)
  ramclustObj <- rc.feature.replace.na(ramclustObj = ramclustObj)
  ramclustObj <- rc.feature.filter.blanks(ramclustObj = ramclustObj, blank.tag = "Blanc")
  ramclustObj <- rc.feature.normalize.qc(ramclustObj = ramclustObj, qc.tag = "QC")
  ramclustObj <- rc.feature.filter.cv(ramclustObj = ramclustObj)
  ramclustObj <- rc.ramclustr(ramclustObj = ramclustObj)
  ramclustObj <- rc.qc(ramclustObj = ramclustObj)
  actual <- do.findmain(ramclustObj = ramclustObj)

  expect_equal(actual, expected)

  setwd(wd)
})
