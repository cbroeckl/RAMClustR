test_that("RAMClustR workflow with xcms works", {
  skip_if_not_installed("xcms")
  wd <- getwd()
  tmp <- tempdir()
  load("testdata/test.rc.ramclustr.fillpeaks")
  expected_rds <- readRDS("testdata/rc.ramclustr.rds")
  setwd(tmp)

  ramclustObj <- rc.get.xcms.data(xcmsObj = xdata)
  ramclustObj <- rc.expand.sample.names(ramclustObj = ramclustObj)
  ramclustObj <- rc.feature.replace.na(ramclustObj = ramclustObj)
  ramclustObj <- rc.feature.filter.blanks(ramclustObj = ramclustObj, blank.tag = "Blanc")
  ramclustObj <- rc.feature.normalize.qc(ramclustObj = ramclustObj, qc.tag = "QC")
  ramclustObj <- rc.feature.filter.cv(ramclustObj = ramclustObj)
  ramclustObj <- rc.ramclustr(ramclustObj = ramclustObj)
  ramclustObj <- rc.qc(ramclustObj = ramclustObj)
  ramclustObj <- do.findmain(ramclustObj = ramclustObj)

  expect_equal(ramclustObj, expected_rds, tolerance = 0.02)

  setwd(wd)
})
