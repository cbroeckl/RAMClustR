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

test_that("RAMClustR workflow comparison test", {
  skip_if_not_installed("xcms")
  wd <- getwd()
  tmp <- tempdir()
  load("testdata/test.rc.ramclustr.fillpeaks")
  expected_rds <- readRDS("testdata/rc.ramclustr.rds")
  setwd(tmp)

  ramclustr_obj <- ramclustR(xcmsObj = xdata, maxt = 20, sr = 0.5)

  ramclustObj <- rc.get.xcms.data(xcmsObj = xdata)
  ramclustObj <- rc.feature.replace.na(ramclustObj = ramclustObj)
  ramclustObj <- rc.feature.normalize.tic(ramclustObj = ramclustObj)
  ramclustObj <- rc.ramclustr(ramclustObj = ramclustObj, sr = 0.5, maxt = 20)
  
  expect_equal_labels(ramclustObj$labels, ramclustr_obj$labels)
  expect_equal(ramclustObj$height, ramclustr_obj$height, tolerance = 0.01)
  expect_equal(ramclustObj$frt, ramclustr_obj$frt)
  expect_equal(ramclustObj$fmz, ramclustr_obj$fmz, tolerance = 0.01)
  expect_equal(ramclustObj$SpecAbund, ramclustr_obj$SpecAbund)

  filepaths <- MSnbase::fileNames(xdata)
  filenames <- basename(filepaths)
  dimnames(ramclustr_obj$MSdata)[[1]]<-filenames

  expect_equal_MSdata(ramclustObj$MSdata, ramclustr_obj$MSdata, tolerance = 0.01)

  setwd(wd)
})

