test_that("RAMClustR workflow comparison test", {
  skip_if_not_installed("xcms")
  wd <- getwd()
  tmp <- tempdir()
  load(file.path("testdata", "test.fillpeaks"))
  setwd(tmp)
  
  ramclustr_obj <- ramclustR(xcmsObj = xdata, maxt = 20, sr = 0.5)

  ramclustObj <- rc.get.xcms.data(xcmsObj = xdata)
  ramclustObj <- rc.feature.replace.na(ramclustObj = ramclustObj)
  ramclustObj <- rc.feature.normalize.tic(ramclustObj = ramclustObj)
  ramclustObj <- rc.ramclustr(ramclustObj = ramclustObj, sr = 0.5, maxt = 20)

  expect_equal_labels(ramclustObj$labels, ramclustr_obj$labels)
  expect_equal(ramclustObj$height, ramclustr_obj$height)
  expect_equal(ramclustObj$frt, ramclustr_obj$frt)
  expect_equal(ramclustObj$fmz, ramclustr_obj$fmz)
  expect_equal(ramclustObj$SpecAbund, ramclustr_obj$SpecAbund, tolerance = 0.01)

  filepaths <- MSnbase::fileNames(xdata)
  filenames <- basename(filepaths)
  dimnames(ramclustr_obj$MSdata)[[1]] <- filenames

  expect_equal(ramclustObj$MSdata, ramclustr_obj$MSdata, tolerance = 0.01)

  setwd(wd)
})
