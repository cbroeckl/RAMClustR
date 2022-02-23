test_that("RAMClustR works", {
  load("testdata/xcmsObj.rdata.xcms.fillpeaks")
  
  ramclustr_obj <- ramclustR(xcmsObj = xdata, blocksize = 200)
  write.msp(ramclustr_obj, one.file = TRUE)
  actual_msp = "spectra/fill.msp"
  expected_msp = "testdata/output.msp"
  
  expect_equal(readLines(expected_msp), readLines(actual_msp))
})