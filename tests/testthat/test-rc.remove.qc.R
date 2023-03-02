test_that("RAMClustR rc.remove.qc", {
  ramclustObj <- readRDS(file.path("testdata", "rc.get.xcms.data.rds"))
  expected <- readRDS(file.path("testdata", "rc.remove.qc.rds"))
  
  actual <- rc.remove.qc(ramclustObj = ramclustObj)

  expect_equal(actual, expected)
})
