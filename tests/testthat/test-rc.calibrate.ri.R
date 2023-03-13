test_that("RAMClustR rc.calibrate.ri", {
  calibrant.data <- file.path("testdata", "calibrant.data.csv")
  ramclustObj <- readRDS(file.path("testdata", "input.rc.calibrate.ri.rds"))
  expected <- readRDS(file.path("testdata", "rc.calibrate.ri.rds"))

  actual <- rc.calibrate.ri(
    ramclustObj = ramclustObj,
    calibrant.data = calibrant.data
  )
  expect_equal(actual, expected)
})
