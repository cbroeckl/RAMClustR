test_that("RAMClustR rc.feature.normalize.qc", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.filter.blanks.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.normalize.qc.rds"))
  metadata <- read_metadata(file.path("testdata", "testMetadata.csv"))

  actual <- rc.feature.normalize.qc(
    ramclustObj = ramclustObj,
    batch = metadata$batch,
    order = metadata$order,
    qc = metadata$qc
  )

  actual$history <- NA
  expected$history <- NA
  actual$params$rc.feature.normalize.qc <- NA
  expected$params$rc.feature.normalize.qc <- NA

  expect_equal(actual, expected)
})
