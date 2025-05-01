test_that("RAMClustR rc.feature.normalize.batch.qc", {
  ramclustObj <- readRDS(file.path("testdata", "rc.feature.filter.blanks.rds"))
  expected <- readRDS(file.path("testdata", "rc.feature.normalize.batch.qc.rds"))
  metadata <- read_metadata(file.path("testdata", "testMetadata.csv"))

  actual <- rc.feature.normalize.batch.qc(
    order = metadata$order,
    batch = metadata$batch,
    qc = metadata$qc,
    ramclustObj = ramclustObj,
    qc.inj.range = 20
  )

  expect_equal(actual, expected)
})
