test_that("RAMClustR rc.get.df.data", {
  phenoData <- readRDS(file.path("testdata", "phenoData_df.rds"))
  expected <- readRDS(file.path("testdata", "rc.get.df.data.rds"))
  ms1_featureDefinitions <- readRDS(file.path("testdata", "ms1_featureDefinitions.rds"))
  ms1_featureValues <- readRDS(file.path("testdata", "ms1_featureValues.rds"))

  actual <- rc.get.df.data(
    ms1_featureDefinitions = ms1_featureDefinitions,
    ms1_featureValues = ms1_featureValues,
    phenoData = phenoData,
    st = 5
  )
  expect_equal(actual, expected)
})
