test_that("RAMClustR with xcms works", {
  skip_if_not_installed("xcms")
  set.seed(123) # to get reproducible results with jitters
  wd <- testthat::test_path()
  tmp <- tempdir()
  load(file.path("testdata", "test.fillpeaks"))
  expected <- readLines(file.path("testdata", "output.msp"))
  
  ramclustr_obj <- ramclustR(xcmsObj = xdata, maxt = 20, sr = 0.5, mzdec = 4, out.dir = tmp)
  write.msp(ramclustr_obj, one.file = TRUE, out.dir = tmp)
  mismatches <- setdiff(expected, readLines(file.path(tmp, "spectra", "fill.msp")))
  
  expect_true(length(mismatches) < 10)
})

test_that("RAMClustR with csv works", {
  tmp <- tempdir()
  filename <- file.path("testdata", "peaks.csv")
  pheno <- file.path("testdata", "phenoData.csv")
  # load(file.path("testdata", "test_csv.Rds"))

  actual <- ramclustR(
    ms = filename,
    pheno_csv = pheno,
    st = 5,
    maxt = 1,
    blocksize = 1000,
    out.dir = tmp
  )
  
  expected <- readRDS(file = file.path("testdata", "test_csv.rds"))
  
  actual$history <- NA
  expected$history <- NA
  
  expect_equal_labels(actual$labels, expected$labels)
  expect_equal_MSdata(actual$MSdata, expected$MSdata)

  actual$labels <- expected$labels <- actual$MSdata <- expected$MSdata <- NA
  
  expect_equal(actual, expected)
})

# devtools::test_active_file()
