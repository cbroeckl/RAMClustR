patrick::with_parameters_test_that("RAMClustR rc.ramclustr", {
    ramclustObj <- readRDS(file.path("testdata", "rc.feature.filter.cv.rds"))
    expected_path <- file.path("testdata", "clustering", paste0(.test_name,"_rc.ramclustr.rds"))
    expected <- readRDS(expected_path)

    actual <- rc.ramclustr(ramclustObj = ramclustObj, rt.only.low.n=rt_only_low_n)
    
    expect_equal(actual, expected)
  },
  rt_only_low_n = list(TRUE, FALSE),
  .test_name = list("only_rt", "normal")
)
