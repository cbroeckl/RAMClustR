expect_equal_labels <- function(actual, expected) {
  actual_features <- split_features(actual)
  expected_features <- split_features(expected)

  expect_equal(actual_features$mz, expected_features$mz)
  expect_equal(actual_features$rt, expected_features$rt)
}

expect_equal_MSdata <- function(actual, expected, tolerance = .Machine$double.eps) {
  actual_labels <- colnames(actual)
  expected_labels <- colnames(expected)
  colnames(actual) <- colnames(expected) <- NULL

  expect_equal_labels(actual_labels, expected_labels)
  expect_equal(actual, expected, tolerance = tolerance)
}

split_features <- function(labels) {
  labels <- strsplit(labels, "_")
  mz_rt_df <- as.data.frame(do.call(rbind, labels))
  mz_rt_df[, 1] <- as.numeric(mz_rt_df[, 1])
  mz_rt_df[, 2] <- as.numeric(mz_rt_df[, 2])
  return(mz_rt_df)
}