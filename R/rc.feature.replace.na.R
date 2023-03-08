#' check_arguments_replace.na
#'
#' check provided arguments
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param replace.int default = 0.1.  proportion of minimum feature value to replace NA (or zero) values with
#' @param replace.noise default = 0.1.  proportion ofreplace.int value by which noise is added via 'jitter'
#' @param replace.zero logical if TRUE, any zero values are replaced with noise as if they were NA values

check_arguments_replace.na <- function(ramclustObj,
                                       replace.int,
                                       replace.noise,
                                       replace.zero) {
  if (is.null(ramclustObj)) {
    stop("please provide a ramclustR Object as input.", "\n")
  }

  if (!is.numeric(replace.int)) {
    stop("replace.int must be numeric", "\n")
  }

  if (!is.numeric(replace.noise)) {
    stop("replace.noise must be numeric", "\n")
  }

  if (!is.logical(replace.zero)) {
    stop("replace.zero must be logical", "\n")
  }
}

#' replace_na
#'
#' add rc.feature.replace.na params in ramclustObj
#'
#' @param data selected data frame to use
#' @param replace.int default = 0.1.  proportion of minimum feature value to replace NA (or zero) values with
#' @param replace.noise default = 0.1.  proportion ofreplace.int value by which noise is added via 'jitter'
#' @param replace.zero logical if TRUE, any zero values are replaced with noise as if they were NA values
#' @return  selected ramclustR data frame with NA and zero values removed.
#' @return  number of features replaced
#' @export

replace_na <- function(data,
                       replace.int,
                       replace.zero,
                       replace.noise) {
  # define a global minimum for the data set to use when all feature values are missing/zero
  min.int.global <- min(data, na.rm = TRUE)

  n.feat.replaced <- 0

  # which values need replacing
  for (i in 1:ncol(data)) {
    rpl <- unique(c(which(is.na(data[, i])), which(is.nan(data[, i])), which(is.infinite(data[, i]))))
    if (replace.zero) {
      rpl.z <- which(data[, i] == 0)
      if (length(rpl.z) > 0) {
        rpl <- unique(c(rpl, rpl.z))
      }
    }
    if (length(rpl) > 0) {
      if (all(is.na(data[, i]))) {
        min.int.local <- min.int.global
      } else {
        min.int.local <- min(data[, i], na.rm = TRUE)
      }
      min.int <- min(min.int.local, min.int.global, na.rm = TRUE)
      rpl.with <- rep((min.int * replace.int), length(rpl))
      rpl.with <- abs(jitter(rpl.with, amount = rpl.with[1] * replace.noise))
      data[rpl, i] <- rpl.with
      n.feat.replaced <- n.feat.replaced + length(rpl)
    }
  }

  return(list(data = data, n.feat.replaced = n.feat.replaced))
}

#' add_params
#'
#' add rc.feature.replace.na params in ramclustObj
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param params vector containing parameters to add
#' @param param_name name of the parameter/step
#' @return  ramclustR object with rc.feature.replace.na params added.
#' @export

add_params <- function(ramclustObj,
                       params,
                       param_name) {
  if (is.null(ramclustObj$params)) {
    ramclustObj$params <- list()
  }
  ramclustObj$params[[param_name]] <- params

  return(ramclustObj)
}

#' rc.feature.replace.na
#'
#' replaces any NA (and optionally zero) values with small signal (20% of minimum feature signal value + 20% random noise)
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param replace.int default = 0.1.  proportion of minimum feature value to replace NA (or zero) values with
#' @param replace.noise default = 0.1.  proportion ofreplace.int value by which noise is added via 'jitter'
#' @param replace.zero logical if TRUE, any zero values are replaced with noise as if they were NA values
#' @param which.data name of dataset
#' @details noise is added by finding for each feature the minimum detected value, multiplying that value by replace.int, then adding (replace.int*replace.noise) noise.  abs() is used to ensure no negative values result.
#' @return  ramclustR object with NA and zero values removed.
#'
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export

rc.feature.replace.na <- function(ramclustObj = NULL,
                                  replace.int = 0.1,
                                  replace.noise = 0.1,
                                  replace.zero = TRUE,
                                  which.data = c("MSdata", "MSMSdata")) {
  check_arguments_replace.na(
    ramclustObj,
    replace.int,
    replace.noise,
    replace.zero
  )

  ########
  # then optionally ensure we have all non-zero values in the dataset.
  # uses a noise addition 'jitter' around minimum values with missing data points.
  # this is mostly necessary for csv input, where other programs may not have used a 'fillPeaks' like step
  # it is important for clustering that variation is present for every feature and MS level.
  n.feat.total <- 0

  for (x in which.data) {
    if (is.null(ramclustObj[[x]])) {
      next
    }

    # select data frame to use
    data <- ramclustObj[[x]]
    n.feat.total <- n.feat.total + (dim(data)[[1]] * dim(data)[[2]])

    replaced_values <- replace_na(
      data,
      replace.int,
      replace.zero,
      replace.noise
    )

    ramclustObj[[x]] <- replaced_values$data
    n.feat.replaced <- replaced_values$n.feat.replaced
  }
  result <- paste(
    "replaced", n.feat.replaced, "of", n.feat.total, "total feature values (",
    round((100 * n.feat.replaced / n.feat.total)), "% )", "\n"
  )
  cat(result)

  ramclustObj$history$replace.na <- {
    paste0(
      "Features with missing values were replaced with small values simulating noise. ",
      "For each feature, the minimum detected value was multiplied by ", replace.int, ". ",
      "Noise was then added using a factor of ", replace.noise, ". ",
      "The absulute value of this value was used as the filled value to ensure that only non-negative values carried forward. ",
      if (replace.zero) {
        "Zero values were treated as missing values."
      }
    )
  }

  ramclustObj <- add_params(
    ramclustObj,
    c(
      "replace.int" = replace.int,
      "replace.noise" = replace.noise,
      "replace.zero" = replace.zero
    ),
    "rc.feature.replace.na"
  )

  return(ramclustObj)
}
