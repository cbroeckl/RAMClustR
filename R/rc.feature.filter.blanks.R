#' check_arguments_filter.blanks
#'
#' check provided arguments
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param sn numeric defines the ratio for 'signal'.  i.e. sn = 3 indicates that signal intensity must be 3 fold higher in sample than in blanks, on average, to be retained.

check_arguments_filter.blanks <- function(ramclustObj, sn) {
  if (is.null(ramclustObj)) {
    stop("please provide a ramclustR Object as input.", "\n")
  }

  if (length(ramclustObj$dist.method) > 0) {
    stop(
      "this data has been clustered.", "\n",
      "       please perform feature filtering before clustering.", "\n"
    )
  }

  if (!is.numeric(sn)) {
    stop("sn must be numeric", "\n")
  }
}

#' mean_signal_intensities
#'
#' calculate MS mean signal intensities
#'
#' @param data MS/MSMS data
#' @param sample sample found using the tag, output of define_samples()
#' @return mean signal intensities
#' @export

mean_signal_intensities <- function(data, sample) {
  if (length(sample) > 1) {
    ms_sample_mean <- apply(data[sample, ], 2, FUN = "mean", na.rm = TRUE)
  } else {
    ms_sample_mean <- data[sample, ]
  }

  return(ms_sample_mean)
}

#' filter_signal
#'
#' filter signal
#'
#' @param ms.qc.mean ms qc mean signal intensities
#' @param ms.blank.mean ms blank mean signal intensities
#' @param sn numeric defines the ratio for 'signal'.  i.e. sn = 3 indicates that signal intensity must be 3 fold higher in sample than in blanks, on average, to be retained.
#' @return union of which signal is at least 3x larger
#' @export

filter_signal <- function(ms.qc.mean, ms.blank.mean, sn) {
  # filters
  # which signal is at least 3x larger in QC
  keep.ms.a <- which(ms.qc.mean / ms.blank.mean >= sn)
  # which signal is present in QC and absent in blank
  keep.ms.b <- which((!is.nan(ms.qc.mean)) & is.nan(ms.blank.mean))
  # union of these two are keepers
  keep <- union(keep.ms.a, keep.ms.b)

  return(keep)
}

#' filter_blanks
#'
#' filter blanks
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param keep union of which signal is at least 3x larger, output of filter_signal()
#' @param d1 MS Data
#' @return ramclustObj object with feature.filter.blanks
#' @export

filter_blanks <- function(ramclustObj, keep, d1) {
  ## define variables that we need to subset
  do <- rep(FALSE, length(ramclustObj))
  rf <- data.frame(index = 1:ncol(d1))
  rf <- rf[, 0]
  for (i in 1:length(do)) {
    if (any(names(ramclustObj)[i] %in% c("MSdata_raw", "MSMSdata_raw"))) {
      next
    }
    if (is.vector(ramclustObj[[i]])) {
      if (length(ramclustObj[[i]]) == ncol(d1)) {
        rf <- data.frame(rf, ramclustObj[[i]])
        names(rf)[ncol(rf)] <- names(ramclustObj[i])
        ramclustObj[[i]] <- ramclustObj[[i]][keep]
      }
    }

    if (is.list(ramclustObj[[i]])) {
      if (length(ramclustObj[[i]]) == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][keep]
      }
    }

    if (is.data.frame(ramclustObj[[i]])) {
      cat("df", names(ramclustObj)[i], "\n")
      if (dim(ramclustObj[[i]])[2] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][, keep]
      }
      if (dim(ramclustObj[[i]])[1] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][keep, ]
      }
    }

    if (is.matrix(ramclustObj[[i]])) {
      cat("ma", names(ramclustObj)[i], "\n")
      if (dim(ramclustObj[[i]])[2] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][, keep]
      }
      if (dim(ramclustObj[[i]])[1] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][keep, ]
      }
    }
  }
  ramclustObj$feature.filter.blanks <- rf[-keep, ]

  return(ramclustObj)
}

#' remove_blanks
#'
#' remove blanks
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param blank blank samples found by define_samples
#' @return ramclustObj object with blanks removed
#' @export

remove_blanks <- function(ramclustObj, blank) {
  ramclustObj$MSdata <- ramclustObj$MSdata[-blank, ]
  if (!is.null(ramclustObj$MSMSdata)) {
    ramclustObj$MSMSdata <- ramclustObj$MSMSdata[-blank, ]
  }
  if (!is.null(ramclustObj$phenoData)) {
    ramclustObj$phenoData <- ramclustObj$phenoData[-blank, ]
  }

  ramclustObj$sample_names <- ramclustObj$sample_names[-blank]

  return(ramclustObj)
}

#' rc.feature.filter.blanks
#'
#' used to remove features which are found at similar intensity in blank samples
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.
#' @param blank.tag see 'qc.tag' , but for blanks to use as background.
#' @param sn numeric defines the ratio for 'signal'.  i.e. sn = 3 indicates that signal intensity must be 3 fold higher in sample than in blanks, on average, to be retained.
#' @param remove.blanks logical. TRUE by default.  this removes any recognized blanks samples from the MSdata and MSMSdata sets after they are used to filter contaminant features.
#' @details This function offers normalization by run order, batch number, and QC sample signal intensity.
#' @details Each input vector should be the same length, and equal to the number of samples in the $MSdata set.
#' @details Input vector order is assumed to be the same as the sample order in the $MSdata set.
#' @return  ramclustR object with normalized data.
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

rc.feature.filter.blanks <- function(ramclustObj = NULL,
                                     qc.tag = "QC",
                                     blank.tag = "blank",
                                     sn = 3,
                                     remove.blanks = TRUE) {
  check_arguments_filter.blanks(ramclustObj, sn)

  if (is.null(ramclustObj$MSMSdata)) {
    msms <- FALSE
  } else {
    msms <- TRUE
    d2 <- ramclustObj$MSMSdata
  }

  ## MS data as reference.  assumes that sample order is the same for MSMS data.
  d1 <- ramclustObj$MSdata

  ## define QC and Blank samples in each set
  qc <- define_samples(ramclustObj, qc.tag)
  blank <- define_samples(ramclustObj, blank.tag)

  ## create logical vector of features to keep
  ## 'good' features should have signal intensity
  ## in pooled QC and/or Samples that is at least
  ## sn higher than the process blank

  #  MS1 mean signal intensities

  ms1.qc.mean <- mean_signal_intensities(d1, qc)
  ms1.blank.mean <- mean_signal_intensities(d1, blank)

  absent.in.blank <- which(is.nan(ms1.blank.mean))

  # filters
  # which signal is at least 3x larger in QC
  keep <- filter_signal(ms1.qc.mean, ms1.blank.mean, sn)

  ## do the same for MS2
  if (msms) {
    ms2.qc.mean <- mean_signal_intensities(d2, qc)
    ms2.blank.mean <- mean_signal_intensities(d2, blank)

    absent.in.blank <- which(is.nan(ms2.blank.mean))

    keep <- filter_signal(ms2.qc.mean, ms2.blank.mean, sn)
  }

  ## union of keep.ms1 and keep.ms2 is what we want to move forward
  length(keep) / ncol(d1)
  cat(round(100 * length(keep) / ncol(d1), digits = 1), "% of features move forward", "\n", sep = "")

  ramclustObj <- filter_blanks(ramclustObj, keep, d1)

  if (remove.blanks) {
    ramclustObj <- remove_blanks(ramclustObj, blank)
  }

  ramclustObj$history$feature.filter.blanks <- {
    paste0(
      "Features which failed to demonstrate signal intensity of at least ",
      sn, " fold greater in QC samples than in blanks were removed from the feature dataset. ",
      ncol(d1) - length(keep), " of ", ncol(d1), " features were removed."
    )
  }

  cat(ramclustObj$history$feature.filter.blanks)

  return(ramclustObj)
}
