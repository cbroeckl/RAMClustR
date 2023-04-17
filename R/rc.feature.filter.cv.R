#' check_arguments_filter.cv
#'
#' check provided arguments
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.

check_arguments_filter.cv <- function(ramclustObj, qc.tag) {
  ## CHECKS
  if (is.null(ramclustObj)) {
    stop(
      "existing ramclustObj required as input", "\n",
      "       see rc.get.xcms.data function for one approach to do so", "\n"
    )
  }

  if (is.null(qc.tag)) {
    stop("qc.tag = NULL; qc.tag must be defined to enable QC variance examination.", "\n")
  }

  if (!is.null(ramclustObj$SpecAbund)) {
    stop(
      "cannot perform feature removal after clustering", "\n",
      "       please run rc.feature.cv.filter before clustering", "\n"
    )
  }
}

#' define_samples
#'
#' define samples in each set
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.
#' @return samples found using the tag
#' @export

define_samples <- function(ramclustObj, tag) {
  ## define samples in each set
  if(length(tag) == 0) {
    stop("no tag provided", "\n")
  }

  samples <- grep(tag[1], ramclustObj$sample_names)
  samples <- samples[which(samples <= nrow(ramclustObj$MSdata))]

  if (length(samples) == 0) {
    stop("no ", tag, " samples found using the tag ", "'", tag, "'", "\n")
  }
  return(samples)
}

#' find_good_features
#'
#' find 'good' features, acceptable CV at either MS or MSMS level results in keeping
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param do.sets select data frame to use.
#' @param max.cv numeric maximum allowable cv for any feature.  default = 0.5
#' @param qc QC samples found by define_samples
#' @return ramclustR object
#' @return features to keep
#' @export

find_good_features <- function(ramclustObj,
                               do.sets,
                               max.cv,
                               qc) {
  ## find 'good' features, acceptable CV at either
  ## MS or MSMS level results in keeping
  keep <- rep(FALSE, ncol(ramclustObj$MSdata))
  par(mfrow = c(1, length(do.sets)))
  for (x in do.sets) {
    td <- ramclustObj[[x]]
    sds <- apply(td[qc, ], 2, FUN = "sd", na.rm = TRUE)
    means <- apply(td[qc, ], 2, FUN = "mean", na.rm = TRUE)
    cvs <- sds / means
    if (x == "MSdata") {
      ramclustObj$qc.cv.feature.msdata <- cvs
    } else {
      ramclustObj$qc.cv.feature.msmsdata <- cvs
      ramclustObj$qc.cv.feature <- pmin(
        ramclustObj$qc.cv.feature.msdata,
        ramclustObj$qc.cv.feature.msmsdata
      )
    }
    hist(cvs, main = x, xlab = "CV")
    keep[which(cvs <= max.cv)] <- TRUE
    cat(x, ": ", length(which(cvs <= max.cv)), "passed the CV filter", "\n")
    hist(cvs[which(cvs <= max.cv)], add = TRUE, col = "gray")
  }

  return(list(ramclustObj = ramclustObj, keep = keep))
}

#' filter_good_features
#'
#' filter to keep only 'good' features
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param keep features to keep. output of find_good_features().
#' @return ramclustR object filtered to keep only 'good' features
#' @export

filter_good_features <- function(ramclustObj, keep) {
  ## filter to keep only 'good' features
  l <- sapply(1:length(ramclustObj), FUN = function(x) length(ramclustObj[[x]]))
  l <- which(l == ncol(ramclustObj$MSdata))
  for (i in l) {
    ramclustObj[[i]] <- ramclustObj[[i]][keep]
  }

  nc <- sapply(1:length(ramclustObj), FUN = function(x) ncol(ramclustObj[[x]]))
  n <- sapply(1:length(nc), FUN = function(x) is.null(nc[[x]]))
  for (i in 1:length(n)) {
    if (n[i]) nc[[i]] <- 0
  }
  nc <- unlist(nc)
  nc <- which(nc == ncol(ramclustObj$MSdata))
  for (i in nc) {
    ramclustObj[[i]] <- ramclustObj[[i]][, keep]
  }

  return(ramclustObj)
}

#' compute_do.sets
#'
#' compute data frame to use in ramclustObj 
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @return vector which is used to select data frame to use in ramclustObj
#' @export

compute_do.sets <- function(ramclustObj) {
  do.sets <- c("MSdata", "MSMSdata")
  if (is.null(ramclustObj$MSMSdata)) {
    do.sets <- do.sets[!(do.sets %in% "MSMSdata")]
  }

  do.sets.rows <- sapply(
    c(do.sets, "phenoData"),
    FUN = function(x) {
      nrow(ramclustObj[[x]])
    }
  )

  if (!all.equal(
    do.sets.rows, do.sets.rows
  )
  ) {
    stop("number of rows in MSdata, SpecAbund, and phenoData sets are not identical.")
  }

  return(do.sets)
}

#' rc.feature.filter.cv
#'
#' extractor for xcms objects in preparation for clustering
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.
#' @param max.cv numeric maximum allowable cv for any feature.  default = 0.5
#' @details This function offers normalization by total extracted ion signal.  it is recommended to first run 'rc.feature.filter.blanks' to remove non-sample derived signal.
#' @return  ramclustR object with total extracted ion normalized data.
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

rc.feature.filter.cv <- function(ramclustObj = NULL,
                                 qc.tag = "QC",
                                 max.cv = 0.5) {
  check_arguments_filter.cv(ramclustObj, qc.tag)

  do.sets <- compute_do.sets(ramclustObj)

  qc <- define_samples(ramclustObj, qc.tag)

  good_features <- find_good_features(
    ramclustObj,
    do.sets,
    max.cv,
    qc
  )

  ramclustObj <- good_features$ramclustObj
  keep <- good_features$keep

  ramclustObj <- filter_good_features(ramclustObj, keep)

  # ramclustObj$MSdata <- ramclustObj$MSdata[,keep]
  # ramclustObj$msint <- ramclustObj$msint[keep]
  # if(!is.null(ramclustObj$MSdata_raw)) {
  #   ramclustObj$MSdata_raw <- ramclustObj$MSdata_raw[,keep]
  # }
  #
  # if(!is.null(ramclustObj$MSMSdata)) {
  #   ramclustObj$MSMSdata <- ramclustObj$MSMSdata[,keep]
  #   ramclustObj$msmsint <- ramclustObj$msmsint[keep]
  #   if(!is.null(ramclustObj$MSMSdata_raw)) {
  #     ramclustObj$MSMSdata_raw <- ramclustObj$MSMSdata_raw[,keep]
  #   }
  # }
  #
  # ramclustObj$frt <- ramclustObj$frt[keep]
  # ramclustObj$fmz <- ramclustObj$fmz[keep]
  # ramclustObj$featnames <- ramclustObj$featnames[keep]
  # ramclustObj$xcmsOrd <- ramclustObj$xcmsOrd[keep]
  #
  # if(!is.null(ramclustObj$qc$MSdata)) {
  #   ramclustObj$qc$MSdata_raw <- ramclustObj$qc$MSdata_raw[,keep]
  # }
  # if(!is.null(ramclustObj$qc$MSMSdata)) {
  #   ramclustObj$qc$MSMSdata_raw <- ramclustObj$qc$MSMSdata_raw[,keep]
  # }

  ramclustObj$history$filter.features <- paste0(
    "Features were filtered based on their qc sample CV values.",
    " Only features with CV vaules less than or equal to ", max.cv,
    if (is.null(ramclustObj$MSMSdata)) {
      " in MSdata set"
    } else {
      " in MS or MSMSdata sets"
    },
    " were retained. ", length(which(!keep)), " of ", length(keep), " features were removed."
  )

  ramclustObj <- add_params(
    ramclustObj,
    c(
      "qc.tag" = qc.tag,
      "max.cv" = max.cv
    ),
    "rc.feature.filter.cv"
  )

  cat(ramclustObj$history$filter.features)

  return(ramclustObj)
}
