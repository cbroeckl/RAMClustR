#' rc.remove.qc
#'
#' summarize quality control for clustering and for quality control sample variation based on compound ($SpecAbund) and feature ($MSdata and $MSMSdata, if present)
#'
#' @param ramclustObj ramclustR object to analyze
#' @param qc.tag qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.
#'
#' @details simply moves QC samples out of the way for downstream processing. moved to a $qc slot.
#' @return new RC object. moves QC samples to new $qc slot from original position.
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
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

rc.remove.qc <- function(ramclustObj = NULL,
                         qc.tag = "QC") {
  params <- c(
    "qc.tag" = qc.tag
  )


  if (is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", "\n")
  }

  if (is.null(qc.tag)) {
    stop("qc.tag = NULL; qc.tag must be defined to enable QC variance examination.", "\n")
  }

  do.sets <- c("MSdata", "MSMSdata", "SpecAbund")
  if (is.null(ramclustObj$MSMSdata)) {
    do.sets <- do.sets[!(do.sets %in% "MSMSdata")]
  }
  if (is.null(ramclustObj$SpecAbund)) {
    do.sets <- do.sets[!(do.sets %in% "SpecAbund")]
  }

  do.sets.rows <- sapply(
    c(do.sets, "phenoData", "sample_names"),
    FUN = function(x) {
      NROW(ramclustObj[[x]])
    }
  )

  if (!all.equal(
    do.sets.rows, do.sets.rows
  )
  ) {
    stop("number of rows in MSdata, SpecAbund, and phenoData sets are not identical.")
  }
  
  ## define QC samples in each set
  qc <- define_samples(ramclustObj, qc.tag)

  ramclustObj$qc <- list()
  for (x in c("phenoData", "sample_names",do.sets)) {
    if (is.vector(ramclustObj[[x]])) {
      ramclustObj$qc[[x]] <- ramclustObj[[x]][qc]
      ramclustObj[[x]] <- ramclustObj[[x]][-qc]
    } else {
      ramclustObj$qc[[x]] <- ramclustObj[[x]][qc, ]
      ramclustObj[[x]] <- ramclustObj[[x]][-qc, ]
    }
  }

  ramclustObj$history$qc.summary <- paste(
    ramclustObj$history$qc.summary,
    "QC samples were removed from the set for downstream processing."
  )

  if (is.null((ramclustObj$params))) {
    ramclustObj$params <- list()
  }
  ramclustObj$params$rc.remove.qc <- params

  return(ramclustObj)
}
