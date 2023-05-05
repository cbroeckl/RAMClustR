#' normalized_data_tic
#'
#' normalize data using TIC
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @return  ramclustR object with total extracted ion normalized data.

normalized_data_tic <- function(ramclustObj = NULL){
  msint <- rowSums(ramclustObj$MSdata, na.rm = TRUE)
  msint.mean <- mean(msint)
  ramclustObj$MSdata <- (ramclustObj$MSdata / msint) * msint.mean

  if (!is.null(ramclustObj$MSMSdata)) {
    msmsint <- rowSums(ramclustObj$MSMSdata, na.rm = TRUE)
    msmsint.mean <- mean(msmsint)
    ramclustObj$MSMSdata <- (ramclustObj$MSMSdata / msmsint) * msmsint.mean
  }

  ramclustObj$history$normalize.tic <- paste0(
    "Features were ",
    if (!is.null(ramclustObj$history$normalize.batch.qc)) {
      "additionally "
    },
    "normalized to total extracted ion signal to account for differences in total solute concentration."
  )

  return(ramclustObj)
}

#' rc.feature.normalize.tic
#'
#' extractor for xcms objects in preparation for clustering
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
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

rc.feature.normalize.tic <- function(ramclustObj = NULL) {
  ## CHECKS
  if (is.null(ramclustObj)) {
    stop(
      "existing ramclustObj required as input", "\n",
      "       see rc.get.xcms.data function for one approach to do so", "\n"
    )
  }

  params <- c()

  ramclustObj <- normalized_data_tic(ramclustObj = ramclustObj)

  ## update msint and optionally msmsint
  global.min <- apply(cbind(ramclustObj$MSdata, ramclustObj$MSMSdata), 2, "min", na.rm = TRUE)

  ramclustObj$msint <- compute_wt_mean(
    ramclustObj$MSdata,
    global.min,
    ramclustObj$fmz,
    TRUE
  )

  if (!is.null(ramclustObj$MSMSdata)) {
    ramclustObj$msmsint <- compute_wt_mean(
      ramclustObj$MSMSdata,
      global.min,
      ramclustObj$fmz,
      TRUE
    )
  }

  if (is.null(ramclustObj$params)) {
    ramclustObj$params <- list()
  }

  ramclustObj$params$rc.feature.normalize.tic <- params

  return(ramclustObj)
}
