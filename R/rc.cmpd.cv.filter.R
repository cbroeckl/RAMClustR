#' rc.cmpd.filter.cv
#'
#' extractor for xcms objects in preparation for clustering  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.
#' @param max.cv numeric maximum allowable cv for any feature.  default = 0.3 
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

rc.cmpd.filter.cv  <- function(
  ramclustObj=NULL,
  qc.tag = "QC",
  max.cv = 0.5
) { 
  
  params <- c(
    qc.tag = "QC",
    max.cv = max.cv
  )
  
  ## CHECKS
  if(is.null(ramclustObj)) {
    stop('existing ramclustObj required as input', '\n', 
         '       see rc.get.xcms.data function for one approach to do so', '\n')
  }
  
  if(is.null(qc.tag)) {
    stop("qc.tag = NULL; qc.tag must be defined to enable QC variance examination.", '\n')
  }
  
  if(is.null(ramclustObj$SpecAbund)) {
    stop('cannot perform compound removal before clustering', '\n', 
         '       please run rc.ramclustr before clustering', '\n')
  }
  

  if(is.null(ramclustObj$SpecAbundAve)) {
    do.sets <- "SpecAbund"
  } else {
    do.sets <- c("SpecAbund", "SpecAbundAve")
  }

  ## define QC samples in each set
  if(length(qc.tag) == 1) {
    qc <- grepl(qc.tag[1], ramclustObj$phenoData$sample.names)
  } 
  if(length(qc.tag) == 2) {
    qc <- grepl(qc.tag[1], ramclustObj$phenoData[[qc.tag[2]]])
  }
  
  if(length(which(qc)) == 0) {
    stop("no QC samples found using the qc.tag ", "'", qc.tag, "'", '\n')
  }
  
  ## find 'good' features, acceptable CV at either
  ## MS or MSMS level results in keeping
  keep <- rep(FALSE, ncol(ramclustObj$SpecAbund))
  # par(mfrow = c(1, length(do.sets)))
  for(x in do.sets) {
    td <- ramclustObj[[x]]
    sds<-apply(td[qc,], 2, FUN="sd", na.rm=TRUE)
    #cat(sds, '\n')
    means<-apply(td[qc,], 2, FUN="mean", na.rm=TRUE)
    cvs<-sds/means
    # hist(cvs, main = x, xlab = "CV")
    keep[which(cvs <= max.cv)] <- TRUE
    # hist(cvs[which(cvs <= max.cv)], add = TRUE, col = "gray")
  }
  
  ## subset all vectors, matrices, data.frames and lists to keep only those of interest
  for(x in names(ramclustObj)) {
    if(is.vector(ramclustObj[[x]]) | (is.list(ramclustObj[[x]])) & !is.data.frame(ramclustObj[[x]])) {
      if(length(ramclustObj[[x]]) == length(cvs)) {
        ramclustObj[[x]] <- ramclustObj[[x]][keep]
      }
    }
    
    if(is.matrix(ramclustObj[[x]]) | is.data.frame(ramclustObj[[x]])) {
      if(dim(ramclustObj[[x]])[[2]] == length(cvs)) {
        ramclustObj[[x]] <- ramclustObj[[x]][,keep]
      }
    }
  }
  
  if(is.null(ramclustObj$params)) {ramclustObj$params <- list()}
  ramclustObj$params$rc.cmpd.cv.filter <- params
  
  ramclustObj$history$filter.cmpds <- paste0(
    "Compounds were filtered based on their qc sample CV values.",
    " All compounds with CV vaules greater than ", max.cv, 
    " in the cluster intensity (SpecAbund) dataset",
    " were removed. ",   length(which(!keep))," of ", length(keep), " compounds were removed."
  )
  
  cat(ramclustObj$history$filter.cmpds)
  
  return(ramclustObj)
}
