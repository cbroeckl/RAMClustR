#' rc.feature.filter.cv
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

rc.feature.filter.cv  <- function(
  ramclustObj=NULL,
  qc.tag = "QC",
  max.cv = 0.5
) { 
  
  ## CHECKS
  if(is.null(ramclustObj)) {
    stop('existing ramclustObj required as input', '\n', 
         '       see rc.get.xcms.data function for one approach to do so', '\n')
  }
  
  if(is.null(qc.tag)) {
    stop("qc.tag = NULL; qc.tag must be defined to enable QC variance examination.", '\n')
  }
  
  if(!is.null(ramclustObj$SpecAbund)) {
    stop('cannot perform feature removal after clustering', '\n', 
         '       please run rc.feature.cv.filter before clustering', '\n')
  }
  
  params <- c(
    "qc.tag" = qc.tag,
    "max.cv" = max.cv
  )
  
  do.sets <- c("MSdata", "MSMSdata")
  if(is.null(ramclustObj$MSMSdata)) {
    do.sets <- do.sets[!(do.sets %in% "MSMSdata")]
  } 
  
  do.sets.rows <- sapply(
    c(do.sets, "phenoData"), 
    FUN = function(x) {
      nrow(ramclustObj[[x]])
    })
  
  if(!all.equal(
    do.sets.rows, do.sets.rows)
  ) {
    stop("number of rows in MSdata, SpecAbund, and phenoData sets are not identical.")
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
  keep <- rep(FALSE, ncol(ramclustObj$MSdata))
  par(mfrow = c(1, length(do.sets)))
  for(x in do.sets) {
    td <- ramclustObj[[x]]
    sds<-apply(td[qc,], 2, FUN="sd", na.rm=TRUE)
    #cat(sds, '\n')
    means<-apply(td[qc,], 2, FUN="mean", na.rm=TRUE)
    cvs<-sds/means
    hist(cvs, main = x, xlab = "CV")
    keep[which(cvs <= max.cv)] <- TRUE
    hist(cvs[which(cvs <= max.cv)], add = TRUE, col = "gray")
  }
  
  ## filter to keep only 'good' features
  
  ramclustObj$MSdata <- ramclustObj$MSdata[,keep]
  ramclustObj$msint <- ramclustObj$msint[keep]
  if(!is.null(ramclustObj$MSdata_raw)) {
    ramclustObj$MSdata_raw <- ramclustObj$MSdata_raw[,keep]
  }
  
  if(!is.null(ramclustObj$MSMSdata)) {
    ramclustObj$MSMSdata <- ramclustObj$MSMSdata[,keep]
    ramclustObj$msmsint <- ramclustObj$msmsint[keep]
    if(!is.null(ramclustObj$MSMSdata_raw)) {
      ramclustObj$MSMSdata_raw <- ramclustObj$MSMSdata_raw[,keep]
    }
  }
  
  ramclustObj$frt <- ramclustObj$frt[keep]
  ramclustObj$fmz <- ramclustObj$fmz[keep]
  ramclustObj$featnames <- ramclustObj$featnames[keep]
  ramclustObj$xcmsOrd <- ramclustObj$xcmsOrd[keep]
  
  if(!is.null(ramclustObj$qc$MSdata)) {
    ramclustObj$qc$MSdata_raw <- ramclustObj$qc$MSdata_raw[,keep]
  }
  if(!is.null(ramclustObj$qc$MSMSdata)) {
    ramclustObj$qc$MSMSdata_raw <- ramclustObj$qc$MSMSdata_raw[,keep]
  }
  
  ramclustObj$history$filter.features <- paste0(
    "Features were filtered based on their qc sample CV values.",
    " Only features with CV vaules less than or equal to ", max.cv, 
    if(is.null(ramclustObj$MSMSdata)) {" in MSdata set"} else {" in MS or MSMSdata sets"},
    " were retained. ",   length(which(!keep))," of ", length(keep), " features were removed."
  )
  
  if(is.null) {ramclustObj$params <- list()}
  ramclustObj$params$rc.feature.filter.cv <- params
  
  cat(ramclustObj$history$filter.features)
  
  return(ramclustObj)
}

