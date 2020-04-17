#' rc.normalize.feature
#'
#' extractor for xcms objects in preparation for clustering  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param mz numerical vector - accurate mass(es) of feature to use for normalization
#' @param rt numerical vector same length as mz - retention times (in seconds, if xcms input) of feature to use for normalization
#' @details This function normalizes all dataset values relative to the summed intensity of all features selected using mz and rt vectors. 
#' @return  ramclustR object with MSdata and MSMSdata normalized to selected features.   
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

rc.normalize.feature  <- function(
  ramclustObj=NULL,
  mz = NULL,
  mztol = 0.01, 
  rt = NULL, 
  rttol = 3
) {
  
  ## CHECKS
  if(is.null(ramclustObj)) {
    stop('existing ramclustObj required as input', '\n', 
         '       see rc.get.xcms.data function for one approach to do so', '\n')
  }
  
  if(!is.numeric(mz)) {
    stop('mz must be numeric', '\n')
  }
  
  if(!is.numeric(rt)) {
    stop('rt must be numeric', '\n')
  }
  
  if(length(rt) != length(mz)) {
    stop('rt must be same length as rt', '\n')
  }
  
  norm <- rep(0, nrow(ramclustObj$MSdata))
  not.found <- vector(length = 0)
  found <- vector(length = 0)
  for(i in 1:length(mz)) {
    keep <- which(
      abs(ramclustObj$fmz - mz[i]) < mztol &
        abs(ramclustObj$frt - rt[i]) < rttol
    )
    if(length(keep) == 0) {
      not.found <- c(not.found,i)
    } else {
      found <- c(found, i)
    }
    
    if(length(keep) > 1) {
      med <- median(ramclustObj$MSdata[,keep])
      keep <- keep[which.max(med)]
    }
    
    norm <- norm + (ramclustObj$MSdata[,keep]/median(ramclustObj$MSdata[,keep]))
  }
  
  if(length(found) > 0) {
    ramclustObj$MSdata <- ramclustObj$MSdata/norm
    if(!is.null(ramclustObj$MSMSdata)) {
      ramclustObj$MSdata <- ramclustObj$MSdata/norm
    }
    if(!is.null(ramclustObj$MSMSdata)) {
      ramclustObj$MSMSdata <- ramclustObj$MSMSdata/norm
      if(!is.null(ramclustObj$MSMSdata)) {
        ramclustObj$MSMSdata <- ramclustObj$MSMSdata/norm
      }
    }
  } else {
    stop('no features found. please consider adjusting tolerances.', '\n')
  }
  
  
  if(length(found) > 1)
    ramclustObj$history$normalize.feature <- paste0(
      "Feature signal intensities were ",
      "normalized to features mz values = ", 
      paste(mz[found], collapse = ", "), 
      " at retention time(s) = ", 
      paste(rt[found], collapse = ", "), 
      " ,resepectively, to account for differences in total solute concentration."
    )
  
  
  if(length(not.found) > 0) {
    warning("Features not found: ", '\n',
            "  mz values = ", 
            paste(mz[not.found], collapse = ", "), '\n',
            "  rt values = ", 
            paste(rt[found], collapse = ", "), '\n')
  }
  
  ## update msint and optionally msmsint
  msint<-rep(0, length(ramclustObj$fmz))
  for(i in 1:ncol(ramclustObj$MSdata)){
    msint[i]<-weighted.mean(ramclustObj$MSdata[,i], ramclustObj$MSdata[,i], na.rm = TRUE)
  }
  ramclustObj$msint <- msint
  
  if(!is.null(ramclustObj$MSMSdata)) {
    msmsint<-rep(0, length(ramclustObj$fmz))
    for(i in 1:ncol(ramclustObj$MSMSdata)){
      msmsint[i]<-weighted.mean(ramclustObj$MSMSdata[,i], ramclustObj$MSMSdata[,i], na.rm = TRUE)
    }
    ramclustObj$msmsint <- msmsint
  }
  
  
  return(ramclustObj)
}

