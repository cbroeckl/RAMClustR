#' rc.feature.normalize.tic
#'
#' extractor for xcms objects in preparation for clustering  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.
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

rc.feature.normalize.tic  <- function(
  ramclustObj=NULL
) { 
  
  ## CHECKS
  if(is.null(ramclustObj)) {
    stop('existing ramclustObj required as input', '\n', 
         '       see rc.get.xcms.data function for one approach to do so', '\n')
  }
  
  ramclustObj$MSdata <- (ramclustObj$MSdata/rowSums(ramclustObj$MSdata))*mean(rowSums(ramclustObj$MSdata), na.rm=TRUE)
  if(!is.null(ramclustObj$MSMSdata)) {
    ramclustObj$MSMSdata<-(ramclustObj$MSMSdata/rowSums(ramclustObj$MSMSdata))*mean(rowSums(ramclustObj$MSMSdata), na.rm=TRUE)
  }
  
  
  ramclustObj$history$normalize.tic <- paste0(
    "Features were ",
    if(!is.null(ramclustObj$history$normalize.batch.qc)) {"additionally "}, 
    "normalized to total extracted ion signal to account for differences in total solute concentration."
  )
  
  
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

