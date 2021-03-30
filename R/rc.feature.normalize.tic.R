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

rc.feature.normalize.tic  <- function(
  ramclustObj=NULL
) { 
  
  ## CHECKS
  if(is.null(ramclustObj)) {
    stop('existing ramclustObj required as input', '\n', 
         '       see rc.get.xcms.data function for one approach to do so', '\n')
  }
  
  params <- c(
  )

  msint <- rowSums(ramclustObj$MSdata, na.rm=TRUE)
  msint.mean <- mean(msint)
  ramclustObj$MSdata <- (ramclustObj$MSdata/msint)*msint.mean

  if(!is.null(ramclustObj$MSMSdata)) {
    msmsint <- rowSums(ramclustObj$MSMSdata, na.rm=TRUE)
    msmsint.mean <- mean(msmsint)
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
  
  if(is.null(ramclustObj$params)) {ramclustObj$params <- list()}
  ramclustObj$params$rc.feature.normalize.tic <- params
  
  return(ramclustObj)
}

