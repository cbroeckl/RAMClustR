#' rc.restore.qc.samples
#'
#' summarize quality control for clustering and for quality control sample variation based on compound ($SpecAbund) and feature ($MSdata and $MSMSdata, if present)
#'
#' @param ramclustObj ramclustR object to analyze
#'
#' @details moves all of $phenoData, $MSdata, $MSMSdata, $SpecAbund back to original positions from $qc slot
#' @return   RC object  
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

rc.restore.qc.samples<-function(
  ramclustObj=NULL
){
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(is.null(ramclustObj$qc)) {
    stop("No qc slot present in ramclustObj.", '\n')
  }
  
  do.sets <- c("MSdata", "MSMSdata", "SpecAbund")
  if(is.null(ramclustObj$MSMSdata)) {
    do.sets <- do.sets[!(do.sets %in% "MSMSdata")]
  } 
  if(is.null(ramclustObj$SpecAbund)) {
    do.sets <- do.sets[!(do.sets %in% "SpecAbund")]
  } 
  
  do.sets.rows <- sapply(
    c(do.sets, "phenoData"), 
    FUN = function(x) {
      nrow(ramclustObj[[x]])
    })
  
  if(!all.equal(
    do.sets.rows, do.sets.rows)
  ) {
    stop("number of rows in sample sets are not identical.")
  }
  
  do.sets <- names(ramclustObj$qc)
  for(x in c(do.sets)) {
    ramclustObj[[x]] <- rbind(ramclustObj[[x]], ramclustObj$qc[[x]])
  }
  
  ramclustObj$qc <- NULL
  
  ord.pheno <- order(as.numeric(row.names(ramclustObj$phenoData)))
  for(x in c(do.sets)) {
    ramclustObj[[x]] <- ramclustObj[[x]][ord.pheno,]
  }
  return(ramclustObj)
}


