#' exportDataset
#'
#' export one of 'SpecAbund', 'SpecAbundAve', 'MSdata' or 'MSMSdata' from an RC object to csv
#' @details Useful for exporting the processed signal intensity matrix to csv for analysis elsewhere.  
#' @param ramclustObj ramclustR object to export from
#' @param which.data name of dataset to export.  SpecAbund, SpecAbundAve, MSdata, or MSMSdata
#' @param label.by either 'ann' or 'cmpd', generally.  name of ramclustObj slot used as csv header for each column (compound)
#' @param appendFactors logical.  If TRUE (default) the factor data frame is appended to the left side of the dataset. 
#' @return nothing is returned.  file exported as csf to 'datasets/*.csv'
#' @concept ramclustR
#' @concept RAMClustR
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @export 

exportDataset<-function(
  ramclustObj=NULL,
  which.data="SpecAbund",
  label.by="ann",
  appendFactors = TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  d <- RAMClustR::getData(
    ramclustObj = ramclustObj, 
    which.data = which.data, 
    cmpdlabel = label.by
  )
  
  # if(filter){
  #   if(!is.null(ramclustObj$cmpd.use)) {
  #     if(length(ramclustObj$cmpd.use == ncol(d[[2]]))) {
  #       cmpd.use <- which(ramclustObj$cmpd.use)
  #     } else {
  #       cmpd.use <- 1:ncol(d[[2]])
  #     }
  #   } else {
  #     cmpd.use <- 1:ncol(d[[2]])
  #   }
  # } else {
  #   cmpd.use <- 1:ncol(d[[2]])
  # }
  check.names = FALSE
  dimnames(d[[2]])[[2]] <- ramclustObj[[label.by]]  
  dimnames(d[[3]])[[2]][((dim(d[[3]])[2] - dim(d[[2]])[2] + 1)):(dim(d[[3]])[2])] <- ramclustObj[[label.by]]
  if(appendFactors) {
    write.csv(d[[3]], paste("datasets/", which.data, ".csv", sep=""), row.names=FALSE)
  } else {
    write.csv(d[[2]], paste("datasets/", which.data, ".csv", sep=""), row.names=TRUE)
  }
}


