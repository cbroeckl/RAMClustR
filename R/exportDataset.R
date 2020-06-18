#' exportDataset
#'
#' export one of 'SpecAbund', 'SpecAbundAve', 'MSdata' or 'MSMSdata' from an RC object to csv
#' @details Useful for exporting the processed signal intensity matrix to csv for analysis elsewhere.  
#' @param ramclustObj ramclustR object to export from
#' @param which.data name of dataset to export.  SpecAbund, SpecAbundAve, MSdata, or MSMSdata
#' @param label.by either 'ann' or 'cmpd', generally.  name of ramclustObj slot used as csv header for each column (compound)
#' @param filter logical, TRUE by default. when $cmpd.use slot is present (from rc.cmpd.filter.cv function), only cmpds that passed cv filtering are used. If you wish to change that behavior, rerun the rc.cmpd.filter.cv function with a really high CV threshold. 
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
  filter = TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  d <- getData(
    ramclustObj = ramclustObj, 
    which.data = which.data, 
    filter = filter,
    cmpdlabel = label.by
  )
  
  if(filter){
    if(!is.null(ramclustObj$cmpd.use)) {
      if(length(ramclustObj$cmpd.use == ncol(d[[2]]))) {
        cmpd.use <- which(ramclustObj$cmpd.use)
      } else {
        cmpd.use <- 1:ncol(d[[2]])
      }
    } else {
      cmpd.use <- 1:ncol(d[[2]])
    }
  } else {
    cmpd.use <- 1:ncol(d[[2]])
  }
  
  dimnames(d[[2]])[[2]]<-ramclustObj[[label.by]][cmpd.use]  
  write.csv(d, paste("datasets/", which.data, ".csv", sep=""), row.names=TRUE)
}


