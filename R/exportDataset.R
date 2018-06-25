#' exportDataset
#'
#' export one of 'SpecAbund', 'SpecAbundAve', 'MSdata' or 'MSMSdata' from an RC object to csv
#' @details Useful for exporting the processed signal intensity matrix to csv for analysis elsewhere.  
#' @param ramclustObj ramclustR object to export from
#' @param which.data name of dataset to export.  SpecAbund, SpecAbundAve, MSdata, or MSMSdata
#' @param label.by either 'ann' or 'cmpd', generally.  name of ramclustObj slot used as csv header for each column (compound)
#' @return nothing is returned.  file exported as csf to 'datasets/*.csv'
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'ramsearch'
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @export 

exportDataset<-function(ramclustObj=RC,
                        which.data="SpecAbund",
                        label.by="ann") {
  temp<-ramclustObj[[which.data]]
  if(!is.null(label.by)) dimnames(temp)[[2]]<-ramclustObj[[label.by]]  
  write.csv(temp, paste("datasets/", which.data, ".csv", sep=""), row.names=TRUE)
}


