#' impRamSearch
#'
#' import ramsearch output for annotating an RC object
#' @details Annotation of ramclustR exported .msp spectra is accomplished using RAMSearch.  Exported ramsearch annotations (.rse) can be imported with this function
#' 
#' @param ramclustObj ramclustR object to annotate
#' @param ramsearchout path to .rse file to import
#' @return returns a ramclustR object.  new slots holding .rse data 
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept RAMSearch
#' @concept xcms
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @export 


impRamSearch<-function(
  ramclustObj=NULL,
  ramsearchout="spectra/results.rse"
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  out<-readLines(ramsearchout)
  ramclustObj$rs.out <- out
  
  ramclustObj$history$ramsearch <- paste0(
    "Annotation was performed in RAMSearch (Broeckling 2016), and annotations imported in the the ramclustR object. "
  )
  
  return(ramclustObj)
}


