#' get.synonyms()
#'
#' Use chemical translation service to retreive synonyms for ramclustR inchikeys
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @details this function uses the chemical translation service (http://cts.fiehnlab.ucdavis.edu/) to look up known names (synonyms for a given inchikey)
#' @return an updated ramclustR object, with the RC$synonyms slot containing a list with length equal to the number of compounds in the dataset.  Each list element represents a character vector of all the synonyms returns (or NA, for compounds with no inchikey)  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references http://cts.fiehnlab.ucdavis.edu/static/download/CTS2-MS2015.pdf 
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms', 'MSFinder', 'chemical translation service', 'cts'
#' @author Corey Broeckling
#' @export 


get.synonyms<-function(ramclustObj = RC
                      ) {
  
  require(jsonlite)
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikeys found in this ramclustR object")
  }
  
  cat("using chemical translation service - requires interet access and may take a few minutes to complete", '\n')
  synonyms<-as.list(rep(NA, length(ramclustObj$ann)))
  names(synonyms)<-ramclustObj$cmpd
  
  require(jsonlite)
  for(i in 1:length(ramclustObj$ann)) {
    Sys.sleep(delay.time)
    if(!is.na(ramclustObj$inchikey[i])) {
      link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/synonyms/", ramclustObj$inchikey[i])
      suppressWarnings(out<-readLines(link))
      syns<-unlist(fromJSON(out))
      if(!is.null(syns)) {
        syns <- syns[order(nchar(syns))]
        synonyms[[i]] <- syns
      }
    }
  }
  
  ramclustObj$synonyms <- synonyms
  return(ramclustObj)
}


