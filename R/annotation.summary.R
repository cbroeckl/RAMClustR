#' annotation.summary()
#'
#' Write a .csv file containing a summary of the annotations in the ramclustR object. 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param outfile file path/name of output csv summary file.  if NULL (default) will be exported to spectra/annotaionSummary.csv
#' @details this function exports a csv file summarizing annotation evidence for each compound
#' @return nothing
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @concept ramlclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept xcms
#' @concept MSFinder
#' @author Corey Broeckling
#' @export 


annotation.summary<-function(ramclustObj = NULL,
                             outfile = NULL
) {
  
  if(!is.null(outfile)) {
    f<-basename(outfile)
    p<-dirname(outfile)
    if(!dir.exists(p)) {
      dir.create(p)
    }
  } else {
    outfile <- paste0(getwd(), "/spectra/annotationSummary.csv")
  }
  
  out<- data.frame("cmpd" = ramclustObj$cmpd,
                   "rt" = ramclustObj$clrt,
                   "annotation" = ramclustObj$ann,
                   "ann.confidence" = ramclustObj$annconf,
                   "median signal" = as.vector(apply(ramclustObj$SpecAbund, 2, "median"))) 
  
  if(any(names(ramclustObj) == "M")) {
    out<- data.frame(out, "inferred M" = ramclustObj$M)
  }
  if(any(names(ramclustObj) == "zmax")) {
    out<- data.frame(out, "zmax" = ramclustObj$zmax)
  }
  if(any(names(ramclustObj) == "msfinder.formula")) {
    out<- data.frame(out, "inferred formula" = ramclustObj$msfinder.formula)
  }
  if(any(names(ramclustObj) == "inchikey")) {
    out<- data.frame(out, "inchikey" = ramclustObj$inchikey)
  }
  if(any(names(ramclustObj) == "inchi")) {
    out<- data.frame(out, "inchi" = ramclustObj$inchi)
  }
  if(any(names(ramclustObj) == "hmdb.url")) {
    # out<- data.frame(out, "hmdb.url" = paste0("<a href='", ramclustObj$hmdb.url, "'>", ramclustObj$hmdb.name, "</a>")) 
    out<- data.frame(out, "hmdb.url" = ramclustObj$hmdb.url)
  }
  if(any(names(ramclustObj) == "lm.url")) {
    out<- data.frame(out, "lm.url" = ramclustObj$lm.url)
  }
  if(any(names(ramclustObj) == "pubchem.url")) {
    out<- data.frame(out, "pubchem.url" = ramclustObj$pubchem.url)
  }
  if(any(names(ramclustObj) == "chebi.url")) {
    out<- data.frame(out, "chebi.url" = ramclustObj$chebi.url)
  }
  if(any(names(ramclustObj) == "synonyms")) {
    out<- data.frame(out, "synonyms" = sapply(1:length(ramclustObj$synonyms), 
                                              FUN = function(x) {
                                                paste(ramclustObj$synonyms[[x]], collapse = " __ ")
                                              }))
  }
  if(any(names(ramclustObj) == "classyfire")) {
    out<- data.frame(out, ramclustObj$classyfire)
  }
  
  write.csv(out, file = outfile, row.names = FALSE)
}

