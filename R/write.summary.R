#' write.summary()
#'
#' Write a .csv file containing a summary of the annotations in the ramclustR object.
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @details This function writes a summary .csv file in the 'spectra' directory, containing details on the annotations for all compounds. 
#' @return nothing is returned, run this without assigning a new object  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Jaeger C, Méret M, Schmitt CA, Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017 Aug 15;31(15):1261-1266. doi: 10.1002/rcm.7905. PubMed PMID: 28499062.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @references http://cts.fiehnlab.ucdavis.edu/static/download/CTS2-MS2015.pdf 
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms', 'MSFinder', 'chemical translation service', 'cts'
#' @author Corey Broeckling
#' @export 


write.summary<-function(ramclustObj = RC
                        ) {
  
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikeys found in this ramclustR object")
  }
  
  write.csv(file="spectra/annotationSummary.csv", data.frame("cmpd" = RC$cmpd,
                                                              "rt" = RC$clrt,
                                                              "annotation" = RC$ann,
                                                              "inchikey" = RC$inchikey,
                                                              "smiles" = RC$smiles,
                                                              "DB_accessions" = RC$dbid,
                                                              "synonyms" = sapply(1:length(RC$synonyms), function(x) {
                                                                paste(RC$synonyms[[x]], collapse =  " ;  ")
                                                              }),
                                                              "ann.confidence" = RC$annconf,
                                                              "inferred M" = RC$M,
                                                              "inferred formula" = RC$msfinder.formula), row.names = FALSE)
  cat("Annotation summary written to 'spectra' directory", '\n')
  
}


