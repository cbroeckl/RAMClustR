#' write.methods
#'
#' write RAMClustR processing methods and citations to text file 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param filename define filename/path to write.  uses 'ramclustr_methods.txt' and the working directory by default.
#' @details this function exports a file called ramclustr_methods.txt which contains the processing history, parameters used, and relevant citations.
#' @return an annotated ramclustR object
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @return nothing - new file written to working director
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept xcms
#' @author Corey Broeckling
#' @export

write.methods <- function (
  ramclustObj = NULL,
  filename = NULL
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC1", '\n')
  }
  
  if(!any(names(ramclustObj) == "history")) {
    stop("no processing history present for this ramclustR object")
  }
  
  if(is.null(filename)) {
    filename <- 'ramclustr_methods.txt' 
  }
  
  cit.list <- c(
    'R Core Team' = paste0(
      citation()$author, 
      " (", citation()$year, "). ",
      citation()$title, ". ",
      citation()$organization, ", ", 
      citation()$address, ", ",
      citation()$url, "."
    ),
    
    '(Broeckling 2012)' = "Broeckling CD, Heuberger A, Prince JA, Ingelsson E, Prenni JE. Assigning precursor-product ion relationships in indiscriminant MS/MS data from non-targeted metabolite profiling studies. Metabolomics 2012. 9(1):33-43.",
    
    '(Broeckling 2014)' = "Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014. 86(14):6812-7.",
    
    '(Broeckling 2016)' = "Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016. 88(18):9226-34.", 
    
    '(Tsugawa 2016)' = "Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016. 88(16):7946-58.",
    
    '(Jaeger 2016)' = "Jaeger C, Meret M, Schmitt C, Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017. 31(15):1261-1266.",
    
    '(Wohlgemuth 2010)' = "Wohlgemuth G, Haldiya PK, Willighagen E, Kind T, Fiehn O. The Chemical Translation Service - a web-based tool to improve standardization of metabolomic reports. Bioinformatics. 2010. 26(20):2647-8.",
    
    '(Wishart 2018)' = "Wishart DS, Feunang YD, Marcu A, Guo AC, Liang K, Vazquez-Fresno R, Sajed T, Johnson D, Li C, Karu N, Sayeeda Z, Lo E, Assempour N, Berjanskii M, Singhal S, Arndt D, Liang Y, Badran H, Grant J, Serra-Cayuela A, Liu Y, Mandal R, Neveu V, Pon A, Knox C, Wilson M, Manach C, Scalbert A. HMDB 4.0: the human metabolome database for 2018. Nucleic Acids Res. 2018. 46(D1):D608-D617.",
    
    '(Fahy 2007)' = "Fahy E, Sud M, Cotter D, Subramaniam S. LIPID MAPS online tools for lipid research. Nucleic Acids Res. 2007. 35:W606-12.",
    
    '(Kim 2019)' = "Kim S, Chen J, Cheng T, Gindulyte A, He J, He S, Li Q, Shoemaker BA, Thiessen PA, Yu B, Zaslavsky L, Zhang J, Bolton EE. PubChem 2019 update: improved access to chemical data. Nucleic Acids Res. 2019. 47(D1):D1102-D1109.",
    
    '(Djoumbou 2016)' = "Djoumbou Feunang Y, Eisner R, Knox C, Chepelev L, Hastings J, Owen G, Fahy E, Steinbeck C, Subramanian S, Bolton E, Greiner R, Wishart DS. ClassyFire:  automated chemical classification with a comprehensive, computable taxonomy. J  Cheminform. 2016. 8:61.",
    
    '(Smith 2006)' = "Smith, C.A. and Want, E.J. and O'Maille, G. and Abagyan,R. and Siuzdak, G.: XCMS: Processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching and identification, Analytical Chemistry, 78:779-787 (2006)",
    
    '(Tautenhahn 2008)' = "Ralf Tautenhahn, Christoph Boettcher, Steffen Neumann: Highly sensitive feature detection for high resolution LC/MS BMC Bioinformatics, 9:504 (2008)"
    
  )

  names(cit.list)[which(names(cit.list) == "R")] <- paste0(
    citation()$author, 
    " (", citation()$year, ")")
  
  # paste0("(", citation()$author, " ",  citation()$year, ")") = paste0(citation()$author)
  sink(filename)
  
  history <- paste(ramclustObj$history, collapse = " " )
  cat(history)
  
  cites <- sapply(1:length(cit.list), FUN = function(x) {
    grepl(names(cit.list[x]), history)
  }
  )
  
  if(any(cites)) {
    cat('\n', '\n')
    cit.list <- cit.list[cites]
    cit.list <- sort(cit.list)
    for(i in 1:length(cit.list)) {
      cat(names(cit.list[i]), ":  ", cit.list[i], '\n', '\n', sep = "")
    }
  }
  
  if(grepl("R Core Team", history)) {
    cat(paste0(
      citation()$author, 
      " (", citation()$year, "). ",
      citation()$title, ". ",
      citation()$organization, ", ", 
      citation()$address, ", ",
      citation()$url, "."
    ))
  }
  sink()

}


