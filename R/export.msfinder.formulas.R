#' export MSFinder formula prediction results in tabular format.  
#'
#' After running MSFinder, results have been imported to the ramclustR object.  This function exports as a .csv file for ease of viewing.
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param export.all logical: default = FALSE.  If TRUE, export all columns, if FALSE, only columns 1: "exactmass"
#' @param output.directory valid path: default = NULL.  If NULL, results are exported to spectra/mat directory. 
#' @details this function exports a .csv file containing all returned MSFinder molecular formula hypotheses. this file is saved (by default) to the working directory spectra/mat/ directory
#' @return an updated ramclustR object, with the RC$ann and RC$ann.conf slots updated to annotated based on output from 1. ramsearch output, 2. msfinder mssearch, 3. msfinder predicted structure, 4. msfinder predicted formula, and 5. interpretMSSpectrum inferred molecular weight, with listed order as priority.  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
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

export.msfinder.formulas <- function(
  ramclustObj = NULL, 
  export.all = FALSE,
  output.directory = NULL) {
  
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(!is.null(output.directory)) {
    if(!dir.exists(output.directory)) {
      stop("output directory ", output.directory, " does not exist")
    }
    output.directory = paste0(output.directory,"/spectra/mat/")
    if(!dir.exists(output.directory)) {
      stop("output directory ", output.directory, " does not exist")
    }
  } else {
    output.directory = paste0(getwd(),"/spectra/mat/")
    if(!dir.exists(output.directory)) {
      stop("output directory ", output.directory, " does not exist")
    }
  }
  
  ### now many rows in MSFinder formula results?
  nrf <- sapply(1:length(ramclustObj$msfinder.formula.details), 
                FUN = function(x) {
                  nrow(ramclustObj$msfinder.formula.details[[x]])
                })
  ncf <- sapply(1:length(ramclustObj$msfinder.formula.details), 
                FUN = function(x) {
                  ncol(ramclustObj$msfinder.formula.details[[x]])
                })
  
  if(sd(ncf)>0) {
    stop("there is a variable number of columns in the MSFinder output - something is wrong" , '\n', 
         "try reimporting MSFinder results... maybe?")
  }
  
  
  out <- data.frame(matrix(nrow = sum(nrf), ncol = (ncf[1]+1)))
  
  for(i in 1:length(nrf)) {
    
    if(i ==1) {
      dimnames(out)[[2]] <- c("cmpd", dimnames(ramclustObj$msfinder.formula.details[[i]])[[2]])
    }
    
    if(nrf[i] == 0) {next}
    
    start.row <- which(is.na(out[,1]))[1]
    tmp <- cbind(rep(as.character(ramclustObj$cmpd[i]), nrf[i]), 
                 ramclustObj$msfinder.formula.details[[i]], stringsAsFactors = FALSE)
    out[start.row:(start.row + nrf[i]-1), ] <- tmp
  }
  
  if(export.all) {
    cols <- 1:ncf[1]
  } else {
    col.end <- grep("exactmass", dimnames(out)[[2]])
    if(length(col.end) == 0) {
      col.end <- 5
    }
    cols <- c(1:col.end)
  }
  
  write.csv(out[,cols], file = paste0(output.directory, "/", "MSFinderFormulas_allResults.csv"), row.names = FALSE)

}

