#' write.gcei.mat
#'
#' Export GC-MS EI spectra for spectral searching in MSFinder
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param out.dir valid directory path describing output directory/file location.
#' @details exports files to a directory called 'spectra'.  a new directory 'spectra/mat' is created to hold the individual mat files.  
#' @return nothing, just exports files to the working directory
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept interpretMSSpectrum
#' @concept xcms
#' @author Corey Broeckling
#' @importFrom methods is
#' @export 
#' 

write.gcei.mat <- function(
    ramclustObj = NULL,
    out.dir = NULL
) {
  
  if(is.null(out.dir)) {
    stop("please provide a valid output directory")
  }
  
  
  if(!is(ramclustObj, "hclus") & 
     ramclustObj$dist.method != "RAMClustR") {
    stop("this is not a RAMClustR object")
  }
  
  if(!dir.exists(paste0(out.dir, '/spectra'))) {
    dir.create(paste0(out.dir, '/spectra'))
  }
  

    if(!dir.exists(paste0(out.dir, '/spectra/mat'))) {
      dir.create(paste0(out.dir, '/spectra/mat'))
    }

  ion.mode <- "Positive"
  
  out.list <- as.list(rep(NA, length(ramclustObj$cmpd)))
  
  for(i in 1:length(ramclustObj$cmpd)) {
    ions <- which(ramclustObj$featclus == i)
    
    spectrum <- data.frame(
      'mz' = round(0.99888*ramclustObj$fmz[ions], digits = 0),
      'int' = round(ramclustObj$msint[ions], digits = 0)
    )
    
    spectrum <- spectrum[order(spectrum[,"mz"], decreasing = FALSE), ]
    
    out <- paste0(
      "NAME: ", ramclustObj$cmpd[i], '\n', 
      "PRECURSORTYPE: [M]+.", '\n',
      "IONMODE: ", ion.mode, '\n',
      "COLLISIONENERGY: 70", '\n',
      "SPECTRUMTYPE: Centroid", '\n',
      "RETENTIONTIME: ", round(ramclustObj$clrt[i], 2), '\n'
    )
    
    if(any(names(ramclustObj)=="clri")) {paste0(
      out,
      "RETENTIONINDEX: ", round(ramclustObj$clri[i], 2),  '\n'
    )
    }
    
    out <- paste0(out, 
                  "PRECURSORMZ: ", round(max(spectrum[,1]),1), '\n' 
    )
    
    out <- paste0(
      out,
      "INSTRUMENTTYPE: GC-EI-Q", '\n',  
      "INSTRUMENT: ", '\n',
      "Authors: ", '\n', 
      "License: ", '\n',
      "FORMULA: ", '\n',
      "ONTOLOGY: ", '\n',
      "SMILES: ", '\n',  
      "INCHIKEY: ", '\n',
      "INCHI: ", '\n', 
      "METABOLITENAME: ", '\n',
      "SCANNUMBER: -1 ", '\n',
      "RETENTIONTIME: 0", '\n',
      "RETENTIONINDEX: 0", '\n',
      "CCS: 0", '\n',
      "INTENSITY: 0", '\n',
      "#Specific field for labeled experiment", '\n',
      "IsMarked: False", '\n'
    )
    
    out <- paste0(out,
                  "MSTYPE: MS1", '\n', "Num Peaks: 0", '\n'
    )
    
    out <- paste0(out,
                  "MSTYPE: MS2", '\n',
                  "Num Peaks: ", nrow(spectrum), '\n'
    )
    for(j in 1:nrow(spectrum)) {
      out <- paste0(out, 
                    spectrum[j,"mz"], 
                    " ", 
                    round(spectrum[j,"int"]),
                    '\n'
      )
    }
    out.list[[i]] <- out
  }
  
  for(i in 1:length(out.list)) {
    writeLines(out.list[[i]], '\n', con = paste0(out.dir, "/spectra/mat/", ramclustObj$cmpd[[i]], ".mat"))
  }

}
