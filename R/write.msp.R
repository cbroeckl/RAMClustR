#' write.msp
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param one.file logical, should all msp spectra be written to one file? If false, each spectrum is an individual file.
#' @details exports files to a directory called 'spectra'.  If one.file = FALSE, a new directory 'spectra/msp' is created to hold the individual msp files. if do.findman has been run, spectra are written as ms2 spectra, else as ms1. 
#' @return nothing, just exports files to the working directory
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept interpretMSSpectrum
#' @concept xcms
#' @author Corey Broeckling
#' @export 
#' 

write.msp <- function(
  ramclustObj = NULL,
  one.file = FALSE
) {
  
  if(class(ramclustObj) != "hclus" & 
     ramclustObj$dist.method != "RAMClustR") {
    stop("this is not a RAMClustR object")
  }
  
  if(is.null(ramclustObj$precursor.mz)) {
    ms2 <- FALSE
  } else {
    ms2 <- TRUE
  }
  
  if(!dir.exists('spectra')) {
    dir.create('spectra')
  }
  
  if(!one.file) {
    if(!dir.exists('spectra/msp')) {
      dir.create('spectra/msp')
    }
  }
  ion.mode <- as.character(ramclustObj$ExpDes[[2]][which(row.names(ramclustObj$ExpDes[[2]]) == "msmode"),1])
  if(toupper(substring(ion.mode, 1, 1)) == "P") {
    ion.mode = "Positive"
  } else {
    ion.mode = "Negative"
  }
  
  out.list <- as.list(rep(NA, length(ramclustObj$cmpd)))
  
  for(i in 1:length(ramclustObj$cmpd)) {
    ions <- which(ramclustObj$featclus == i)
    
    if(ms2) {
      spectrum <- data.frame(
        'mz' = ramclustObj$fmz[ions],
        'int' = ramclustObj$msmsint[ions]
      )
    } else {
      spectrum <- data.frame(
        'mz' = ramclustObj$fmz[ions],
        'int' = ramclustObj$msint[ions]
      )
    }
    
    spectrum <- spectrum[order(spectrum[,"int"], decreasing = TRUE), ]
    
    out <- paste0(
      "NAME:", ramclustObj$cmpd[i], '\n', 
      "IONMODE:", ion.mode, '\n',
      "SPECTRUMTYPE:Centroid", '\n',
      "RETENTIONTIME:", round(ramclustObj$clrt[i], 2), '\n'
    )
    
    if(any(names(ramclustObj)=="clri")) {paste0(
      out,
      "RETENTIONINDEX:", round(ramclustObj$clri[i], 2),  '\n'
    )
    }
    
    if(ms2) {
      out <- paste0(out, 
                    "PRECURSORMZ:", ramclustObj$precursor.mz[i],'\n', 
                    "PRECURSORTYPE:", ramclustObj$precursor.type[i], '\n'
      )
    }
    out <- paste0(out,
                  "Num Peaks:", nrow(spectrum), '\n'
                  # m/z intensity pair (tab, comma, space can be used as the delimiter.)
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
  
  if(one.file) {
    out <- vector(mode = "character")
    for(i in 1:length(out.list)) {
      out <- paste0(out, out.list[[i]], '\n')
    }
    exp.name <- ramclustObj$ExpDes[[1]][which(row.names(ramclustObj$ExpDes[[1]]) == "Experiment"),1]
    if(nchar(exp.name) == 0) {
      exp.name <- "spectra"
    }
    sink(paste0("spectra/", exp.name, ".msp"))
    cat(out)
    sink()
  } else {
    for(i in 1:length(out.list)) {
      sink(paste0("spectra/msp/", ramclustObj$cmpd[[i]], ".msp"))
      cat(out.list[[i]], '\n')
      sink()
    }
  }
  
}
