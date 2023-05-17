#' rc.export.msp.rc
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param one.file logical, should all msp spectra be written to one file? If false, each spectrum is an individual file.
#' @param mzdec integer.  Number of decimal points to export mass values with.
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
#' @importFrom methods is
#' @export 
#' 
rc.export.msp.rc <- function(
  ramclustObj = NULL,
  one.file = TRUE,
  mzdec = 1
) {
  
  if(!is(ramclustObj, "hclus") & 
     ramclustObj$dist.method != "RAMClustR") {
    stop("this is not a RAMClustR object")
  }
  
  if(is.null(ramclustObj$precursor.mz)) {
    ms2 <- FALSE
  } else {
    ms2 <- TRUE
  }
  
  if(is.null(ramclustObj$MSMSdata)) {
    mslev <- 1
  } else {
    mslev <- 2
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
    ion.mode = "P"
  } else {
    ion.mode = "N"
  }
  
  ExpDes <- ramclustObj$ExpDes
  out.list <- as.list(rep(NA, length(ramclustObj$cmpd)*mslev))
  
  for (m in 1:as.numeric(mslev)){
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
      spectrum$mz <- round(spectrum$mz, digits = mzdec)
      spectrum$int <- round(spectrum$int, digits = 1)
      
      tmp <- paste0(
        "Name: ", ramclustObj$cmpd[i], '\n',
        "SYNON: $:00", if(!ms2){"in-source"} else {"MS2"}, '\n',
        "SYNON: $:05", if(m==1) {
          ExpDes[[2]]["CE1", 1]
        } else {
          ExpDes$instrument["CE2", "InstVals"]
        }, 
        '\n',
        "SYNON: $:06", ExpDes[[2]]["mstype", 1], '\n',           #mstype
        "SYNON: $:07", ExpDes[[2]]["msinst", 1], '\n',           #msinst
        "SYNON: $:09", ExpDes[[2]]["chrominst", 1], '\n',        #chrominst
        "SYNON: $:10", ExpDes[[2]]["ionization", 1],   '\n',      #ionization method
        "SYNON: $:11", ExpDes[[2]]["msmode", 1], '\n',           #msmode
        "SYNON: $:14", ExpDes[[2]]["msscanrange", 1], '\n',      #ms scanrange
        "Comment: Rt=", round(ramclustObj$clrt[i], digits=2), 
        "  Contributor=", ExpDes[[1]]["Contributor", 1], 
        "  Study=", ExpDes[[1]]["Experiment", 1], '\n',
        "Num Peaks:", nrow(spectrum), '\n',
        paste(sapply(1:nrow(spectrum), 
                     FUN = function(x) {
                       paste(spectrum$mz[x], spectrum$int[x])
                     }
        ), collapse = '\n'), '\n', '\n'
      )
      
      if(m == 1) {out.list[[i]] <- tmp}
      if(m == 2) {out.list[[i + length(ramclustObj$cmpd)]] <- tmp}
    }
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
    sink(paste0("spectra/", exp.name, ".rc.msp"))
    cat(out)
    sink()
  } else {
    for(i in 1:length(out.list)) {
      sink(paste0("spectra/msp/", ramclustObj$cmpd[[i]], ".rc.msp"))
      cat(out.list[[i]], '\n')
      sink()
    }
  }
  
}

