#' import.msfinder.formulas
#'
#' After running MSFinder on .mat or .msp files, import the formulas that were predicted and their scores 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir optional path to .mat directory
#' @param msp.dir optional path to .msp directory
#' @details this function imports the output from the MSFinder program to support annotation of the ramclustR object
#' @return new slot at $msfinder.formula.details
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

import.msfinder.formulas <- function (ramclustObj = NULL, 
                                      mat.dir = NULL, 
                                      msp.dir = NULL) 
{
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  home.dir <- getwd()
  
  r <- grep("msfinder.formula", names(ramclustObj))
  if (length(r) > 0) {
    warning("removed previously assigned MSFinder formulas and structures", 
            "\n")
    ramclustObj <- ramclustObj[-r]
    r <- grep("msfinder.structure", names(ramclustObj))
    if(length(r)>0) {
      ramclustObj <- ramclustObj[-r]
    }
    rm(r)
  }
  if (is.null(mat.dir)) {
    mat.dir = paste0(getwd(), "/spectra/mat")
  }
  if (is.null(msp.dir)) {
    msp.dir = paste0(getwd(), "/spectra/msp")
  }
  usemat = TRUE
  usemsp = TRUE
  if (!dir.exists(mat.dir)) {
    usemat = FALSE
  }
  if (!dir.exists(msp.dir)) {
    usemsp = FALSE
  }
  if (!usemsp & !usemat) {
    stop("neither of these two directories exist: ", "\n", 
         "  ", mat.dir, "\n", "  ", msp.dir, "\n")
  }
  if (usemsp & usemat) {
    msps <- list.files(msp.dir, recursive = TRUE)
    mats <- list.files(mat.dir, recursive = TRUE)
    if (length(mats) > length(msps)) {
      usemsp <- FALSE
    }
    if (length(msps) > length(mats)) {
      usemat <- FALSE
    }
    if (length(msps) == length(mats)) {
      feedback <- readline(prompt = "Press 1 for .mat or 2 for .msp to continue")
      if (feedback == 1) {
        usemsp <- FALSE
      }
      if (feedback == 2) {
        usemat <- FALSE
      }
    }
  }
  mat.dir <- c(mat.dir, msp.dir)[c(usemat, usemsp)]
  do <- list.files(mat.dir, pattern = ".fgt", full.names = TRUE)
  
  ### retrieve parameter file from mat directory and parse to save with results.
  params <- list.files(mat.dir, pattern = "batchparam", full.names = TRUE)
  if(length(params) == 0) {
    params <- list.files(mat.dir, pattern = "MSFINDER.INI", full.names = TRUE)
  }
  if(length(params) > 0) {
    mtime <- rep(NA, length(params))
    for(i in 1:length(mtime)) {
      mtime[i] <- format(file.info(params[i])$mtime, format = '%y%m%d%H%M%S')
    }
    params <- params[which.max(mtime)]
    params <- readLines(params)
    breaks <- which(nchar(params)==0)
    
    ## formula inference parameters
    st <- grep ("Formula finder parameters", params)+1
    end <- breaks[which(breaks > st)[1]]-1
    if(end <= st) {stop('parsing of parameter file has failed')}
    tmp <- strsplit(params[st:end], "=")
    nms <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][1]})
    vals <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][2]})
    names(vals) <- nms
    ramclustObj$msfinder.formula.parameters <- vals
    
    ## DB used record
    st <- grep ("Data source", params)+1
    end <- breaks[which(breaks > st)[1]]-1
    if(end <= st) {stop('parsing of parameter file has failed')}
    tmp <- strsplit(params[st:end], "=")
    nms <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][1]})
    vals <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][2]})
    names(vals) <- nms
    if(grepl("F", vals["IsUserDefinedDB"])) {
      vals <- vals[1:(length(vals)-1)]
      nms <- nms[1:(length(nms)-1)]
    }
    vals <- as.logical(vals)
    names(vals) <- nms
    vals <- names(which(vals))
    vals <- vals[!grepl("NeverUse", vals)]
    vals <- gsub("OnlyUseForNecessary", "", vals)
    vals <- gsub("Allways", "", vals)
    vals <- unique(vals)
    ramclustObj$msfinder.dbs <- vals
  } 
    
  
  
  cmpd <- gsub(".fgt", "", basename(do))
  specres <- 0
  allres <- 0
  for(i in 1:min(10,length(do))) {
    tmp <- readLines(do[[i]])
    allres <- allres + length(which(grepl("NAME:", tmp, ignore.case = TRUE)))
    specres <- specres + length(which(grepl("NAME: Spectral DB search", tmp, ignore.case = TRUE)))
    rm(tmp)
  }
  if (specres == allres) {
    stop("these MSFinder contain only spectral search results; please use 'import.MSFinder.search' function instead", 
         "\n")
  }
  tags <- c("NAME: ", "EXACTMASS: ", "ISSELECTED: ", "MASSDIFFERENCE: ", 
            "TOTALSCORE: ", "ISOTOPICINTENSITY[M+1]: ", "ISOTOPICINTENSITY[M+2]: ", 
            "ISOTOPICDIFF[M+1]: ", "ISOTOPICDIFF[M+2]: ", "MASSDIFFSCORE: ", 
            "ISOTOPICSCORE: ", "PRODUCTIONSCORE: ", "NEUTRALLOSSSCORE: ", 
            "PRODUCTIONPEAKNUMBER: ", "PRODUCTIONHITSNUMBER: ", 
            "NEUTRALLOSSPEAKNUMBER: ", "NEUTRALLOSSHITSNUMBER: ", 
            "RESOURCENAMES: ", "RESOURCERECORDS: ", "ChemOntDescriptions: ", 
            "ChemOntIDs: ", "ChemOntScores: ", "ChemOntInChIKeys: ", 
            "PUBCHEMCIDS: ")
  names(tags) <- tolower(gsub(": ", "", tags))
  fill <- matrix(nrow = 0, ncol = length(tags))
  dimnames(fill)[[2]] <- names(tags)
  fill <- data.frame(fill, check.names = FALSE, stringsAsFactors = FALSE)
  msfinder.formula <- as.list(rep(NA, length(ramclustObj$cmpd)))
  names(msfinder.formula) <- ramclustObj$cmpd
  for (i in 1:length(do)) {
    tmp <- readLines(do[[i]])
    starts <- grep("NAME: ", tmp)
    if (length(starts) < 1) {
      msfinder.formula[[cmpd[i]]] <- fill
      next
    }
    if (length(starts) > 1) {
      stops <- c(((starts[2:length(starts)]) - 1), (length(tmp) - 
                                                      1))
    }
    else {
      stops <- length(tmp) - 1
    }
    out <- matrix(nrow = length(starts), ncol = length(tags))
    dimnames(out)[[2]] <- names(tags)
    for (j in 1:length(starts)) {
      d <- tmp[starts[j]:stops[j]]
      vals <- sapply(1:length(tags), FUN = function(x) {
        m <- grep(tags[x], d, fixed = TRUE)
        if (length(m) == 0) {
          NA
        }
        else {
          if (length(m) > 1) {
            m <- m[1]
          }
          gsub(tags[x], "", d[m])
        }
      })
      out[j, ] <- vals
    }
    out <- data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
    if (any(out[, "name"] == "Spectral DB search")) {
      out <- out[-grep("Spectral DB search", out[, "name"]), 
                 , drop = FALSE]
    }
    if (nrow(out) == 0) {
      msfinder.formula[[cmpd[i]]] <- fill
    }
    else {
      msfinder.formula[[cmpd[i]]] <- out
    }
  }
  
  missing <- which(is.na(msfinder.formula))
  if(length(missing > 0)) {
    for(x in missing) {
      msfinder.formula[[x]] <- fill
    }
  }
  
  ramclustObj$msfinder.formula.details <- msfinder.formula
  setwd(home.dir)
  
  if(is.null(ramclustObj$history)) {
    ramclustObj$history <- ""
  }
  
  ramclustObj$history$msfinder <- paste(
    "MSFinder (Tsugawa 2016) was used for spectral matching,",
    "formula inference, and tentative structure assignment,",
    "and results were imported into the RAMClustR object.")
  
  if(is.null(ramclustObj$msfinder.dbs)) {
    dbs <- sapply(1:length(ramclustObj$ann), 
                          FUN = function(x) {
                            tmp <- paste(ramclustObj$msfinder.formula.details[[x]]$resourcenames, collapse = ",")
                            tmp <- strsplit(tmp, ",", fixed = TRUE)
                            tmp <- tmp[which(nchar(tmp)>0)]
                            return(tmp)
                          }
    )
    dbs <- unique(unlist(dbs))
    dbs <- dbs[which(nchar(dbs)>0)]
    ramclustObj$msfinder.dbs <- dbs
  }
  
  return(ramclustObj)
}
