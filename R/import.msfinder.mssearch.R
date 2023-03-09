#' import.MSFinder.mssearch
#'
#' After running MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir optional path to .mat directory
#' @param msp.dir optional path to .msp directory
#' @param dir.extension optional directory name code specifying subset of results to use.  Useful if running MSFinder from the command line for both spectral searching and interpretation. 
#' @details this function imports the output from the MSFinder program to annotate the ramclustR object
#' @return an updated ramclustR object, with new slots at $msfinder.mssearch.details and $msfinder.mssearch.scores
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

import.msfinder.mssearch <- function (
  ramclustObj = NULL,
  mat.dir = NULL,
  msp.dir = NULL,
  dir.extension = ".mssearch"
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  home.dir <-getwd()
  
  if(is.null(mat.dir)) {
    mat.dir = paste0(getwd(), "/spectra/mat")
  }
  
  if(is.null(msp.dir)) {
    msp.dir = paste0(getwd(), "/spectra/msp")
  }
  
  usemat = TRUE
  usemsp = TRUE
  
  if(!dir.exists(mat.dir)) {
    usemat = FALSE
  }
  
  if(!dir.exists(msp.dir)) {
    usemsp = FALSE
  }
  
  if(!usemsp & !usemat) {
    stop("neither of these two directories exist: ", '\n',
         "  ", mat.dir, '\n',
         "  ", msp.dir, '\n')
  }
  
  if(usemsp & usemat) {
    msps<-list.files(msp.dir, recursive  = TRUE)
    mats<-list.files(mat.dir, recursive = TRUE)
    if(length(mats) > length(msps)) {usemsp <- FALSE}
    if(length(msps) > length(mats)) {usemat <- FALSE}
    if(length(msps) == length(mats)) {
      feedback<-readline(prompt="Press 1 for .mat or 2 for .msp to continue")
      if(feedback == 1) {usemsp <- FALSE}
      if(feedback == 2) {usemat <- FALSE}
    }
  } 
  
  mat.dir <- c(mat.dir, msp.dir)[c(usemat, usemsp)]
  
  do<-list.dirs(mat.dir, recursive = FALSE, full.names = FALSE)
  cmpds <- do
  if(!is.null(dir.extension)) {
    keep <- grep(dir.extension, do)
    do<-do[keep] 
    cmpds <- cmpds[keep]
    cmpds <- gsub(dir.extension, "", cmpds)
  }
  
  if(length(cmpds) == 0) {
    stop("  --  no spectral match results to return")
  }

  
  ### retrieve parameter file from mat directory and parse to save with results.
  params <- list.files(mat.dir, pattern = "batchparam", full.names = TRUE)
  if(length(params) > 0) {
    mtime <- rep(NA, length(params))
    for(i in 1:length(mtime)) {
      mtime[i] <- format(file.info(params[i])$mtime, format = '%y%m%d%H%M%S')
    }
    params <- params[which.max(mtime)]
    params <- readLines(params)
    breaks <- which(nchar(params)==0)
    
    ## spectral search parameters
    st <- grep ("Spectral database search", params)+1
    end <- breaks[which(breaks > st)[1]]-1
    if(end <= st) {stop('parsing of parameter file has failed')}
    tmp <- strsplit(params[st:end], "=")
    nms <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][1]})
    vals <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][2]})
    names(vals) <- nms
    ramclustObj$msfinder.mssearch.parameters <- vals
    
  } 
  
  
  ncmpd <- max(as.integer(gsub("C", "", cmpds)))
  
  tags<-c(
    "NAME: ",
    "ID: ",
    "IsSpectrumSearch: ",
    "INCHIKEY: ",
    "SMILES: ",
    "RESOURCES: ",
    "SubstructureInChIKeys: ",
    "RETENTIONTIME: ",
    "RETENTIONINDEX: ",
    "TotalBondEnergy: ",
    "TotalScore: ",
    "TotalHrRulesScore: ",
    "TotalBondCleavageScore: ",
    "TotalMassAccuracyScore: ",
    "TotalFragmentLinkageScore: ",
    "TotalBondDissociationEnergyScore: ",
    "DatabaseScore: ",
    "SubstructureAssignmentScore: ",
    "RtSimilarityScore: ",
    "RiSimilarityScore: ",
    "Num Peaks: "
  )
  
  names(tags)<-tolower(gsub(": ", "", tags))
  
  msfinder.mssearch<-as.list(rep(NA, length(ramclustObj$cmpd)))
  names(msfinder.mssearch)<-ramclustObj$cmpd
  msfinder.mssearch.details<-as.list(rep(NA, length(ramclustObj$cmpd)))
  names(msfinder.mssearch.details)<-ramclustObj$cmpd
  template<-as.list(c("summary" = NA, "spectra" = NA))
  for(i in 1:length(msfinder.mssearch)) {
    msfinder.mssearch[[i]]<-template
    msfinder.mssearch.details[[i]]<-template
  }
  
  for(i in 1:length(ramclustObj$cmpd)) {
    tar.dir <- paste0(mat.dir, "/", ramclustObj$cmpd[i])
    if(!is.null(dir.extension)) tar.dir <- paste0(tar.dir, dir.extension)
    if(!dir.exists(tar.dir)) next
    setwd(tar.dir)
    if(!file.exists("Spectral DB search.sfd")) {next}
    tmp<-readLines("Spectral DB search.sfd")
    if(length(tmp) == 0) {next}
    starts <- grep("NAME: ", tmp)
    if(length(starts) > 1) {
      stops <- c(((starts[2:length(starts)])-1), (length(tmp)-1))
    } else {
      stops <- length(tmp)-1
    }
    
    out<-matrix(nrow = length(starts), ncol = length(tags))
    dimnames(out)[[2]]<-names(tags)
    spectra<-as.list(rep(NA, length(starts)))
    for(j in 1:length(starts)) {
      d<-tmp[starts[j]:stops[j]]
      vals<-sapply(1:length(tags), FUN = function(x) {
        m<-grep(tags[x], d, fixed = TRUE)
        if(length(m)==0) {
          NA
        } else {
          if(length(m)>1) {
            m<-m[1]
          }
          gsub(tags[x], "", d[m])
        }
      }
      )
      out[j,]<-vals
      spec<-d[(grep("Num Peaks: ", d)+1):length(d)]
      if(any(nchar(spec) == 0)) { 
        spec <- spec[-(which(nchar(spec) == 0))]
      }
      if(length(spec) > 0) {
        spec<-data.frame(matrix(unlist(strsplit(spec, "\t")) , nrow = length(spec), byrow = TRUE), stringsAsFactors = FALSE)
        names(spec)<-c("mz", "int")
        for(k in 1:2) {
          spec[,k]<-as.numeric(spec[,k])
        }
        spectra[[j]]<-(spec)
      } else {
        spec<-data.frame("mz" = NA, "int" = NA)
        spec<-spec[0,]
        spectra[[j]]<-spec
      }
      
      
    }
    out<-data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
    if(any(out$name == "Spectral DB search")) {
      out <- out[-grep("Spectral DB search", out$name),]
    }
    for(k in 1:ncol(out)) {
      suppressWarnings(an<-as.numeric(out[,k]))
      if(is.na(an[1])) {next}
      if(is.numeric(an[1])) {
        out[,k]<-an
      }
    }
    msfinder.mssearch.details[[i]]$summary<-out
    msfinder.mssearch.details[[ramclustObj$cmpd[i]]]$spectra<-spectra
  }
  
  
  for(i in 1:length(msfinder.mssearch.details ) ){
    if(is.null(nrow(msfinder.mssearch.details[[i]]$summary))) {
      msfinder.mssearch.details[[i]]$summary <- data.frame(t(tags))[0,]
    }
  }
  
  
  ramclustObj$msfinder.mssearch.details<-msfinder.mssearch.details
  setwd(home.dir)
  
  ramclustObj$history$msfinder <- paste(
    "MSFinder (Tsugawa 2016) was used for spectral matching,",
    "formula inference, and tentative structure assignment,",
    "and results were imported into the RAMClustR object.")
  
  return(ramclustObj)
  
}

