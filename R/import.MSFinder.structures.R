#' write.methods
#'
#' write RAMClustR processing methods and citations to text file 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir directory in which to look for mat file MSFinder output - by default the /spectra/mat in the working directory
#' @param msp.dir directory in which to look for msp file MSFinder output - by default the /spectra/msp in the working directory
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
#' @concept MSFinder
#' @concept xcms
#' @concept MSFinder
#' @author Corey Broeckling
#' @export

import.msfinder.structures <- function (
  ramclustObj = NULL,
  mat.dir = NULL,
  msp.dir = NULL
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  home.dir <-getwd()
  
  if(is.null(mat.dir)) {
    mat.dir = paste0(home.dir, "/spectra/mat")
  }
  
  if(is.null(msp.dir)) {
    msp.dir = paste0(home.dir, "/spectra/msp")
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
  
  do <-  ramclustObj$cmpd
  
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
    
    ## structure inference parameters
    st <- grep ("Structure finder parameters", params)+1
    end <- breaks[which(breaks > st)[1]]-1
    if(end <= st) {stop('parsing of parameter file has failed')}
    tmp <- strsplit(params[st:end], "=")
    nms <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][1]})
    vals <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][2]})
    names(vals) <- nms
    ramclustObj$msfinder.structure.parameters <- vals
    
    # ## DB used record  THIS IS ALREADY DONE IN FORMULA IMPORT
    # st <- grep ("Data source", params)+1
    # end <- breaks[which(breaks > st)[1]]-1
    # if(end <= st) {stop('parsing of parameter file has failed')}
    # tmp <- strsplit(params[st:end], "=")
    # nms <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][1]})
    # vals <- sapply(1:length(tmp), FUN = function(x) {tmp[[x]][2]})
    # names(vals) <- nms
    # if(grepl("F", vals["IsUserDefinedDB"])) {
    #   vals <- vals[1:(length(vals)-1)]
    #   nms <- nms[1:(length(nms)-1)]
    # }
    # vals <- as.logical(vals)
    # names(vals) <- nms
    # vals <- names(which(vals))
    # vals <- vals[!grepl("NeverUse", vals)]
    # vals <- gsub("OnlyUseForNecessary", "", vals)
    # vals <- gsub("Allways", "", vals)
    # vals <- unique(vals)
    # ramclustObj$msfinder.formula.dbs <- vals
    
  } 
  
  
  
  
  
  # 
  # if(grepl("Spectral DB search", tmp[2])) {
  #   stop("these MSFinder contain only spectral search results; please use 'import.MSFinder.search' function instead", '\n')
  # }
  # 
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
    "RiSimilarityScore: "
  )
  
  
  names(tags)<-tolower(gsub(": ", "", tags))
  
  fill<-matrix(nrow = 0, ncol = length(tags))
  dimnames(fill)[[2]]<-names(tags)
  fill<-data.frame(fill, check.names = FALSE, stringsAsFactors = FALSE)
  
  msfinder.structure<-as.list(rep("", length(ramclustObj$cmpd)))
  names(msfinder.structure)<-ramclustObj$cmpd
  
  fillfrag<-structure(list(`M/Z` = character(0), Intensity = character(0), 
                           MatchedExactMass = character(0), SaturatedExactMass = character(0), 
                           Formula = character(0), RearrangedHydrogen = character(0), 
                           PPM = character(0), MassDiff_mDa = character(0), IsEvenElectronRule = character(0), 
                           IsHrRule = character(0), IsSolventAdductFragment = character(0), 
                           AssignedAdductMass = character(0), AdductString = character(0), 
                           BondDissociationEnergy = character(0), TreeDepth = character(0), 
                           SMILES = character(0), TotalScore = character(0), HrLikelihood = character(0), 
                           BcLikelihood = character(0), MaLikelihood = character(0), 
                           FlLikelihood = character(0), BeLikelihood = character(0)), .Names = c("M/Z", 
                                                                                                 "Intensity", "MatchedExactMass", "SaturatedExactMass", "Formula", 
                                                                                                 "RearrangedHydrogen", "PPM", "MassDiff_mDa", "IsEvenElectronRule", 
                                                                                                 "IsHrRule", "IsSolventAdductFragment", "AssignedAdductMass", 
                                                                                                 "AdductString", "BondDissociationEnergy", "TreeDepth", "SMILES", 
                                                                                                 "TotalScore", "HrLikelihood", "BcLikelihood", "MaLikelihood", 
                                                                                                 "FlLikelihood", "BeLikelihood"), row.names = integer(0), class = "data.frame")
  
  for(i in 1:length(do)) {
    if(!dir.exists(paste0(mat.dir, "/", do[i]))) {next}
    setwd(paste0(mat.dir, "/", do[i]))
    res<-list.files(pattern = ".sfd")
    if(length(res) == 0) {next}
    structures<-as.list(rep(NA, length(res)))
    # for(j in 1:length(formulas)) {
    #   formulas[[j]]<-c('structures' = "structures", 'fragments' =  "fragments")
    # }
    names(structures)<-gsub(".sfd", "", res)
    for(j in 1:length(structures)) {
      structures[[j]]<-as.list(rep(NA, 2))
      names(structures[[j]])<-c("structures", "fragments")
    }
    
    for(j in 1:length(structures)) {
      tmp<-readLines(res[j])
      if(length(tmp)==0) {
        suppressWarnings(structures[[j]]$structures<-fill)
        next
      } else {
        starts <- grep("NAME: ", tmp)
        
        if(length(starts) < 1) {
          msfinder.formula[[ramclustObj$cmpd[i]]]<- fill
          next
        }
        if(length(starts) > 1) {
          stops <- c(((starts[2:length(starts)])-1), (length(tmp)-1))
        } else {
          stops <- length(tmp)-1
        }
        out<-matrix(nrow = length(starts), ncol = length(tags))
        dimnames(out)[[2]]<-names(tags)
        for(k in 1:length(starts)) {
          d<-tmp[starts[k]:stops[k]]
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
          out[k,]<-vals
        }
        out<-data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
        if(any(out[,"name"] == "Spectral DB search")) {
          out <- out[-grep("Spectral DB search", out[,"name"]), ,drop = FALSE]
        }
        if(nrow(out) == 0) {
          suppressWarnings(structures[[j]]$structures<-fill)
        } else {
          suppressWarnings(structures[[j]]$structures<-out)
        }
        numfrag<-grep("Num Fragment", tmp)
        fragments<-as.list(rep(NA, length(numfrag)))
        if(length(fragments > 0)) {
          for(k in 1:length(numfrag)) {
            fragdat<-tmp[numfrag[k]:stops[k]]
            if(nchar(fragdat[length(fragdat)])==0) {
              fragdat = fragdat[1:(length(fragdat)-1)]
            }
            if(length(fragdat)==1) {
              fragments[[k]] <- fillfrag
              next
            }
            fragdat<-lapply(1:length(fragdat), FUN = function(x) {unlist(strsplit(fragdat[x], "\t")) })
            fragtab<-data.frame(matrix(nrow = length(fragdat)-1, ncol = length(fragdat[[2]])))
            names(fragtab)<- unlist(strsplit(unlist(strsplit(unlist(strsplit(fragdat[[1]], "\\(" ))[2], "\\)" ))[1], " "))
            for(l in 2:length(fragdat)) {
              fragtab[l-1,]<-fragdat[[l]]
            }
            fragments[[k]] <- data.frame(fragtab)
          }
          names(fragments)<-structures[[j]]$structures[,"id"]
        }
        
        suppressWarnings(structures[[j]]$fragments<-fragments)
        
      }
      
    }
    msfinder.structure[[do[i]]]<-structures
    setwd(mat.dir)
  }
  
  ramclustObj$msfinder.structure.details<-msfinder.structure
  
  setwd(home.dir)
  
  if(is.null(ramclustObj$history)) {
    ramclustObj$history <- ""
  }
  
  ramclustObj$history$msfinder <- paste(
    "MSFinder (Tsugawa 2016) was used for spectral matching,",
    "formula inference, and tentative structure assignment,",
    "and results were imported into the RAMClustR object.")
  
  return(ramclustObj)
  
}


