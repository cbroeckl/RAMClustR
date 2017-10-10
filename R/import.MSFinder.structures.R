#' import.MSFinder.structures
#'
#' After running MSFinder on .mat or .msp files, import the structure and scores that were predicted
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir optional path to .mat directory
#' @param msp.dor optional path to .msp directory
#' @details this function imports the output from the MSFinder program to annotate the ramclustR object
#' @return an annotated ramclustR object
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @return ramclustR object with new $ slots for: 
#' @return - msfinder.structure: data frame of one row describing the best structure match for the selected best formula match
#' @return - msfinder.structure$fragments:  data frame the fragment ion interpretations supporting the best structure match
#' @return - msfinder.structure$details: list containing the best structure match for all formulas.  Each list element represents one formula (with the formula as the name), and is composed of a list:
#' @return    - structures: table of most likely structures
#' @return    - fragments: list of data frames, with one element for each from from 'structures'.  names are assigned based on 'id' column from structures.   
#' @author Corey Broeckling
#' @export

import.MSFinder.structures <- function (
  ramclustObj = RC,
  mat.dir = NULL,
  msp.dir = NULL,
  MSFinder.dir = "C:\\MSFinder\\MS-FINDER_2.20"
) {
  
  if(is.null(ramclustObj$msfinder.formula)) {
    warning("trying to run 'import.MSFinder.formulas' first")
    ramclustObj<-import.MSFinder.formulas(ramclustObj = RC, mat.dir = NULL, msp.dir = NULL)
  }
  
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
  
  
  do <- list.dirs(mat.dir, full.names = FALSE, recursive = FALSE)
  
  do <- do[do %in% ramclustObj$cmpd]
  
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
  data.frame(fill, check.names = FALSE, stringsAsFactors = FALSE)
  
  msfinder.structure<-as.list(rep("", length(ramclustObj$cmpd)))
  names(msfinder.structure)<-ramclustObj$cmpd
  
  for(i in 1:length(do)) {
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
          msfinder.formula[[cmpd[i]]]<- fill
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
        for(k in 1:length(numfrag)) {
          fragdat<-tmp[numfrag[k]:stops[k]]
          if(nchar(fragdat[length(fragdat)])==0) {
            fragdat = fragdat[1:(length(fragdat)-1)]
          }
          if(length(fragdat)==1) {
            fragments[[k]] <- NULL
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
        suppressWarnings(structures[[j]]$fragments<-fragments)
        
      }
      
    }
    msfinder.structure[[do[i]]]<-structures
    setwd(mat.dir)
  }
  
  ## need to fix below!!!
  
  ramclustObj$msfinder.structure.details<-msfinder.structure
  
  ramclustObj$msfinder.structure<-lapply(1:length(msfinder.structure), FUN = function(x) {
    # for (x in 1:length( msfinder.structure)) {
    if(!is.na(ramclustObj$msfinder.formula[x])) {
      if(nrow(ramclustObj$msfinder.structure.details[[x]][[ramclustObj$msfinder.formula[x]]]$structures)>0) {
        # if(nrow(ramclustObj$msfinder.structure[[x]][[ramclustObj$msfinder.formula[x]]]$structures)>0) {
        best<-which.max(ramclustObj$msfinder.structure.details[[x]][[ramclustObj$msfinder.formula[x]]]$structures[,"totalscore"])
        ramclustObj$msfinder.structure.details[[x]][[ramclustObj$msfinder.formula[x]]]$structures[best,]
      } else {NA}
    } else {NA}
  }
  )
  
  ramclustObj$msfinder.structure.fragments<-lapply(1:length(msfinder.structure), FUN = function(x) {
    # for (x in 1:length( msfinder.structure)) {
    if(!is.na(ramclustObj$msfinder.structure[[x]][1])) {
      ramclustObj$msfinder.structure.details[[x]][[ramclustObj$msfinder.formula[x]]][["fragments"]][[ramclustObj$msfinder.structure[[x]][1,"id"]]]
    } else {
      NA
    }
  }
  )
  
  return(ramclustObj)
  
}
