#' import.MSFinder.formulas
#'
#' After running MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir optional path to .mat directory
#' @param msp.dor optional path to .msp directory
#' @details this function imports the output from the MSFinder program to annotate the ramclustR object
#' @return an annotated ramclustR object
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @author Corey Broeckling
#' #@export

import.MSFinder.search <- function (
  ramclustObj = RC,
  mat.dir = NULL,
  msp.dir = NULL
) {

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
  
  dirs<-list.dirs(mat.dir, recursive = FALSE, full.names = FALSE) 

  do<-list.files(mat.dir, pattern = "Spectral DB search.sfd", full.names = TRUE, recursive = TRUE)

    
  if(!identical(sort(ramclustObj$cmpd), sort(dirs))) {
    stop("spectral match directory names do not match compound names: please ensure you have the correct pair of ramclust object and directory")
  } else {
    cmpd<-dirs
  }
  # 
  # tags<-c(
  #   "NAME: ",
  #   "INCHIKEY: ",
  #   "SMILES: ",
  #   RESOURCES: 
  #     SubstructureInChIKeys: 
  #     RETENTIONTIME: -1
  #   RETENTIONINDEX: 0
  #   TotalBondEnergy: -1
  #   TotalScore: 3.9117
  #   TotalHrRulesScore: 0
  #   TotalBondCleavageScore: 0
  #   TotalMassAccuracyScore: 0
  #   TotalFragmentLinkageScore: 0
  #   TotalBondDissociationEnergyScore: 0
  #   DatabaseScore: 0
  #   SubstructureAssignmentScore: 0
  #   RtSimilarityScore: 0
  #   RiSimilarityScore: 0
  #   Num Peaks: 289
  # )
  # names(tags)<-tolower(gsub(": ", "", tags))
  # 
  # msfinder.formula<-as.list(rep(NA, length(ramclustObj$cmpd)))
  # names(msfinder.formula)<-ramclustObj$cmpd
  # 
  # tmp<-readLines(do[[i]])
  # 
  # for(i in 1:length(do)) {
  #   tmp<-readLines(do[[i]])
  #   starts <- grep("NAME: ", tmp)
  #   stops <- c(((starts[2:length(starts)])-1), (length(tmp)-1))
  #   out<-matrix(nrow = length(starts), ncol = length(tags))
  #   dimnames(out)[[2]]<-names(tags)
  #   for(j in 1:length(starts)) {
  #     d<-tmp[starts[j]:stops[j]]
  #     vals<-sapply(1:length(tags), FUN = function(x) {
  #       m<-grep(tags[x], d, fixed = TRUE)
  #       if(length(m)==0) {
  #         NA
  #       } else {
  #         if(length(m)>1) {
  #           m<-m[1]
  #         }
  #         gsub(tags[x], "", d[m])
  #       }
  #     }
  #     )
  #     out[j,]<-vals
  #   }
  #   out<-data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  #   if(any(out$name == "Spectral DB search")) {
  #     out <- out[-grep("Spectral DB search", out$name),]
  #   }
  #   msfinder.formula[[cmpd[i]]]<-out
  # }
  # 
  # ramclustObj$msfinder.formula<-sapply(1:length(msfinder.formula), FUN = function(x) {
  #   if(nrow(msfinder.formula[[x]])>0) {
  #     msfinder.formula[[x]][1,"name"]
  #   } else {
  #     NA
  #   }
  # }
  # )
  # 
  # ramclustObj$msfinder.formula.score<-as.numeric(sapply(1:length(msfinder.formula), FUN = function(x) {
  #   if(nrow(msfinder.formula[[x]])>0) {
  #     msfinder.formula[[x]][1,"totalscore"]
  #   } else {
  #     NA
  #   }
  # }
  # ))
  # 
  # ramclustObj$msfinder.formula.details<-msfinder.formula
  # 
  # return(ramclustObj)
  # 
}
