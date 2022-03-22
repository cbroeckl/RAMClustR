#' assign.z
#'
#' infer charge state of features in ramclustR object. 
#' @details Annotation of ramclustR spectra. looks at isotope spacing for clustered features to infer charge state for each feature and a max charge state for each compound 
#' 
#' @param ramclustObj ramclustR object to annotate
#' @param chargestate integer vector. vector of integers of charge states to look for.  default = c(1:5)
#' @param mzError numeric. the error allowed in charge state m/z filtering.  absolute mass units
#' @param nEvents integer. the number of isotopes necessary to assign a charnge state > 1.  default = 2.
#' @param minPercentSignal numeric.  the ratio of isotope signal (all isotopes) divided by total spectrum signal * 100 much be greater than minPercentSignal to evaluate charge state. Value should be between 0 and 100.  
#' @param assume1 logical.  when TRUE, m/z values for which no isotopes are found are assumed to be at z = 1.
#' @return returns a ramclustR object.  new slots holding: 
#' @return zmax. vector with length equal to number of compounds.  max charge state detected for that compound
#' @return fm. vector of inferred 'm', m/z value * z value
#' @return fz. vector of inferred 'z' values based on analysis of isotopes in spectrum.  
#' @concept ramlclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept charge state
#' @concept feature
#' @concept xcms
#' @concept MSFinder
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @export 


assign.z<-function(
  ramclustObj=NULL,
  chargestate=c(1:5),
  mzError=0.02,
  nEvents=2,
  minPercentSignal = 10,
  assume1=TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  ramclustObj$zmax<-rep(if(assume1) {1} else {NA}, length(ramclustObj$clrt))
  ramclustObj$fm<-rep(if(assume1) {ramclustObj$fmz} else {NA}, length(ramclustObj$fmz))
  ramclustObj$fz<-rep(if(assume1) {1} else {NA}, length(ramclustObj$fmz))
  ppmdifs<- 1.007276/chargestate
  for(i in 1:max(ramclustObj$featclus)) {
    members<-which(ramclustObj$featclus==i)
    mzvals<-ramclustObj$fmz[members]
    mzint<-ramclustObj$msint[members]
    mzdifs<-abs(outer(mzvals, mzvals, "-"))
    deltamz<-as.list(rep(NA, 2))
    deltamz[[1]]<-ppmdifs - mzError
    deltamz[[2]]<-ppmdifs + mzError
    for(j in 1:length(ppmdifs)) {
      check<-which(abs(mzdifs - ppmdifs[j])<= mzError, arr.ind=TRUE)
      if(nrow(check)==0) {break}
      if(nrow(check)>0) {
        check<-check[which(check[,2]>check[,1]),, drop=FALSE]
      }
      if(nrow(check)==0) {break}
      if(nrow(check)>=nEvents) {
        check<-unique(c(check))
      }
      if((100*sum(mzint[check])/ sum(mzint)) < minPercentSignal) {break}
      ramclustObj$fz[members[check]] <- j
      ramclustObj$fm[members[check]] <- ramclustObj$fmz[members[check]]*j
    }
    ramclustObj$zmax[i] <- max(ramclustObj$fz[members])
  }
  ramclustObj$m<-ramclustObj$fmz*ramclustObj$fz
  
  ramclustObj$history <- paste(ramclustObj$history, 
                               " Charge state detection was performed using the assign.z function ", 
                               " using parameters: ",
                               " chargestate = ", chargestate,
                               ", mzError = ", mzError,
                               ", nEvents = ", nEvents,
                               ", minPercentSignal = ", minPercentSignal,
                               ", and assume1 = ", assume1, ".", sep = "")
    

  return(ramclustObj)
}

