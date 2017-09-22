#' do.findmain
#'
#' cluster annotation function 
#'
#' This is the Details section
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param cmpd integer: vector defining compound numbers to annotated.  if NULL (default), all compounds
#' @param mode character: "positive" or "negative"
#' @param mzabs.error numeric: absolute mass deviation allowd
#' @param ppm.error numeric: ppm mass error added to mzabs.error
#' @param ads character: vector of allowed adducts, i.e. c("[M+H]+")
#' @param nls  character: vector of allowed neutral losses, i.e. c("[M+H-H2O]+")
#' @param adwts numeric: vector of weights for adducts. length should match that of 'ads', else first value will be repeated.  if NULL (default), value of 1 is assigned to all. Only used for ramclustR scoring of M.
#' @param nlwts numeric: vector of weights for neutral losses. length should match that of 'nls', else first value will be repeated.  if NULL (default), value of 0.1 is assigned to all.  Only used for ramclustR scoring of M.
#' @param plot.findmain logical: should pdf polts be generated for evaluation?
#' @param writeMat logical: should indidual .mat files (for MSFinder) be generated in a 'mat' subdirectory in the 'spectra' folder?

#' @return an annotated ramclustR object.  base structure is that of a standard R heirarchical clustering output
#' @author Corey Broeckling
#' @export


do.findmain<-function(ramclustObj = RC, 
                      cmpd = NULL,
                      mode = "positive", 
                      mzabs.error = 0.01, 
                      ppm.error = 10, 
                      ads = NULL, 
                      nls = NULL, 
                      adwts = NULL, 
                      nlwts=NULL, 
                      plot.findmain=TRUE, 
                      writeMat=TRUE) {
  require(InterpretMSSpectrum)
  if(is.null(ads)) {
    ads<-c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[2M+H]+", "[2M+Na]+", 
           "[2M+K]+", "[2M+NH4]+", "[3M+H]+", "[3M+Na]+", "[3M+K]+", "[3M+NH4]+")
  }
  
  if(is.null(nls)) {
    nls<-c("[M+H-NH3]+", "[M+H-H2O]+", "[M+H-COCH2]+", "[M+H-CO2]+", "[M+H-NH3-CO2]+", 
           "[M+H-NH3-HCOOH]+", "[M+H-NH3-H2O]+", "[M+H-NH3-COCH2]+", "[M+H-S]+", 
           "[M+H-S-NH3-HCOOH]+", "[M+H-H4O2]+", "[M+H-CH2]+", "[M+H-O]+", 
           "[M+H-C2H2]+", "[M+H-C2H4]+", "[M+H-CO]+", "[M+H-C3H6]+", "[M+H-C2H4O]+", 
           "[M+H-C4H6]+", "[M+H-C3H4O]+", "[M+H-C4H8]+", "[M+H-C5H8]+", 
           "[M+H-C4H6O]+", "[M+H-C5H10]+", "[M+H-C6H12]+", "[M+H-C4H8O2]+", 
           "[M+H-H2O-HCOOH]+", "[M+H-CH4]+", "[M+H-CH2O]+", "[M+H-C2H6]+", 
           "[M+H-CH3OH]+", "[M+H-C3H4]+", "[M+H-C3H6O]+", "[M+H-CO2-C3H6]+", 
           "[M+H-SO3]+", "[M+H-SO3-H2O]+", "[M+H-SO3-H2O-NH3]+", "[M+H-NH3-C3H4]+", 
           "[M+H-H2O-CO2]+", "[M+H-H2O-H2O-C2H4O]+", "[M+H-NH3-CO-CO]+", 
           "[M+H-NH3-CO-COCH2]+", "[M+H-C8H6O]+", "[M+H-C8H6O-NH3]+", "[M+H-C8H6O-H2O]+", 
           "[M+H-C2H2O2]+", "[M+H-C2H4O2]+", "[M+H-C5H8O]+", "[M+H-NH3-CO2-CH2O]+", 
           "[M+H-NH3-CO2-NH3-H2O]+", "[M+H-NH3-CO2-C3H4O]+", "[M+H-NH3-CO2-C5H8]+", 
           "[M+H-HCOOH-HCOOH]+", "[M+H-C2H4-CO2]+", "[M+H-C2H4-HCOOH]+", 
           "[M+H-NH3-H2O-H2O]+", "[M+H-H2O-C2H2O2]+", "[M+H-COCH2-C4H8]+", 
           "[M+H-NH3-NH3-C3H4]+", "[M+H-C2H4O2-CH3OH]+", "[M+H-C3H6O-CH3OH]+", 
           "[M+H-NH3-CO-COCH2-C4H6O]+", "[M+H-C4H6-H2O]+", "[M+H-C4H6-C2H4]+", 
           "[M+H-C4H6-NH3-H2O]+", "[M+H-C4H6-COCH2]+", "[M+H-C4H6-C4H6O]+", 
           "[M+H-C3H4O-C4H6]+", "[M+H-C3H4O-C4H8O2]+", "[M+H-C4H8-C4H6]+", 
           "[M+H-NH3-HCOOH-CH3OH]+", "[M+H-NH3-C2H6]+", "[M+H-NH3-C8H6O-CH2]+", 
           "[M+H-NH3-C3H4-COCH2]+", "[M+H-C3H9N]+", "[M+H-C3H9N-C2H4O2]+", 
           "[M+H-C6H10O7]+", "[M+H-C6H10O7+H2O]+", "[M+H-C6H10O7-(H2O)2]+", 
           "[M+H-C6H12O6]+", "[M+H-C6H12O6+H2O]+", "[M+H-C6H12O6-H2O]+", 
           "[M+H-C5H10O5]+", "[M+H- C5H10O5+H2O]+", "[M+H- C5H10O5-H2O]+", 
           "[M+H-C2H8NO4P]+", "[M+H-(H2O)3-CO]+", "[M+H-C6H13NO2]+", "[M+H-C5H11NO2]+", 
           "[M+H-CH3S]+", "[M+H-C8H8O2]+", "[M+H-C12H22O11]+", "[M+H-C12H24O12]+", 
           "[M+H-C12H20O10]+")
  }
  
  if(length(adwts)>0) {
    if(length(adwts) != length(ads)) {
      if(length(adwts) > 1) {
        warning("adduct weight length not equal to adduct length: assigning weight ", adwts[1], " to all adduct weights")
      }
      adwts<-rep(adwts[1], length(ads)) 
    } 
  } else  {
    adwts = rep(1, length(ads))
  }
  
  
  if(length(nlwts)>0) {
    if(length(nlwts) != length(nls)) {
      if(length(nlwts) > 1) {
        warning("adduct weight length not equal to adduct length: assigning weight ", nlwts[1], " to all adduct weights")
      }
      nlwts<-rep(nlwts[1], length(nls)) 
    } 
  } else  {
    nlwts = rep(0.1, length(nls))
  }
  
  adnlwts<-c(adwts, nlwts); names(adnlwts)<-c(ads, nls)
  
  # cat(names(adnlwts),'\n')
  # cat(adnlwts,'\n')

  M.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ppm.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ann.findmain <- as.list(rep(NA, max(ramclustObj$featclus)))
  M.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.ppm.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.rank.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.ann.ramclustr <- as.list(rep(NA, max(ramclustObj$featclus)))
  
  
  if(is.null(cmpd)) {cmpd <- (1:max(ramclustObj$featclus)) }
  
  for(cl in cmpd){
    s<-data.frame("mz"=ramclustObj$fmz[which(ramclustObj$featclus==cl)], "int"=ramclustObj$msint[which(ramclustObj$featclus==cl)])
    
    out<-findMAIN(
      s, 
      adductmz = NULL, 
      ionmode = mode,
      adducthyp = ads, 
      ms2spec = NULL, 
      rules = c(ads, nls),
      mzabs = mzabs.error,
      ppm = ppm.error, 
      mainpkthr = 0.1, 
      collapseResults = FALSE)
    
    summarytable<-summary(out)
    
    M.findmain[cl]<-summarytable[1, "neutral_mass"]
    M.ppm.findmain[cl]<-summarytable[1, "medppm"]
    M.ann.findmain[[cl]]<-out[[1]]
    
    ## set all major adduct ppm errors to half ppm.error, or 5
    for(y in 1:length(out)) {
      keep<-which(!is.na(out[[y]][,"adduct"]))  ## which are annotated peaks
      out[[y]][keep[which(is.na(out[[y]][keep,"ppm"]))],"ppm"] <- ppm.error/2
    }
    
    summaryscores<-sapply(
      1:length(out), 
      FUN=function(x) 
      {
        keep<-which(!is.na(out[[x]][,"adduct"]))  ## which are annotated peaks
        wt<- adnlwts[out[[x]][keep,"adduct"]]
        int<-(out[[x]][keep, "int"])^0.1   ## use square root of intensity to prevent bias toward base peak
        mzerr<-out[[x]][keep, "ppm"] 
        # mzerr[is.na(mzerr)]<-0
        mzerr<-round(exp(-mzerr^2/(2*(ppm.error^2)) ), digits = 4)  ## this is a sigmoid function to weight high ppm error peaks lower
        #massorder<-order(out[[x]][keep, "mz"])
        massorder<-sqrt(order(out[[x]][keep, "mz"]))
        massorder<-(massorder)/max((massorder))
        return(sum((massorder *int * mzerr * wt), na.rm=TRUE) + (0.1*length(keep)))  ## product of the intensity and the scaled mzerror
        # data.frame(keep, wt, int, mzerr, 'ppm'=out[[x]][keep, "ppm"], massorder, 'score'=(massorder *int * mzerr * wt))
      }
    )
    
    ## the result with the best summaryscore can be selected 
    ## and tagged as such somehow
    
    best<-which.max(summaryscores)
    
    M.ramclustr[cl]<-summarytable[best, "neutral_mass"]
    M.ppm.ramclustr[cl]<-summarytable[best, "medppm"]
    M.rank.ramclustr[cl]<-best
    M.ann.ramclustr[[cl]]<-out[[best]]
    
    if(10*round(cl/10, digits=0) == cl) {
      cat(cl, "of", max(ramclustObj$featclus), '\n')
    }
  } 
  
  ramclustObj$M.ramclustr <- M.ramclustr
  ramclustObj$M.ppm.ramclustr <- M.ppm.ramclustr
  ramclustObj$M.rank.ramclustr <- M.rank.ramclustr
  ramclustObj$M.ann.ramclustr <- M.ann.ramclustr
  ramclustObj$M.findmain <- M.findmain
  ramclustObj$M.ppm.findmain <- M.ppm.findmain
  ramclustObj$M.ann.findmain <- M.ann.findmain
  ramclustObj$use.findmain <- rep(TRUE, length(M.ppm.findmain))
  
  resolve <- which(abs(ramclustObj$M.ramclustr - ramclustObj$M.findmain)  > (2*mzabs.error))
  for(i in resolve) {
    if(ramclustObj$M.ramclustr[i] >  ramclustObj$M.findmain[i]) {
      ramclustObj$use.findmain[i] <-FALSE
    }
  }
  
  ramclustObj$M<-ramclustObj$M.findmain
  ramclustObj$M[!ramclustObj$use.findmain]<-ramclustObj$M.ramclustr[!ramclustObj$use.findmain]
  ramclustObj$M.ann <- ramclustObj$M.ann.findmain
  resolved <- which(!ramclustObj$use.findmain)
  for(i in resolved) {
    ramclustObj$M.ann[[i]]<-ramclustObj$M.ann.ramclustr[[i]]
  }
    
  if(plot.findmain) {
    cat("plotting findmain annotation results", '\n')
    pdf("spectra/findmainPlots.pdf", width=10, height = 4.6)
    par(mfrow=c(1,2))
    for(cl in cmpd) {
      PlotSpec(x=ramclustObj$M.ann.ramclustr[[cl]], txt=ramclustObj$M.ann.ramclustr[[cl]][,c("mz","adduct")])
      title(main=list(paste(
        cl, ":",
        "M.ramclustr =", 
        round(ramclustObj$M.ramclustr[cl], digits=4), "( +/-", round(ramclustObj$M.ppm.ramclustr[cl], digits=1),
        "ppm )"), font = if(ramclustObj$use.findmain[cl]) {1} else {2}, col = if(ramclustObj$use.findmain[cl]) {1} else {2} ))
      PlotSpec(x=ramclustObj$M.ann.findmain[[cl]], txt=ramclustObj$M.ann.findmain[[cl]][,c("mz","adduct")])
      title(main = list(paste(
        cl, ":",
        "M.findmain =", 
        round(ramclustObj$M.findmain[cl], digits=4), "( +/-", round(ramclustObj$M.ppm.findmain[cl], digits=1), 
        "ppm )"  ),
            font = if(ramclustObj$use.findmain[cl]) {2} else {1}, col = if(ramclustObj$use.findmain[cl]) {2} else {1}  ))
    }
    dev.off()
  }
  
  if(writeMat) {
    dir.create('spectra/mat')
    for(cl in cmpd) {
      ms<-ramclustObj$M.ann[[cl]]
      prcr<-which(ms[,"adduct"] %in% ads)
      prcr<-prcr[which.max(ms[prcr,"int"])]
      prcmz<-ms[prcr,"mz"]
      prctype<-ms[prcr,"adduct"]
      
      if(!is.null(ramclustObj$msmsint)) {
        do<-which(ramclustObj$featclus == cl)
        msms<-cbind(
          'mz' = ramclustObj$fmz[do],
          'int' = ramclustObj$msmsint[do])
        msms<-msms[which(msms[,"mz"] <= (prcmz + 3)),]
        msms<-msms[order(msms[,"int"], decreasing=TRUE)]
      }
      
      sendToMSF(x = ms, 
                precursormz = prcmz, 
                precursortype = prctype, 
                outfile=paste0(getwd(), "/spectra/mat/", ramclustObj$cmpd[cl], ".mat")# ,
               # MSFexe = "K:/software/MSFinder/MS-FINDER program ver. 2.20"
                )
    }
  }
  cat("finished", '\n')
  return(ramclustObj)
  
}

# rcbackup<-RC
# RC2<-do.findmain(ramclustObj = RC, cmpd = c(1:15, 625) )
