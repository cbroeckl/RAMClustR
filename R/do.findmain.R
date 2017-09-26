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
#' @param ads character: vector of allowed adducts, i.e. c("[M+H]+"). if NULL, default values of H+, Na+, K+, and NH4+, as monomer, dimer, and trimer, are assigned.
#' @param nls  character: vector of allowed neutral losses, i.e. c("[M+H-H2O]+").  if NULL, an extensive list derived from CAMERA's will be used. 
#' @param adwts numeric: vector of weights for adducts. length should match that of 'ads', else first value will be repeated.  if NULL (default), value of 1 is assigned to all. Only used for ramclustR scoring of M.
#' @param nlwts numeric: vector of weights for neutral losses. length should match that of 'nls', else first value will be repeated.  if NULL (default), value of 0.1 is assigned to all.  Only used for ramclustR scoring of M.
#' @param plot.findmain logical: should pdf polts be generated for evaluation?
#' @param writeMat logical: should indidual .mat files (for MSFinder) be generated in a 'mat' subdirectory in the 'spectra' folder?

#' @return a partially annotated ramclustR object.  base structure is that of a standard R heirarchical clustering output, with additional slots described in ramclustR documentation (?ramclustR).  New slots added after using the interpretMSSpectrum functionality include:
#' @return more slots to describe. 
#' @return M The inferred molecular weight of the compound giving rise to the each spectrum
#' @return M.ppm The ppm error of all the MS signals annotated, high error values should be considered 'red flags'
#' @return M.ann The annotated spectrum supporting the intepretation of M
#' @return use.findmain  Logical vector indicating whether findmain scoring (TRUE) or ramclustR scoring (FALSE) was used to support inference of M.  By default, findmain scoring is used.  When ramclustR scoring differs from findmain scoring, the scoring metric which predicts higher M is selected. 
#' @return M.ramclustr  M selected using ramclustR scoring
#' @return M.ppm.ramclustr ppm error of M selected using ramclustR scoring
#' @return M.ann.ramclustr annotated spectrum supporing M using ramclustR scoring
#' @return M.findmain M selected using findmain scoring
#' @return M.ppm.findmain ppm error of M selected using findmain scoring
#' @return M.ann.findmain annotated spectrum supporing M using findmain scoring
#' 
#' @references Jaeger C, Méret M, Schmitt CA, Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017 Aug 15;31(15):1261-1266. doi: 10.1002/rcm.7905. PubMed PMID: 28499062.
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' 
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
  
  #############
  ## if adducts and neutral losses and weights are undefined, set to defaults
  
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
  
  #############
  ##  create structure in ramclustObj for holding results
  
  M.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ppm.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ann.findmain <- as.list(rep(NA, max(ramclustObj$featclus)))
  M.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.ppm.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.rank.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.ann.ramclustr <- as.list(rep(NA, max(ramclustObj$featclus)))
  
  #############
  ## define (if null) full compound range for analysis
  
  if(is.null(cmpd)) {cmpd <- (1:max(ramclustObj$featclus)) }
  
  #############
  ## do findMain for all compounds in 'cmpd', append results to ramclustObj
  
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
    
    
    #############
    ## calculate ramclustR scoring
    
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

    best<-which.max(summaryscores)
    
    M.ramclustr[cl]<-summarytable[best, "neutral_mass"]
    M.ppm.ramclustr[cl]<-summarytable[best, "medppm"]
    M.rank.ramclustr[cl]<-best
    M.ann.ramclustr[[cl]]<-out[[best]]
    
    #############
    ## progress reporting
    if(10*round(cl/10, digits=0) == cl) {
      cat(cl, "of", max(ramclustObj$featclus), '\n')
    }
  } 
  
  #############
  ## append all results to ramclustObj
  ramclustObj$M.ramclustr <- M.ramclustr
  ramclustObj$M.ppm.ramclustr <- M.ppm.ramclustr
  ramclustObj$M.rank.ramclustr <- M.rank.ramclustr
  ramclustObj$M.ann.ramclustr <- M.ann.ramclustr
  ramclustObj$M.findmain <- M.findmain
  ramclustObj$M.ppm.findmain <- M.ppm.findmain
  ramclustObj$M.ann.findmain <- M.ann.findmain
  ramclustObj$use.findmain <- rep(TRUE, length(M.ppm.findmain))
  
  #############
  ## determine inconsistencies between two scoring methods, select higher MW
  resolve <- which(abs(ramclustObj$M.ramclustr - ramclustObj$M.findmain)  > (2*mzabs.error))
  for(i in resolve) {
    if(ramclustObj$M.ramclustr[i] >  ramclustObj$M.findmain[i]) {
      ramclustObj$use.findmain[i] <-FALSE
    }
  }
  
  ramclustObj$M <- ramclustObj$M.findmain
  ramclustObj$M[!ramclustObj$use.findmain] <- ramclustObj$M.ramclustr[!ramclustObj$use.findmain]
  ramclustObj$M.ann <- ramclustObj$M.ann.findmain
  resolved <- which(!ramclustObj$use.findmain)
  for(i in resolved) {
    ramclustObj$M.ann[[i]] <- ramclustObj$M.ann.ramclustr[[i]]
  }
  
  #############
  ## (optionally) plot spectra to pdf
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
  
  #############
  ## (optionally) write spectra to .mat format in spectra/mat directory
  if(writeMat) {
    if(!dir.exists("spectra")) {dir.create("spectra")}
    list.dirs()
    dir.create('spectra/mat')
    for(cl in cmpd) {
      ms<-ramclustObj$M.ann[[cl]] # ms <- matrix(round(rnorm(30, mean = 700, sd = 200), digits=3), ncol = 3)
      prcr<-which(ms[,"adduct"] %in% ads) # prcr = sample (1:10, 1)
      prcr<-prcr[which.max(ms[prcr,"int"])]
      prcmz<-ms[prcr,"mz"]  # prcmz <- ms[prcr, 1]
      prctype<-ms[prcr,"adduct"] # prctype = "[M+H]+"
      
      out<-paste(
        "NAME: ", ramclustObj$cmpd[cl], '\n',
        "RETENTIONTIME: ", round(ramclustObj$clrt[cl], digits=2), '\n',
        "PRECURSORMZ: ", prcmz, '\n',
        "PRECURSORTYPE: ", prctype, '\n',
        "IONTYPE: ", mode, '\n', 
        "SPECTRUMTYPE: Centroid", '\n',
        if((!is.null(ramclustObj$msmsint))) {paste("COLLISIONENERGY: ", as.character(ramclustObj$ExpDes[[2]][which(row.names(ramclustObj$ExpDes[[2]]) == "CE2"),1]), '\n', sep="")},
        "MSTYPE: ", "MS1", '\n',
        "Num Peaks: ", nrow(ms), '\n', sep="")
      for(i in 1:nrow(ms)) {
        out<-paste(out, ms[i,1], " ", ms[i,2], '\n', sep="")
      }
      
      
      if(!is.null(ramclustObj$msmsint)) {
        do<-which(ramclustObj$featclus == cl)
        msms<-cbind(
          'mz' = ramclustObj$fmz[do],
          'int' = ramclustObj$msmsint[do])
        msms<-msms[which(msms[,"mz"] <= (prcmz + 3)),, drop = FALSE]
        msms<-msms[order(msms[,"int"], decreasing=TRUE),]
        if(nrow(msms)>0) {
          out <- paste (out,
            "MSTYPE:" , "MS2", '\n' ,
            "Num Peaks: ", nrow(msms), '\n', 
            sep="")
          for(i in 1:nrow(msms)) {
            out<-paste(out, msms[i,1], " ", msms[i,2], '\n', sep="")
          }
        }
      }
      write(out, file=paste0("spectra/mat/", ramclustObj$cmpd[cl], ".mat"))
      
      
      #   sendToMSF(x = ms, 
      #             precursormz = prcmz, 
      #             precursortype = prctype, 
      #             outfile=paste0(getwd(), "/spectra/mat/", ramclustObj$cmpd[cl], ".mat")# ,
      #            # MSFexe = "K:/software/MSFinder/MS-FINDER program ver. 2.20"
      #             )
    }
  }
  
  cat("finished", '\n')
  return(ramclustObj)
  
}
