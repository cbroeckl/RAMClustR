#' import.MSFinder.formulas
#'
#' After running MSFinder on .mat or .msp files, import the formulas and scores that were predicted
#' 

import.MSFinder.formulas <- function (
  ramclustObj = RC,
  mat.dir = NULL,
  rt.form.filter = TRUE
) {
  
  if(rt.form.filter) {
    library(CHNOSZ)
  }
    
  if(is.null(mat.dir)) {
    mat.dir = paste0(getwd(), "/spectra/mat")
  }
  if(!dir.exists(mat.dir)) {
    stop(paste("there is no directory called:", '\n', getwd()))
  }

  do<-list.files(mat.dir, pattern = ".fgt", full.names = TRUE)
  
  cmpd<-gsub(".fgt", "", basename(do))
  
  tags<-c(
    "NAME: ",
    "EXACTMASS: ",
    "ISSELECTED: ",
    "MASSDIFFERENCE: ",
    "TOTALSCORE: ",
    "ISOTOPICINTENSITY[M+1]: ",
    "ISOTOPICINTENSITY[M+2]: ",
    "ISOTOPICDIFF[M+1]: ",
    "ISOTOPICDIFF[M+2]: ",                                                                                                                                 
    "MASSDIFFSCORE: ",                                                                                                                                      
    "ISOTOPICSCORE: ",                                                                                                                                      
    "PRODUCTIONSCORE: ",                                                                                                                                    
    "NEUTRALLOSSSCORE: ",                                                                                                                                   
    "PRODUCTIONPEAKNUMBER: ",                                                                                                                               
    "PRODUCTIONHITSNUMBER: ",                                                                                                                               
    "NEUTRALLOSSPEAKNUMBER: ",                                                                                                                              
    "NEUTRALLOSSHITSNUMBER: ",                                                                                                                              
    "RESOURCENAMES: ",                                                                                                                                       
    "RESOURCERECORDS: ",                                                                                                                                    
    "ChemOntDescriptions: ",                                                                                                                                 
    "ChemOntIDs: ",                                                                                                                                          
    "ChemOntScores: ",                                                                                                                                       
    "ChemOntInChIKeys: ",                                                                                                                                    
    "PUBCHEMCIDS: "
  )
  names(tags)<-tolower(gsub(": ", "", tags))
  
  msfinder.formula<-as.list(rep(NA, length(ramclustObj$cmpd)))
  names(msfinder.formula)<-ramclustObj$cmpd
  
  
  
  for(i in 1:length(do)) {
    tmp<-readLines(do[[i]])
    starts <- grep("NAME: ", tmp)
    stops <- c(((starts[2:length(starts)])-1), (length(tmp)-1))
    out<-matrix(nrow = length(starts), ncol = length(tags))
    dimnames(out)[[2]]<-names(tags)
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
    }
    out<-data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
    if(any(out$name == "Spectral DB search")) {
      out <- out[-grep("Spectral DB search", out$name),]
    }
    msfinder.formula[[cmpd[i]]]<-out
  }
  
  if(rt.form.filter) {
    df<-data.frame(
      'cmpd' = rep(ramclustObj$cmpd[1], nrow(msfinder.formula[[1]])),
      'rt' = rep(ramclustObj$clrt[1], nrow(msfinder.formula[[1]])), 
      msfinder.formula[[1]], stringsAsFactors = FALSE, check.names = FALSE)
    for(i in 2:length(msfinder.formula)) {
      tmp<-data.frame(
        'cmpd' = rep(ramclustObj$cmpd[i], nrow(msfinder.formula[[i]])),
        'rt' = rep(ramclustObj$clrt[i], nrow(msfinder.formula[[i]])), 
        msfinder.formula[[i]], stringsAsFactors = FALSE, check.names = FALSE)
      df<-rbind(df, tmp, deparse.level = 1)
    }
    df<-as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    makes<-lapply(df$name, FUN = makeup)
    elemnames<-unique(names(unlist(makes)))
    elems<-matrix(ncol = length(elemnames), nrow = length(makes))
    elems[,]<-0
    dimnames(elems)[[2]]<-elemnames
    for(i in 1:length(makes)) {
      elems[i,names(makes[[i]])]<-makes[[i]]
    }
    dbe<-elems[,"C"] - (elems[,"H"]/2) + (elems[,"N"]/2) + 1
    O<-elems[,"O"] 
    rt<-df$rt
    wt<-as.numeric(df$totalscore)
    plot(rt, O, pch = 19, cex = 0.1*wt)
    fit<-lm(rt~ dbe + elems +I(dbe^2) + I(elems^2) , weights = wt^10); summary(fit)
    anova(fit)
    summary(fit)
    predrt<-predict(fit)
    plot(rt, predrt, pch = 19, cex= 0.4, col = 4)
    abline(smooth.spline(predrt~rt))
  }

  ramclustObj$msfinder.formula<-msfinder.formula

  
}
