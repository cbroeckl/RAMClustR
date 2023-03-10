#' mergeRCobjects
#'
#' merge two ramclustR objects
#' @details  Two ramclustR objects are merged with this function, mapping features between them.  The first (ramclustObj.1) object use used as the template - all data in it is retained.  ramclustObj.2 is mapped to ramclustObj.1 feature by feature - only mapped features are retained.  A new ramlcustObj is returned, with a new SpecAbund dataset with the same column number as the ramclustObj.1$SpecAbund set. 
#' 
#' @param ramclustObj.1 ramclustR object 1: this object will be the base for the new object.  That is all the features from ramclustObj.1 will be retained.
#' @param ramclustObj.2 ramclustR object 2: this object will mapped and appended to racmlustObj1.  That is only features which appear consistent with those from ramclustObj.1 will be retained.
#' @param mztol numeric: absolute mass tolerance around mz
#' @param rttol numeric: feature retention time tolerance.  Value set by this option will be used during the initial anchor mapping phase.  Two times the standard error of the rt loess correction will be used for the full mapping.
#' @param course.rt.adj numeric: default = NULL. optional approximate retention time shift between ramclustObj.1 and ramclustObj.2.  i.e if the retention time of ramclustObj.1 is on average 15 seconds longer than that of ramclustobj.2, enter '15'.  if 1 is less than 2, enter a negative number.  This is applied before mapping to enable a smaller 'rttol' value to be used. 
#' @param mzwt numeric: when mapping features, weighting value used for similarities between feature mass values (see rtwt, intwt)
#' @param rtwt numeric: when mapping features, weighting value used for similarities between feature retention time values (see mzwt, intwt)
#' @param intwt numeric: when mapping features, weighting value used for similarities between ranked signal intensity values (see rtwt, mzwt)
#' @return returns a ramclustR object.  All values from ramclustObj.1 are retained.  SpecAbund dataset from ramclustObj.1 is moved to RC$SpecAbund.1, where RC is the new ramclustObj.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @author Corey Broeckling
#' @export 


mergeRCobjects <- function(
    ramclustObj.1 = NULL,  
    ramclustObj.2 = NULL,  
    mztol = 0.02, 
    rttol = 30,
    course.rt.adj = NULL,
    mzwt = 2,
    rtwt = 1,
    intwt = 3
) {
  
  if(is.null(ramclustObj.1)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj.1 = RC1", '\n')
  }
  if(is.null(ramclustObj.2)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj.2 = RC2", '\n')
  }
  
  ## make new RC object and add slots for new data
  
  newRC <- ramclustObj.1
  
  ## set variables if need be
  
  if(is.null(mztol)) {
    if(is.na(newRC$ExpDes[[2]]["mzdifftof",])) {
      mztol <- 0.75
      cat(" -- setting mztol to 0.75 Da", '\n')
    } else {
      mztol <- as.numeric(as.character(newRC$ExpDes[[2]]["mzdifftof",]))
    }
  }
  
  
  
  if(is.null(rttol)) {
    rttol <- round(100*median(newRC$clrtsd), digits = 2)
    cat (" -- rttol set to ", rttol, '\n')
  }
  
  ## generate a selectivity score based on intensity and isolation (how close is the nearest similar mass in rt)
  
  int.score <- (rank(newRC$msint)) / max(rank(newRC$msint))
  isol.score <- sapply(1:length(newRC$frt), FUN = function(x) {
    do <- which(abs(newRC$fmz - newRC$fmz[x]) <= mztol)
    do <- do[!do %in% x]
    if(length(do) == 0) {
      out <- 1
    } else {
      out <- exp(- 1/min(abs(newRC$frt[x] - newRC$frt[do])) )
    }
    return(out)
  }
  )
  isol.score<-(rank(isol.score, ties.method = "random"))/ max((rank(isol.score, ties.method = "random")))
  sel.score<-int.score * isol.score
  # hist(sel.score)
  
  ## map the most selective features across the two ramclustObjects
  
  map<-rep(NA, length(newRC$frt))
  if(length(map) > 8000) {qu <- 4}
  if(length(map) > 2000 & length(map) <= 7999) {qu <- 3}
  if(length(map) < 1000) {qu <- 2}
  
  if(!is.null(course.rt.adj)) {
    ramclustObj.2$frt <- ramclustObj.2$frt + course.rt.adj
  }
  
  
  marks<-which(sel.score >= quantile(sel.score)[[qu]])
  for(i in marks) {
    keep <- which(
      abs(newRC$frt[i] - ramclustObj.2$frt) <= rttol & 
        abs(newRC$fmz[i] - ramclustObj.2$fmz) <= mztol
    )
    if(length(keep) == 0) {
      next 
    }
    if(length(keep) == 1 ) {
      map[i] <- keep 
    }
    if(length(keep) > 1 ) {
      rtsim<-abs(newRC$frt[i] - ramclustObj.2$frt[keep])/rttol
      mzsim<-abs(newRC$fmz[i] - ramclustObj.2$fmz[keep])/mztol
      sims<-sapply(1:length(keep), FUN = function(x) {
        weighted.mean(c(rtsim[x], mzsim[x]), weights = c(rtwt, mzwt))
      }
      )
      # cat ('i = ', i, '\n')
      # print(data.frame("matches"=keep, rtsim, mzsim, sims))
      
      keep<-keep[which(sims == min(sims))]
      if(length(keep) > 1) {
        keep<-keep[which.max(ramclustObj.2$msint[keep])]
      }
      # cat("kept: ", keep, '\n', '\n')
      map[i]<-keep
    }
    rm(keep)
  }
  keep<-which(!is.na(map))
  if(length(keep) < 500) {
    warning(cat('only',  length(keep) ,'matching features: use output with caution', '\n'))
  }
  par(mfrow = c(2,2))
  plot(newRC$frt[keep], (newRC$frt[keep] - ramclustObj.2$frt[map[keep]]), xlab = "rt", ylab = "correction", 
       col = "gray", pch = 19, main = "anchor map", lwd = 3)
  # points(smooth.spline((newRC$frt[keep] - ramclustObj.2$frt[map[keep]]) ~ newRC$frt[keep], spar = 1.5), type = "l", col = 2)
  # fit<-lm((newRC$frt[keep] - ramclustObj.2$frt[map[keep]]) ~ 
  #           newRC$frt[keep] + I(newRC$frt[keep]^2), 
  #         weights =  log10(newRC$msin[keep]))
  
  
  fit<- loess((newRC$frt[keep] - ramclustObj.2$frt[map[keep]]) ~ 
                newRC$frt[keep] + I(newRC$frt[keep]^2), 
              span = rttol/2, degree = 2, surface = "direct")
  pred<-fitted(fit)
  se<-predict(fit, se = TRUE)
  points(newRC$frt[keep], pred, type = "l", col = 2, lwd = 3)
  plot(log10(newRC$msint[keep]) ,  log10(ramclustObj.2$msint[map[keep]]), 
       xlab = "log10(intensity.1)", ylab = "log10(intensity.2", 
       main = paste("anchor int r^2 = ", 
                    round(cor(
                      log10(ramclustObj.2$msint[map[keep]]), 
                      log10(newRC$msint[keep])), digits = 2), sep = ""),
       col = "gray", pch = 19)
  abline(lm(log10(ramclustObj.2$msint[map[keep]]) ~ log10(newRC$msint[keep])), col = 2, lwd = 3)
  sp<-(range(ramclustObj.2$frt)[2] - range(ramclustObj.2$frt)[1])/100
  new.rt <- newRC$frt[keep] - pred
  pred.frt <- sapply(1:length(ramclustObj.2$frt), FUN = function(x) {
    do<-which(abs(new.rt - ramclustObj.2$frt[x]) <= sp)
    return(mean(pred[do])+ramclustObj.2$frt[x])
  })
  
  
  ## map all features across the two ramclustObjects using new pred.frt
  
  rttol<- fit$s  # set new rttol to 1.5 times the standard error of the fit from selective features
  
  map<-rep(NA, length(newRC$frt))
  marks<-1:length(newRC$frt)
  for(i in marks) {
    keep <- which(
      abs(newRC$frt[i] - pred.frt) <= rttol & 
        abs(newRC$fmz[i] - ramclustObj.2$fmz) <= mztol
    )
    if(length(keep) == 0) { next }
    if(length(keep) == 1 ) {
      map[i] <- keep 
    }
    if(length(keep >1 )) {
      rtsim <- rank(abs(newRC$frt[i] - pred.frt[keep]))
      mzsim <- rank(abs(newRC$fmz[i] - ramclustObj.2$fmz[keep]))
      intsim <-  rank(abs(rank(newRC$msint)[i] - rank(ramclustObj.2$msint)[keep])) ##rank(RC$msint)
      sims<-sapply(1:length(keep), FUN = function(x) {
        weighted.mean(c(rtsim[x], mzsim[x], intsim[x]), weights = c(rtwt, mzwt, intwt))
      }
      )
      keep<-keep[which(sims == min(sims))]
      if(length(keep) > 1) {
        keep<-keep[which.max(ramclustObj.2$msint[keep])]
      }
      map[i]<-keep
    }
  }
  keep<-which(!is.na(map))
  cat(" --  ", 100*round(length(keep)/length(map), digits = 3), "% of features mapped", '\n')
  
  
  # par(mfrow = c(1,2))
  plot(newRC$frt[keep], (newRC$frt[keep] - pred.frt[map[keep]]), xlab = "rt", ylab = "correction", 
       col = "gray", pch = 19, main = "all rt vs pred rt")
  # points(smooth.spline((newRC$frt[keep] - ramclustObj.2$frt[map[keep]]) ~ newRC$frt[keep], spar = 1.5), type = "l", col = 2)
  # fit<-lm((newRC$frt[keep] - ramclustObj.2$frt[map[keep]]) ~ 
  #           newRC$frt[keep] + I(newRC$frt[keep]^2), 
  #         weights =  log10(newRC$msin[keep]))
  
  
  fit<- loess((newRC$frt[keep] - pred.frt[map[keep]]) ~ 
                newRC$frt[keep] + I(newRC$frt[keep]^2), 
              span = rttol/2, degree = 2, surface = "direct")
  pred<-fitted(fit)
  se<-predict(fit, se = TRUE)
  points(newRC$frt[keep], pred, type = "l", col = 2, lwd = 3)
  plot(log10(newRC$msint[keep]) ,  log10(ramclustObj.2$msint[map[keep]]), 
       xlab = "log10(intensity.1)", ylab = "log10(intensity.2", 
       main = paste("all r^2 = ", 
                    round(cor(
                      log10(ramclustObj.2$msint[map[keep]]), 
                      log10(newRC$msint[keep])), digits = 2), sep = ""),
       col = "gray", pch = 19)
  abline(lm(log10(ramclustObj.2$msint[map[keep]]) ~ log10(newRC$msint[keep])), col = 2, lwd = 3)
  
  
  newRC$MSdata <- matrix(nrow = (nrow(ramclustObj.1$MSdata) + nrow(ramclustObj.2$MSdata)), ncol = ncol(ramclustObj.1$MSdata))
  newRC$MSMSdata <- matrix(nrow = (nrow(ramclustObj.1$MSdata) + nrow(ramclustObj.2$MSdata)), ncol = ncol(ramclustObj.1$MSdata))
  newRC$SpecAbund <- matrix(nrow = nrow(newRC$MSdata), ncol = ncol(ramclustObj.1$SpecAbund))
  
  newRC$MSdata[1:nrow(ramclustObj.1$MSdata),] <- ramclustObj.1$MSdata
  newRC$MSdata[(nrow(ramclustObj.1$MSdata)+1):nrow(newRC$MSdata),which(!is.na(map))] <- ramclustObj.2$MSdata[,map[which(!is.na(map))]]
  newRC$MSMSdata[1:nrow(ramclustObj.1$MSMSdata),] <- ramclustObj.1$MSMSdata
  newRC$MSMSdata[(nrow(ramclustObj.1$MSMSdata)+1):nrow(newRC$MSMSdata),which(!is.na(map))] <- ramclustObj.2$MSMSdata[,map[which(!is.na(map))]]
  
  
  newRC$SpecAbund.1<-ramclustObj.1$SpecAbund
  newRC$SpecAbund.2<-ramclustObj.2$SpecAbund
  cat(" --   collapsing spectra", '\n' )
  
  wts <- colSums(ramclustObj.1$MSdata, na.rm = TRUE)
  for (ro in 1:nrow(newRC$MSdata)) {
    for (co in 1:ncol(newRC$SpecAbund.1)) {
      newRC$SpecAbund[ro, co] <- weighted.mean(
        newRC$MSdata[ro, which(ramclustObj.1$featclus == co)], 
        wts[which(ramclustObj.1$featclus ==  co)],
        na.rm = TRUE)
    }
  }
  
  # wts <- colMeans(newRC$MSdata, na.rm = TRUE)
  # for (ro in 1:nrow(newRC$MSdata)) {
  #   for (co in 1:ncol(newRC$SpecAbund)) {
  #     if(length(which(map == co)) > 0) {
  #       newRC$SpecAbund[ro, co] <- weighted.mean(ramclustObj.2$MSdata[ro, which(map == co)], 
  #                                                                               wts[which(map ==  co)])
  #     } else {
  #       # cat(co, '\n')
  #       newRC$SpecAbund[ro, co] <- 0
  #     }
  #     
  #   }
  # }
  
  dimnames(newRC$SpecAbund)[[2]] <- dimnames(ramclustObj.1$SpecAbund)[[2]]
  rnames <- c(dimnames(ramclustObj.1$SpecAbund)[[1]], dimnames(ramclustObj.2$SpecAbund)[[1]])
  if(length(rnames) == nrow(newRC$SpecAbund) ) {
    dimnames(newRC$SpecAbund)[[1]] <- rnames
  } else {
    rnames <- c(dimnames(ramclustObj.1$MSdata)[[1]], dimnames(ramclustObj.2$MSdata)[[1]])
    dimnames(newRC$SpecAbund)[[1]] <- rnames
  }
  
  
  newRC$pred.frt2 <- pred.frt
  newRC$frt2 <- rep(NA, length(newRC$frt))
  newRC$frt2[!is.na(map)] <- ramclustObj.2$frt[map[!is.na(map)]]
  
  newRC$fmz2 <- rep(NA, length(newRC$frt))
  newRC$fmz2[!is.na(map)] <- ramclustObj.2$fmz[map[!is.na(map)]]
  
  cat(" --   finished", '\n' )
  
  return(newRC)
}

