#' rc.feature.replace.na.impute.knn
#'
#' replaces any NA (and optionally zero) values with small signal (20% of minimum feature signal value + 20% random noise)
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param replace.int default = 0.2.  proportion of minimum feature value to replace NA (or zero) values with
#' @param replace.noise default = 0.2.  proportion ofreplace.int value by which noise is added via 'jitter'
#' @param replace.zero logical if TRUE, any zero values are replaced with noise as if they were NA values 
#' @details noise is added by finding for each feature the mimimum detected value, multiplying that value by replace.int, then adding (replace.int*replace.noise) noise.  abs() is used to ensure no negative values result.  
#' @return  ramclustR object with NA and zero values removed.     
#'  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
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

rc.feature.replace.na.impute.knn  <- function(
  ramclustObj=NULL,
  replace.zero = TRUE,
  samp.max.missing = 0.8,
  feat.max.missing = 0.5,
  remove.samples = TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("please provide a ramclustR Object as input.", '\n')
  }
  
  if(!is.logical(replace.zero)) {
    stop("replace.zero must be logical",'\n')
  }
  
  ## Step 1. optionally replace zeros with NA so they can be treated identically.
  if(replace.zero) {
    ramclustObj$MSdata[which(ramclustObj$MSdata == 0, arr.ind = TRUE)] <- NA
    if(!is.null(ramclustObj$MSMSdata)) {
      ramclustObj$MSMSdata[which(ramclustObj$MSMSdata == 0, arr.ind = TRUE)] <- NA
    }
  }
  
  
  d <- ramclustObj$MSdata
  if(!is.null(ramclustObj$MSMSdata)) {
    d <- rbind(d, (ramclustObj$MSMSdata))
  }
  
  ## correlate predictor testing
  # for(i in 1:ncol(d)) {
  #   rval <- as.vector(cor(d[,i], d[,1:ncol(d)], use = "pairwise.complete.obs")^2)
  #   rval[which(is.na(rval))] <- 0
  #   use.feat <- which(order(rval, decreasing = TRUE) <=10)
  # }
  
  ## Step 2. optionally remove samples with too many missing values
  if(remove.samples) {
    missing <- sapply(1:nrow(ramclustObj$MSdata), FUN = function(x) {
      length(which(is.na(ramclustObj$MSdata[x,])))/ncol(ramclustObj$MSdata)
    })
    remove <- which(missing >= samp.max.missing)
    if(length(remove)> 0) {
      cat("removed samples:", '\n')
      cat(paste("  ", as.character(ramclustObj$phenoData$sample.names)[remove], collapse = '\n'))
      d <- d[,-remove]
      ramclustObj$phenoData <- ramclustObj$phenoData[-remove,]
    }
    ramclustObj$MSdata <- ramclustObj$MSdata[-remove,]
    if(!is.null(ramclustObj$MSMSdata)) {ramclustObj$MSMSdata <- ramclustObj$MSMSdata[-remove,]}
    ramclustObj$phenoData <- ramclustObj$phenoData[-remove,]
  }
  
  
  ## step 3.  for MS data and then (optionally) MSMSdata, 
  ## fill any features that have > samp.max.missing values
  ## with low level noise. 
  ## low level noise is defined as the minimum of the area 
  ##that feature value in either MS or
  ## MSMSdata * 0.5 * rnorm(n, mean = 1, sd = 0.05)
  
  feat.min <- suppressWarnings(apply(ramclustObj$MSdata, 2, FUN = "min", na.rm = TRUE))
  if(!is.null(ramclustObj$MSMSdata)) {
    feat.min2 <- suppressWarnings(apply(ramclustObj$MSMSdata, 2, FUN = "min", na.rm = TRUE))
    feat.min <- pmin(feat.min, feat.min2, na.rm = TRUE)
  }
  
  ## MSdata
  missing <- sapply(1:ncol(ramclustObj$MSdata), FUN = function(x) {
    length(which(is.na(ramclustObj$MSdata[,x])))/nrow(ramclustObj$MSdata)
  })
  fill <- which(missing > feat.max.missing)
  for(i in fill) {
    fill.current <- which(is.na(ramclustObj$MSdata[,i]))
    ramclustObj$MSdata[fill.current,i] <- {
      feat.min[i] * 0.5 * rnorm(length(fill.current), mean = 1, sd = 0.05)
    }
  }
  
  if(!is.null(ramclustObj$MSMSdata)) {
    missing <- sapply(1:ncol(ramclustObj$MSMSdata), FUN = function(x) {
      length(which(is.na(ramclustObj$MSMSdata[,x])))/nrow(ramclustObj$MSMSdata)
    })
    fill <- which(missing > feat.max.missing)
    for(i in fill) {
      fill.current <- which(is.na(ramclustObj$MSMSdata[,i]))
      ramclustObj$MSMSdata[fill.current,i] <- {
        feat.min[i] * 0.5 * rnorm(length(fill.current), mean = 1, sd = 0.05)
      }
    }
  }
  
  ## Step 4 - impute remaining missing values.  
  require(impute)
  # ramclustObj = RC
  MSnames <- dimnames(ramclustObj$MSdata)
  MSMSnames <- dimnames(ramclustObj$MSMSdata)
  d <- ramclustObj$MSdata
  if(!is.null(ramclustObj$MSMSdata)) {
    tmp <- ramclustObj$MSMSdata
    names(tmp) <- paste0(names(tmp), "_MSMS")
    d <- cbind(d, tmp)
    rm(tmp)
  }
  
  d <- impute::impute.knn(
    data = t(d), 
    k = 10,
    rowmax = feat.max.missing,
    colmax = samp.max.missing,
    maxp = 20000,
    rng.seed = 362436069
  )$data
  
  if(!is.null(ramclustObj$MSMSdata)) {
    ms <- t(d[1:(nrow(d)/2),])
    msms <- t(d[((1+(nrow(d)/2)):nrow(d)),])
    dimnames(ms) <- MSnames
    dimnames(msms) <- MSMSnames
    ramclustObj$MSdata <- ms
    ramclustObj$MSMSdata <- msms
  } else {
    ms <- t(d[1:(nrow(d)),])
    dimnames(ms) <- MSnames
    ramclustObj$MSdata <- ms
  }
  
  return(ramclustObj)
}



