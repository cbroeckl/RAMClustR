#' rc.normalize.batch.qc
#'
#' extractor for xcms objects in preparation for clustering  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param normalize character: either "none", "TIC", "quantile", or "batch.qc" normalization of feature intensities.  see batch.qc overview in details. 
#' @param qc.inj.range integer: how many injections around each injection are to be scanned for presence of QC samples when using batch.qc normalization?  A good rule of thumb is between 1 and 3 times the typical injection span between QC injections.  i.e. if you inject QC ever 7 samples, set this to between 7 and 21.  smaller values provide more local precision but make normalization sensitive to individual poor outliers (though these are first removed using the boxplot function outlier detection), while wider values provide less local precision in normalization but better stability to individual peak areas.
#' @param batch integer vector with length equal to number of injections in xset or csv file
#' @param order integer vector with length equal to number of injections in xset or csv file
#' @param qc logical vector with length equal to number of injections in xset or csv file.  
#' @details This function offers normalization by run order, batch number, and QC sample signal intensity.
#' @details Each input vector should be the same length, and equal to the number of samples in the $MSdata set.
#' @details Input vector order is assumed to be the same as the sample order in the $MSdata set.  
#' @return  ramclustR object with normalized data.   
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

rc.normalize.batch.qc  <- function(ramclustObj=NULL,
                                   qc.inj.range = 20,
                                   order = NULL,
                                   batch = NULL,
                                   qc = NULL
) {
  
  ## CHECKS
  if(is.null(ramclustObj)) {
    stop('existing ramclustObj required as input', '\n', 
         '       see rc.get.xcms.data function for one approach to do so', '\n')
  }
  
  if(is.null(order)) {
    stop("'order' must be set as a vector of integer values (i.e. c(10,13,7)) for this normalization approach", '\n')
  }
  
  if(is.null(batch)) {
    warning("batch = NULL; assumption of a single batch experiment is applied", '\n')
    batch <- rep(1, nrow(ramclustObj$MSdata))
  }
  
  if(is.null(qc)) {
    warning("qc = NULL; run order correction will be applied without QC samples.", '\n',
            "     This approach assumes randomized run order. lack of QC samples", '\n', 
            "     reduces confidence that run order trends are purely analtyical.", '\n')
    qc <- rep(TRUE, nrow(ramclustObj$MSdata))
  }
  
  if(!is.logical(qc)) {
    stop("qc must be a logical vector", '\n')
  }
  
  if(!all.equal(length(batch), length(qc), length(order), nrow(ramclustObj$MSdata))) {
    stop("all lengths must be identical and are not: ", '\n',
         "  length(batch) = ", length(batch), '\n', 
         "  length(order) = ", length(order), '\n',
         "  length(qc) = ", length(qc), '\n',
         "  number of injections = ", nrow(data1), '\n')
  }
  
  
  data1 <- ramclustObj$MSdata
  mslev <- 1
  if(!is.null(ramclustObj$MSMSdata)) {
    data2 <- ramclustObj$MSMSdata
    mslev <- 2
  }
  
  ndf <- data.frame(batch, order, qc)
  new.ord <- order(ndf$order)
  ndf <- ndf[new.ord,]
  data1 <- data1[new.ord,]
  new.ord <- order(ndf$batch)
  ndf <- ndf[new.ord,]
  data1 <- data1[new.ord,]
  batch <- ndf[,"batch"]
  qc <- ndf[,"qc"]
  order <- ndf[,"order"]
  
  
  data1raw <- data1
  if(mslev > 1) {data2raw <- data2}
  
  ramclustObj$history$normalization.batch.qc <- paste(
    " Features were normalized to nearby QC samples on a feature-by-feature basis using the 'batch.qc' option ",
    "with qc.inj.range = ", qc.inj.range, ".", 
    " QC normalization was applied to ", length(order), " injections in ", length(unique(batch)), 
    " bactches, and normization was based on ", length(which(qc)), " recognized QC samples.", 
    sep = "") 
  
  pdf(file = "norm.plots.pdf", height = 4, width = 9)
  max.ratio <- 4
  for(z in 1:ncol(data1)) {
    # z <- sample(1:ncol(data1), 1)
    tmp <- data1[,z]  
    featmed <- mean(tmp[qc])
    tmpn <- tmp
    
    for( i in unique(batch)) {
      
      do <- which (batch == i)
      doqc <- which(batch == i & qc)
      
      ## use only 'typical' QC sample values from the given batch
      ## outliers are detected using the standard boxplot definition (1.5 * the interquartile range)
      # out <- boxplot(tmp[doqc], plot = FALSE)$out
      sds <- 1.96
      lcl <- mean(tmp[doqc]) - (sds*sd(tmp[doqc]))
      ucl <- mean(tmp[doqc]) + (sds*sd(tmp[doqc]))
      keep <- which(tmp[doqc] > lcl & tmp[doqc] < ucl)
      #if(length(out)>0) doqc <- doqc[!(names(doqc) %in% names(out))]
      if(length(keep)>0) doqc <- doqc[keep]
      
      batchmed <- mean(tmpn[doqc])
      f <- batchmed / featmed
      if(is.na(f)) next
      if(abs(log2(f)) > max.ratio)  {
        if(f > 1) {f <- max.ratio}
        if(f < 1) {f <- 1/max.ratio}
      }
      tmpn[do] <- tmp[do]/f
      
      tmpnqc <- tmpn
      
      # cat("batch:", i, " raw   CV =", sd(tmp[doqc])/mean(tmp[doqc]), '\n')
      # tmpna <- tmpn
      ## local QC adjustment here: 
      
      for(x in do) {
        
        batchmed <- mean(tmpn[doqc])
        ## try injection spacing weighted mean instead
        space <- abs(x - doqc)
        use <- which(space <= qc.inj.range)
        wts <- 1/space[use]
        wts <- wts/sum(wts)
        localmed <- weighted.mean(tmpnqc[doqc[use]], weights = wts)
        
        # localmed <- median(tmpn[doqc[which(abs(x - doqc) <= qc.inj.range*y)]])
        if(is.na(localmed)){
          for(y in 2:5) {
            # median
            # localmed <- median(tmpn[doqc[which(abs(x - doqc) <= qc.inj.range*y)]])
            mean
            use <- which(space <= qc.inj.range*y)
            wts <- 1/space[use]
            wts <- wts/sum(wts)
            localmed <- weighted.mean(tmpnqc[doqc[use]], weights = wts)
            if(!is.na(localmed)) break
          }
        }
        if(is.na(localmed)) {localmed <- batchmed}
        f <- localmed/batchmed
        if(is.na(f)) next
        
        if(abs(log2(f)) > max.ratio) {
          # f <- batchmed / featmed
          if(f > 1) {f <- max.ratio}
          if(f < 1) {f <- 1/max.ratio}
        }
        
        tmpn[x] <- tmpnqc[x] / f
        rm(localmed); rm(f)
      }
      # cat("batch:", i, " normb CV =", sd(tmpn[doqc])/mean(tmpn[doqc]), '\n')
    }
    data1[,z] <- tmpn
    par(mfrow = c(1,2))
    tmp <- as.vector(tmp)
    tmpn <- as.vector(tmpn)
    if(all(is.na(tmp))) {
      plot(0,0, pch = "", main = "no peaks detected")
      plot(0,0, pch = "", main = "no peaks detected")
    } else {
      plot(tmp, col = batch, cex = (qc + 1)/2, 
           ylim = c(0.9,1.11)*range(tmp, na.rm = TRUE), 
           main = paste("all:", round(sd(tmp, na.rm = TRUE)/mean(tmp, na.rm = TRUE), digits = 2), '\n',
                        "qc:", round(sd(tmp[qc], na.rm = TRUE)/mean(tmp[qc], na.rm = TRUE), digits = 2)))
      plot(tmpn, col = batch, pch = 19, cex = (qc + 1)/2, 
           ylim = c(0.9,1.11)*range(tmp, na.rm = TRUE), 
           main = paste("all:", round(sd(tmpn, na.rm = TRUE)/mean(tmpn, na.rm = TRUE), digits = 2), '\n',
                        "qc:", round(sd(tmpn[qc], na.rm = TRUE)/mean(tmpn[qc], na.rm = TRUE), digits = 2)))
      
    }
    rm(tmp); rm(tmpn); rm(tmpnqc); gc()
    
  }
  
  if(mslev == 2) {
    for(z in 1:ncol(data2)) {
      # z <- sample(1:ncol(data2), 1)
      tmp <- data2[,z]
      featmed <- mean(tmp[qc])
      tmpn <- tmp
      
      for( i in unique(batch)) {
        
        do <- which (batch == i)
        doqc <- which(batch == i & qc)
        names(doqc) <- names(tmp[doqc])
        
        ## use only 'typical' QC sample values from the given batch
        ## outliers are detected using the standard boxplot definition (1.5 * the interquartile range)
        #out <- boxplot(tmp[doqc], plot = FALSE)$out
        #if(length(out)>0) doqc <- doqc[!(names(doqc) %in% names(out))]
        lcl <- mean(tmp[doqc]) - (sds*sd(tmp[doqc]))
        ucl <- mean(tmp[doqc]) + (sds*sd(tmp[doqc]))
        keep <- which(tmp[doqc] > lcl & tmp[doqc] < ucl)
        if(length(keep)>0) doqc <- doqc[keep]
        
        batchmed <- mean(tmpn[doqc])
        f <- batchmed / featmed
        if(is.na(f)) next
        if(abs(log2(f)) > max.ratio)  {
          if(f > 1) {f <- max.ratio}
          if(f < 1) {f <- 1/max.ratio}
        }
        tmpn[do] <- tmp[do]/f
        
        tmpnqc <- tmpn
        
        # cat("batch:", i, " raw   CV =", sd(tmp[doqc])/mean(tmp[doqc]), '\n')
        # tmpna <- tmpn
        ## local QC adjustment here: 
        
        for(x in do) {
          
          batchmed <- mean(tmpn[doqc])
          ## try injection spacing weighted mean instead
          space <- abs(x - doqc)
          use <- which(space <= qc.inj.range)
          wts <- 1/space[use]
          wts <- wts/sum(wts)
          localmed <- weighted.mean(tmpnqc[doqc[use]], weights = wts)
          
          # localmed <- median(tmpn[doqc[which(abs(x - doqc) <= qc.inj.range*y)]])
          if(is.na(localmed)){
            for(y in 2:5) {
              # median
              # localmed <- median(tmpn[doqc[which(abs(x - doqc) <= qc.inj.range*y)]])
              mean
              use <- which(space <= qc.inj.range*y)
              wts <- 1/space[use]
              wts <- wts/sum(wts)
              localmed <- weighted.mean(tmpnqc[doqc[use]], weights = wts)
              if(!is.na(localmed)) break
            }
          }
          if(is.na(localmed)) {localmed <- batchmed}
          f <- localmed/batchmed
          if(is.na(f)) next
          if(abs(log2(f)) > max.ratio) {
            # f <- batchmed / featmed
            if(f > 1) {f <- max.ratio}
            if(f < 1) {f <- 1/max.ratio}
          }
          tmpn[x] <- tmpnqc[x] / f
          rm(localmed); rm(f)
        }
        # cat("batch:", i, " normb CV =", sd(tmpn[doqc])/mean(tmpn[doqc]), '\n')
      }
      data2[,z] <- tmpn
      par(mfrow = c(1,2))
      tmp <- as.vector(tmp)
      tmpn <- as.vector(tmpn)
      if(all(is.na(tmp))) {
        plot(0,0, pch = "", main = "no peaks detected")
        plot(0,0, pch = "", main = "no peaks detected")
      } else {
        plot(tmp, col = batch, cex = (qc + 1)/2, 
             ylim = c(0.9,1.11)*range(tmp, na.rm = TRUE), 
             main = paste("all:", round(sd(tmp, na.rm = TRUE)/mean(tmp, na.rm = TRUE), digits = 2), '\n',
                          "qc:", round(sd(tmp[qc], na.rm = TRUE)/mean(tmp[qc], na.rm = TRUE), digits = 2)))
        plot(tmpn, col = batch, pch = 19, cex = (qc + 1)/2, 
             ylim = c(0.9,1.11)*range(tmp, na.rm = TRUE), 
             main = paste("all:", round(sd(tmpn, na.rm = TRUE)/mean(tmpn, na.rm = TRUE), digits = 2), '\n',
                          "qc:", round(sd(tmpn[qc], na.rm = TRUE)/mean(tmpn[qc], na.rm = TRUE), digits = 2)))
      }
      rm(tmp); rm(tmpn); rm(tmpnqc); gc()
      
    }
    
  }
  
  dev.off()
  
  
  
  return(ramclustObj)
}

