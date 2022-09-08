#' rc.cmpd.filter.blanks
#'
#' used to remove compounds which are found at similar intensity in blank samples.  Only applied after clustering.  see also rc.feature.filter.blanks for filtering at the feature level (only done before clustering).
#'
#' @param ramclustObj ramclustObj containing SpecAbund dataframe.
#' @param qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default. 
#' @param blank.tag see 'qc.tag' , but for blanks to use as background.
#' @param sn numeric defines the ratio for 'signal'.  i.e. sn = 3 indicates that signal intensity must be 3 fold higher in sample than in blanks, on average, to be retained. 
#' @param remove.blanks logical. TRUE by default.  this removes any recognized blanks samples from the SpecAbund sets after they are used to filter contaminant compounds
#' @details This function removes compounds which contain signal in QC samples comparable to blanks. 
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

rc.cmpd.filter.blanks  <- function(ramclustObj=NULL,
                                      qc.tag = "QC",
                                      blank.tag = "blank",
                                      sn = 3,
                                      remove.blanks = TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("please provide a ramclustR Object as input.", '\n')
  }
  
  if(is.null(ramclustObj$SpecAbund)) {
    
  }
  
  if(!is.numeric(sn)) {
    stop("sn must be numeric",'\n')
  }
  
  if(is.null(ramclustObj$SpecAbund)) {
    stop("this data does not contain a SpecAbund dataset", '\n')
  } else {
    d1 <- ramclustObj$SpecAbund
  }
  
  
  ## define QC and Blank samples in each set
  if(length(qc.tag) == 1) {
    qc <- grep(qc.tag[1], ramclustObj$phenoData$sample.names)
    qc <- qc[which(qc <= nrow(d1))]
  } 
  
  if(length(qc.tag) == 2) {
    qc <- grep(qc.tag[1], ramclustObj$phenoData[[qc.tag[2]]])
    qc <- qc[which(qc <= nrow(d1))]
  }
  
  if(length(qc) == 0) {
    stop("no qc samples found. ", '\n')
  }
  
  if(length(blank.tag) == 1) {
    blank <- grep(blank.tag[1], ramclustObj$phenoData$sample.names)
    blank <- blank[which(blank <= nrow(d1))]
  } 
  
  if(length(blank.tag) == 2) {
    blank <- grep(blank.tag[1], ramclustObj$phenoData[[blank.tag[2]]])
    blank <- blank[which(blank <= nrow(d1))]
  }
  
  if(length(blank) == 0) {
    stop("no blank samples found. ", '\n')
  }
  
  ## create logical vector of features to keep
  ## 'good' features should have signal intensity
  ## in pooled QC and/or Samples that is at least 
  ## sn higher than the process blank
  
  #  MS1 mean signal intensities
  if(length(qc) > 1) {
    ms1.qc.mean <- apply(d1[qc,], 2, FUN = "mean", na.rm = TRUE)
  } else {
    ms1.qc.mean <- d1[qc,]
  } 
  if(length(blank)>1) {
    ms1.blank.mean <- apply(d1[blank,], 2, FUN = "mean",  na.rm = TRUE)
  } else {
    ms1.blank.mean <- d1[blank,]
  }
  absent.in.blank <- which(is.nan(ms1.blank.mean))
  
  # filters
  # which signal is at least 3x larger in QC
  keep.ms1.a <- which(ms1.qc.mean/ms1.blank.mean >= sn)
  # which signal is present in QC and absent in blank
  keep.ms1.b <- which((!is.nan(ms1.qc.mean)) & is.nan(ms1.blank.mean))
  # union of these two are keepers
  keep <- union(keep.ms1.a, keep.ms1.b)
  
  length(keep)/ncol(d1)
  cat(round(100*length(keep)/ncol(d1), digits = 1), "% of compounds move forward", '\n', sep = "")
  
  ## define variables that we need to subset
  do <- rep(FALSE, length(ramclustObj))
  rf <- data.frame(index = 1:ncol(d1)); rf <- rf[,0]
  for(i in 1:length(do)) {
    if(any(names(ramclustObj)[i] %in% c("MSdata_raw","MSMSdata_raw"))) {
      next
    }
    if(is.vector(ramclustObj[[i]]) & !is.list(ramclustObj[[i]])) {
      if(length(ramclustObj[[i]]) == ncol(d1)) {
        rf <- data.frame(rf, ramclustObj[[i]])
        names(rf)[ncol(rf)] <- names(ramclustObj[i])
        ramclustObj[[i]] <- ramclustObj[[i]][keep]
      }
    }
    
    if(is.list(ramclustObj[[i]])) {
      if(length(ramclustObj[[i]]) == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][keep]
      }
    }
    
    if(is.data.frame(ramclustObj[[i]])) {
      cat("df", names(ramclustObj)[i], '\n')
      if(dim(ramclustObj[[i]])[2] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][,keep]
      }
      if(dim(ramclustObj[[i]])[1] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][keep,]
      }
    }
    
    if(is.matrix(ramclustObj[[i]])) {
      cat("ma", names(ramclustObj)[i], '\n')
      if(dim(ramclustObj[[i]])[2] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][,keep]
      }
      if(dim(ramclustObj[[i]])[1] == ncol(d1)) {
        ramclustObj[[i]] <- ramclustObj[[i]][keep,]
      }
    }
    
    if(is.matrix(ramclustObj[[i]])) {
      
    }
  }
  ramclustObj$feature.filter.blanks <- rf[-keep,]
  
  if(remove.blanks){
    ramclustObj$SpecAbund <- ramclustObj$SpecAbund[-blank,]
    if(!is.null(ramclustObj$phenoData)) {
      ramclustObj$phenoData <- ramclustObj$phenoData[-blank,]
    }
  }
  
  ramclustObj$history$cmpd.filter.blanks <- {
    paste0(
      "Compounds which failed to demonstrate signal intensity of at least ",
      sn, " fold greater in QC samples than in blanks were removed from the SpecAbund dataset. ", 
      ncol(d1) - length(keep)," of ", ncol(d1), " compounds were removed."
    )
  }  
  
  cat(ramclustObj$history$cmpd.filter.blanks)
  
  return(ramclustObj)
}
