#' rc.expand.sample.names
#'
#' turn concatenated sample names into factors  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param delim what delimiter should be used to separate names into factors?  '-' by default
#' @param factor.names logical or character vector.  if TRUE, user will enter names one by on in console.  If character vector (i.e. c("trt", "time")) names are assigned to table
#' @details THis function only works on newer format ramclustObjects with a $phenoData slot.
#' @details This function will split sample names by a delimiter, and enable users to name factors
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

rc.expand.sample.names <- function(
  ramclustObj = NULL,
  delim = "-",
  factor.names = TRUE
) {
  
  params <- c(
    "delim" = delim
  )


  if(!is.null(ramclustObj$phenoData$sample.names)) {
    sn <- as.character(ramclustObj$phenoData$sample.names)
  }
  if(!is.null(ramclustObj$phenoData$sample.names.sn)) {
    sn <- as.character(ramclustObj$phenoData$sample.names.sn)
  }
  
  if(!any(ls()=="sn")) {
    stop('missing sample names in phenoData slot', '\n')
  }
  des <-  strsplit(sn, delim)
  l <- sapply(1:length(des), FUN = function(x) {
    length(des[[x]])
  })
  if(length(table(l)) != 1) {
    cat("delimited sample names have variable lengths ranging from: ", '\n',
        range(l)[1], " to ", range(l)[2], '\n')
    ch <- which(l != median(l))
    stop("please fix sample names:", '\n', paste(" ", sn[ch], sep = '\n'))
  }
  des <- data.frame(t(data.frame(des, check.names = FALSE)), 
                    stringsAsFactors = FALSE, check.names = FALSE)
  
  if(is.logical(factor.names)) {
    if(factor.names){
      for(x in 1:ncol(des)) {
        cat(
          "column",x, "variables:",'\n',
          unique(des[,x]), '\n')
        fn <- readline(prompt=cat("Type name and press [enter] to continue:", '\n'))
        dimnames(des)[[2]][x] <- fn
      }
    }
  }
  
  if(is.character(factor.names)) {
    if(length(factor.names) != ncol(des)) {
      stop(length(factor.names), "factor names and", ncol(des), "factors - please correct", '\n')
    }
    dimnames(des)[[2]] <- factor.names
  }
  rn <- row.names(ramclustObj$phenoData)
  ramclustObj$phenoData <- cbind(ramclustObj$phenoData, des[,1:ncol(des)])
  
  if(is.null(ramclustObj$params)) {ramclustObj$params <- list()}
  ramclustObj$params$rc.expand.sample.names <- params
  
  row.names(ramclustObj$phenoData) <- rn
  return(ramclustObj)
  
}

