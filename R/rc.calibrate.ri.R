#' rc.cmpd.filter.cv
#'
#' extractor for xcms objects in preparation for clustering  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param calibrant.data character vector defining the file path/name to a csv file containing columns including 'rt', and 'ri'.  Alternatively, a data.frame with those columnn names (case sensitive)
#' @param poly.order integer default = 3.  polynomical order used to fit rt vs ri data, and calculate ri for all feature and metabolite rt values. poly.order should be apprciably smaller than the number of calibrant points. 
#' @details This function generates a new slot in the ramclustR object for retention index.  Calibration is performed using a polynomial fit of order poly.order.  It is the user's responsibility to ensure that the number and span of calibrant points is sufficient to calibrate the full range of feature and compound retention times.  i.e. if the last calibration point is at 1000 seconds, but the last eluting peak is at 1300 seconds, the calibration will be very poor for the late eluting compound. 
#' @return  ramclustR object with retention index assigned for features ($fri) and compounds ($clri).   
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

rc.calibrate.ri  <- function(
  ramclustObj=NULL,
  calibrant.data = "",
  poly.order = 3
) { 
  
  if(is.null(ramclustObj)) {
    stop("RAMClustR object must be provided")
  }
  
  if(!is.numeric(poly.order)) {
    stop("poly.order must be an integer value - default = 2")
  }
  
  if(!is.data.frame(calibrant.data)) {
    if(file.exists(calibrant.data)) {
      tmp <- read.csv(calibrant.data, 
                      header = TRUE, 
                      check.names = FALSE, 
                      stringsAsFactors = FALSE,
                      encoding = "UTF-8")
      if(!all(c("ri", "rt") %in% names(tmp))) {
        stop("calibrant data must have column names minimally containing 'rt' and 'ri'")
      }
      calibrant.data <- tmp
    }
  } else {
    if(!all(c("ri", "rt") %in% names(tmp))) {
      stop("calibrant data must have column names minimally containing 'rt' and 'ri'")
    }
  }
  
  calibrant.data$rt <- as.numeric(calibrant.data$rt)
  calibrant.data$ri <- as.numeric(calibrant.data$ri)
  
  pl.data <- calibrant.data[,c("rt", "ri")]
  fit <- lm(ri ~ poly(rt, poly.order), data = pl.data)
  mean.error <- round(mean(100*abs(fit$fitted.values - pl.data$ri) / pl.data$ri), digits = 2)
  nd <- data.frame("rt" = ramclustObj$clrt,
                   "ri" = rep(NA, length(ramclustObj$clrt)))
  clri <- as.vector(stats::predict.lm(fit, newdata = nd))
  type <- c(rep("cal", nrow(pl.data)), rep("cmpd", length(clri)))
  pl.data <- rbind(pl.data, nd)
  ri.fit <- c(fit$fitted.values, rep(NA, length(clri)))
  pl.data <- data.frame(type, pl.data, ri.fit)
  pl.data[((nrow(calibrant.data)+1):nrow(pl.data)),"ri"] <- clri
  pl.data$type <- factor(pl.data$type, levels = c("cmpd", "cal"))
  suppressWarnings(p <- ggplot2::ggplot(pl.data, ggplot2::aes(x=pl.data$rt, y=pl.data$ri)) + 
                     ggplot2::geom_point(ggplot2::aes(shape = type, size = type)) + 
                     ggplot2::geom_line(ggplot2::aes(x=pl.data$rt[1:nrow(calibrant.data)], y=ri.fit[1:nrow(calibrant.data)]), 
                                        data = pl.data[1:nrow(calibrant.data),],
                                        size = 0.8, col = 2) +
                     ggplot2::ggtitle(paste("Calibration: polynomical order =", poly.order)) + 
                     ggplot2::geom_text(x = 1.15*min(pl.data$rt),
                               y = 0.925*max(pl.data$ri), 
                               hjust=0,
                               label=paste("Mean RI error: ", mean.error, "%")) +
                     ggplot2::theme_bw()
  )
  suppressWarnings( print(p) )
  if(max(ramclustObj$clrt) > (max(calibrant.data$rt)* 1.1) | 
     min(ramclustObj$clrt) < (min(calibrant.data$rt)* 0.9)
  ) {
    warning(" -- calibration range cannot accurantely assign RI values outside the calibrant range.")
  }
  
  if(poly.order*3 > nrow(calibrant.data)
  ) {
    warning(" -- poly.order should be set to a value appreciably less than the number of calibrant points.")
  }
  
  ramclustObj$clri <- clri
  
  ## same for features
  nd <- data.frame("rt" = ramclustObj$frt,
                   "ri" = rep(NA, length(ramclustObj$frt)))
  fri <- as.vector(stats::predict.lm(fit, newdata = nd))
  ramclustObj$fri <- fri
  
  ramclustObj$ri.calibrant <- calibrant.data
  return(ramclustObj)
}
