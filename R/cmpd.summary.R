#' cmpd.summary
#'
#' a bit of reporting for compounds, quick access summary and plot (if available)
#' @details Reports name, annotation, retention time, number of features in spectrum, median and mean signal intensity, and if interpretMSSpectrum (do.findmain) has been run, plots an annotated MS level spectrum. 
#' 
#' @param ramclustObj ramclustR object to annotate
#' @param cmpd integer. compound number to report.  i.e. 459.  
#' @concept ramlclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @export 


cmpd.summary <- function(ramclustObj = NULL, 
                         cmpd = 1) {
  
  do.feat <- which(ramclustObj$featclus == cmpd)
  # r <- median(
  #   pmax(
  #     cor(  ramclustObj$MSdata_raw[,do.feat],  ramclustObj$MSdata[,do.feat]),
  #     cor(  ramclustObj$MSdata[,do.feat],ramclustObj$MSMSdata[,do.feat]),
  #     cor(ramclustObj$MSMSdata[,do.feat],ramclustObj$MSMSdata[,do.feat])
  #   )
  # )
  
  cat("Name             :",  ramclustObj$cmpd[cmpd], '\n')
  cat("Annotation       :",  ramclustObj$ann[cmpd], '\n')
  cat("retention time   :",  ramclustObj$clrt[cmpd], '\n')
  cat("n features       :", ramclustObj$nfeat[cmpd], '\n')
  cat("median signal    :", median(ramclustObj$SpecAbund[,cmpd], na.rm = TRUE), '\n')
  cat("mean signal      :", mean(ramclustObj$SpecAbund[,cmpd], na.rm = TRUE), '\n')
  #cat("median feature r :", round(r, digits = 3), '\n')
  
  tmp <- ramclustObj$M.ann[[cmpd]]
  plot(tmp[,"mz"], tmp[,"int"], type = "h", 
       ylim = c(0,1.15)*range(tmp[,"int"]),
       main = ramclustObj$cmpd[cmpd],
       xlab = "mz", ylab = "int")
  abline(h = 0, col = grDevices::gray(0.2))
  graphics::text(x = tmp[,1], 
       y = tmp[,2],
       labels = tmp[,"adduct"],
       pos = 3, offset = 1)
  tmp <- tmp[which(tmp$int > 0.03*max(tmp$int)),]
  graphics::text(x = tmp[,1], 
       y = tmp[,2],
       labels = as.character(tmp[,"mz"]),
       pos = 3, offset = 0.2, cex = 0.7, col = grDevices::gray(0.2))
}
