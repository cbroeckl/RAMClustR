#' write.single.msp
#'
#' writes single msp files compatible with MSfinder or other spectral interpretation tools.   
#'
#' @param ramclustObj your ramclustR output R object. 
#' @param mzdec integer: number of decimal places used in printing m/z values
#' 
#' @details writes individual msp files for each compound to a new 'spectra/msp' directory
#' @return   - 'spectra/msp' directory is created in the working directory.  MS level 1 spectra are labelled with compound name, while MS2 spectra have an appended "MS2" in the filename.
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @author Corey Broeckling
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms'
#' @export


write.single.msp<-function(
  ramclustObj = RC,
  mzdec = 4
) {
  if (!dir.exists("spectra")) {
    dir.create("spectra")
  }
  dir.create("spectra/msp")
  
  ########
  # define ms levels, used several times below
  mslev <- as.integer(as.numeric(as.character(ExpDes[[2]][which(row.names(ExpDes[[2]]) == "MSlevs"),1])))
  
  
  ########
  # write msp spectra
  cat(paste("writing msp formatted spectra...", '\n'))
  for (m in 1:as.numeric(mslev)){
    for (j in 1:max(ramclustObj$featclus)) {
      if(m == 1) {
        libName<-paste("spectra/msp/", ramclustObj$cmpd[j], ".msp", sep="")
      } else {
        libName<-paste("spectra/msp/", ramclustObj$cmpd[j], "_MS", m, ".msp", sep="")
      }

      #print(paste(j,"_", sep=""))
      sl<-which(ramclustObj$featclus==j)
      wm<-vector(length=length(sl))
      if(m==1) {wts<-rowSums(ramclustObj$MSdata[,sl, drop=FALSE])
      for (k in 1:length(sl)) {     
        wm[k]<-weighted.mean(ramclustObj$MSdata[,sl[k]], wts)
      }}
      if(m==2) {wts<-rowSums(ramclustObj$MSMSdata[,sl, drop=FALSE])
      for (k in 1:length(sl)) {    
        wm[k]<-weighted.mean(ramclustObj$MSMSdata[,sl[k]], wts)
      }}
      mz<-round(ramclustObj$fmz[sl][order(wm, decreasing=TRUE)], digits=mzdec)
      rt<-ramclustObj$frt[sl][order(wm, decreasing=TRUE)]
      wm<-round(wm[order(wm, decreasing=TRUE)])
      mrt<-mean(rt)
      npeaks<-length(mz)
      specdat<-""
      for (l in 1:length(mz)) {
        specdat<-paste(specdat, mz[l], " ", wm[l], '\n', sep="")
      }
      cat(
        paste("Name: ", ramclustObj$cmpd[j], sep=""), '\n',
        paste("SYNON: $:00in-source", sep=""), '\n',
        paste("SYNON: $:04", sep=""), '\n', 
        paste("SYNON: $:05", if(m==1) {ExpDes[[2]]["CE1", 1]} else {ExpDes$instrument["CE2", "InstVals"]}, sep=""), '\n',
        paste("SYNON: $:06", ExpDes[[2]]["mstype", 1], sep=""), '\n',           #mstype
        paste("SYNON: $:07", ExpDes[[2]]["msinst", 1], sep=""), '\n',           #msinst
        paste("SYNON: $:09", ExpDes[[2]]["chrominst", 1], sep=""), '\n',        #chrominst
        paste("SYNON: $:10", ExpDes[[2]]["ionization", 1], sep=""),  '\n',      #ionization method
        paste("SYNON: $:11", ExpDes[[2]]["msmode", 1], sep=""), '\n',           #msmode
        if(any(row.names(ExpDes[[2]])=="colgas")) {
          paste("SYNON: $:12", ExpDes[[2]]["colgas", 1], '\n', sep="") 
        },          #collision gas
        paste("SYNON: $:14", ExpDes[[2]]["msscanrange", 1], sep=""), '\n',      #ms scanrange
        if(any(row.names(ExpDes[[2]])=="conevolt")) {
          paste("SYNON: $:16", ExpDes[[2]]["conevolt", 1], '\n', sep="")
        },         #conevoltage
        paste("Comment: Rt=", round(mrt, digits=2), 
              "  Contributor=", ExpDes[[1]]["Contributor", 1], 
              "  Study=", ExpDes[[1]]["Experiment", 1], 
              sep=""), '\n',
        paste("Num Peaks:", npeaks), '\n',
        paste(specdat), '\n', sep="", file=libName, append= TRUE)
    }
  }
  cat(paste('\n', "msp files complete", '\n')) 
}