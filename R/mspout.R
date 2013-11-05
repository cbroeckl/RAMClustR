##the script is working, need to incorporate it into the RAMClustR.R function.
##must figure out how to pass parameters
#source("UPLC_C18params.R")

mspoutfun<- function (j, m) {
  sl<-which(RC$featclus==j)
  wm<-vector(length=length(sl))
  if(m==1) {wts<-rowSums(RC$MSdata[,sl])
            for (k in 1:length(sl)) {     
              wm[k]<-weighted.mean(RC$MSdata[,sl[k]], wts)
            }}
  if(m==2) {wts<-rowSums(RC$MSMSdata[,sl])
            for (k in 1:length(sl)) {    
              wm[k]<-weighted.mean(RC$MSMSdata[,sl[k]], wts)
            }}
  mz<-RC$fmz[sl][order(wm, decreasing=TRUE)]
  rt<-RC$frt[sl][order(wm, decreasing=TRUE)]
  wm<-wm[order(wm, decreasing=TRUE)]
  mrt<-mean(rt)
  npeaks<-length(mz)
  for (l in 1:length(mz)) {
    ion<- paste(round(mz[l], digits=4), round(wm[l]))
    if(l==1) {specdat<-ion} 
    if(l>1)  {specdat<-c(specdat, " ", ion)}
  }
  cat(
    paste("Name: C", j, sep=""), '\n',
    paste("SYNON: $:00in-source", sep=""), '\n',
    paste("SYNON: $:04", sep=""), '\n', 
    paste("SYNON: $:05", if(m==1) {CE1} else {CE2}, sep=""), '\n',
    paste("SYNON: $:06", mstype, sep=""), '\n',
    paste("SYNON: $:07", msinst, sep=""), '\n',
    paste("SYNON: $:09", chrominst, sep=""), '\n',
    paste("SYNON: $:10", ionization, sep=""),  '\n',
    paste("SYNON: $:11", msmode, sep=""), '\n',
    paste("SYNON: $:12", colgas, sep=""), '\n',
    paste("SYNON: $:14", msscanrange, sep=""), '\n',
    paste("SYNON: $:16", conevolt, sep=""), '\n',
    paste("Comment: Rt=", round(mrt, digits=2), 
          "  Contributor=\"Colorado State University Proteomics and Metabolomics Facility\"", 
          "  Study=", Experiment, 
          sep=""), '\n',
    paste("Num Peaks:", npeaks), '\n',
    paste(specdat), '\n', '\n', sep="", file=libName, append= TRUE)
} 

for(m in 1:2){mspout}
