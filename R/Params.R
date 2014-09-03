
defineExperiment<-function()
{

  if(file.exists(system.file('params/paramsets.Rdata', package = "RAMClustR"))) {
    
	load(system.file('params/paramsets.Rdata', package = "RAMClustR"))} else {
	load(system.file('params/defparamsets.Rdata', package = "RAMClustR"))   }
 

 
  platforms<-names(paramsets) 
  
  ExpVals<-c(Experiment=".",
             Species=".",
             Sample=".",
             Contributor=".",
             platform="undefined")
  
  VarDesc<-c("experiment name, no spaces",
             "species name",
             "sample type",
             "individual and/or organizational affiliation",
             paste(platforms, sep=" ", collapse=" "))
  
  Experiment<-data.frame(ExpVals,VarDesc, stringsAsFactors=FALSE)
  
  suppressWarnings(design<-edit(Experiment))
  
  platform<-platforms[grep(as.character(design["platform",1]), platforms, ignore.case=TRUE)]

  suppressWarnings(instrument<-edit(paramsets[[as.character(platform)]]))
  
  if(!(instrument["saveAs",1] %in% platforms) & instrument["saveAs",1] != "") 
                                  {newParamset<-list(instrument)
                                  names(newParamset)<-instrument["saveAs",1]
                                  paramsets<-c(paramsets, newParamset)
                                  save(paramsets, file=paste(libPaths(), "/RAMClustR/params/paramsets.Rdata", sep=""))
                                  }
  names(instrument)<-"InstVals"
  exp.pars<-list(design, instrument)
  names(exp.pars)<-c("design", "instrument")
  return(exp.pars)
}
