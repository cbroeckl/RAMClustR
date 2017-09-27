#' ramclustR
#'
#' Main clustering function for grouping features based on their analytical behavior.  
#'
#' @param xcmsObj xcmsObject: containing grouped feature data for clustering by ramclustR
#' @param ms filepath: optional csv input. Features as columns, rows as samples. Column header mz_rt
#' @param idmsms filepath: optional idMSMS / MSe csv data.  same dim and names as ms required
#' @param idMSMStag character: character string in 'taglocation' to designat idMSMS / MSe files e.g. "02.cdf"
#' @param taglocation character: "filepaths" by default, "phenoData[,1]" is another option. referse to xcms slot
#' @param featdelim character: how feature mz and rt are delimited in csv import column header e.g. ="-"
#' @param timepos integer: which position in delimited column header represents the retention time (csv only)
#' @param st numeric: sigma t - time similarity decay value 
#' @param sr numeric: sigma r - correlational similarity decay value
#' @param maxt numeric: maximum time difference to calculate retention similarity for - all values beyond this are assigned similarity of zero
#' @param deepSplit logical: controls how agressively the HCA tree is cut - see ?cutreeDynamicTree
#' @param blocksize integer: number of features (scans?) processed in one block  =1000,
#' @param mult numeric: internal value, can be used to influence processing speed/ram usage
#' @param hmax numeric: precut the tree at this height, default 0.3 - see ?cutreeDynamicTree
#' @param sampNameCol integer: which column from the csv file contains sample names?
#' @param collapse logical: reduce feature intensities to spectrum intensities?
#' @param usePheno logical: tranfer phenotype data from XCMS object to SpecAbund dataset?
#' @param mspout logical: write msp formatted specta to file?
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param normalize character: either "none", "TIC", or "quantile" normalization of feature intensities
#' @param minModuleSize integer: how many features must be part of a cluster to be returned? default = 2
#' @param linkage character: heirarchical clustering linkage method - see ?hclust
#' @param mzdec integer: number of decimal places used in printing m/z values
#' @param cor.method character: which correlational method used to calculate 'r' - see ?cor
#'
#' @details Main clustering function output - see citation for algorithm description of vignette('RAMClustR') for a walk through
#' @return   featclus: integer vector of cluster membership for each feature
#' @return   frt: feature retention time, in whatever units were fed in (xcms uses seconds, by default)
#' @return   fmz: feature retention time, reported in number of decimal points selected in ramclustR function
#' @return   xcmsOrd: the original XCMS (or csv) feature order for cross referencing, if need be
#' @return   clrt: cluster retention time
#' @return   clrtsd: retention time standard deviation of all the features that comprise that cluster
#' @return   nfeat: number of features in the cluster
#' @return   nsing: number of 'singletons' - that is the number of features which clustered with no other feature
#' @return   ExpDes: the experimental design object used when running ramclustR.  List of two dataframes. 
#' @return   cmpd: compound name.  C#### are assigned in order of output by dynamicTreeCut.  Compound with the most features is classified as C0001...
#' @return   ann: annotation.  By default, annotation names are identical to 'cmpd' names.  This slot is a placeholder for when annotations are provided
#' @return   MSdata:  the MSdataset provided by either xcms or csv input
#' @return   MSMSdata: the (optional) MSe/idMSMS dataset provided be either xcms or csv input
#' @return   SpecAbund: the cluster intensities after collapsing features to clusters
#' @return   SpecAbundAve: the cluster intensities after averaging all samples with identical sample names
#' @return   - 'spectra' directory is created in the working directory.  In this directory a .msp is (optionally) created, which contains the spectra for all compounds in the dataset following clustering.  if MSe/idMSMS data are provided, they are listed witht he same compound name as the MS spectrum, with the collision energy provided in the ExpDes object provided to distinguish low from high CE spectra. 
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @author Corey Broeckling
#' @export

ramclustR<- function(  xcmsObj=NULL,
                       ms=NULL, 
                       idmsms=NULL,
                       taglocation="filepaths",
                       MStag=NULL,
                       idMSMStag=NULL, 
                       featdelim="_", 
                       mzpos=1, 
                       timepos=2, 
                       st=NULL, 
                       sr=NULL, 
                       maxt=NULL, 
                       deepSplit=FALSE, 
                       blocksize=2000,
                       mult=5,
                       hmax=NULL,
                       sampNameCol=1,
                       collapse=TRUE,
                       usePheno=TRUE,
                       mspout=TRUE, 
                       ExpDes=NULL,
                       normalize="TIC",
                       minModuleSize=2,
                       linkage="average",
                       mzdec=4,
                       cor.method="pearson"
) {
  
  ########
  # load packages
  require(xcms, quietly=TRUE)
  require(ff, quietly=TRUE)
  require(fastcluster, quietly=TRUE)
  require(dynamicTreeCut, quietly=TRUE)
  
  ########
  # define ms levels, used several tims below
  mslev <- as.integer(as.numeric(as.character(ExpDes[[2]][nrow(ExpDes[[2]]),1])))
  
  
  ########
  # do some checks to make sure we have everything we need before proceeding
  if(is.null(xcmsObj) & is.null(ms))  {
    stop("you must select either 
          1: an MS dataset with features as columns 
             (__)one column may contain sample names, defined by sampNameCol) 
          2: an xcmsObj. If you choose an xcms object, set taglocation: 'filepaths' by default
            and both MStag and idMSMStag")}
  
  if(!is.null(xcmsObj) & mslev==2 & any(is.null(MStag), is.null(idMSMStag), is.null(taglocation)))
  {stop("you must specify the the MStag, idMSMStag, and the taglocations")}
  
  if( normalize!="none"  & normalize!="TIC" & normalize!="quantile") {
    stop("please selected either 'none', 'TIC', or 'quantile' for the normalize setting")}
  

  a<-Sys.time()   
  
  ########
  # set default parameters of not defined in function call
  
  if(is.null(hmax)) {hmax<-0.3}
  
  
  ########
  # if csv input of MS data, do this: 
  if(!is.null(ms)){
    if(is.null(st)) stop("please specify st: 
      a recommended starting point is half the value of 
      your average chromatographic peak width at half max (seconds)")
    if(is.null(sr)) sr<-0.5
    if(is.null(maxt)) maxt<-60
    MSdata<-read.csv(file = ms, header=TRUE, check.names=FALSE)
    if(!is.null(idmsms)){
      MSMSdata<-read.csv(file = idmsms, header=TRUE, check.names=FALSE)}
    if(is.null(idmsms)) { MSMSdata<-MSdata}
    if(is.null(sampNameCol)) {featcol<-1:ncol(MSdata)} else {
      featcol<-setdiff(1:(ncol(MSdata)), sampNameCol)}
    if(is.null(sampNameCol)) {featcol<-1:ncol(MSdata)} else {
      featcol<-setdiff(1:(ncol(MSdata)), sampNameCol)}
    sampnames<-MSdata[,sampNameCol]
    data1<-as.matrix(MSdata[,featcol])
    dimnames(data1)[[1]]<-MSdata[,sampNameCol]
    dimnames(data1)[[2]]<-names(MSdata[,featcol])
    data2<-as.matrix(MSMSdata[,featcol])
    dimnames(data2)[[1]]<-MSMSdata[,sampNameCol]
    dimnames(data2)[[2]]<-names(MSMSdata[,featcol])
    if(!all(dimnames(data1)[[2]]==dimnames(data2)[[2]])) 
    {stop("the feature names of your MS and idMSMS data are not identical")}
    
    if(!all(dimnames(data1)[[1]]==dimnames(data2)[[1]])) 
    {stop("the order and names of your MS and idMSMS data sample names are not identical")}
    
    rtmz<-matrix(
      unlist(
        strsplit(dimnames(data1)[[2]], featdelim)
      ), 
      byrow=TRUE, ncol=2)
    times<-as.numeric(rtmz[,timepos])
    mzs<-as.numeric(rtmz[,which(c(1:2)!=timepos)])
    rm(rtmz)
  }
  
  ########
  # if xcms object as input, do this:
  if(!is.null(xcmsObj)){
    if(!class(xcmsObj)=="xcmsSet")
    {stop("xcmsObj must reference an object generated by XCMS, of class 'xcmsSet'")}
    
    if(is.null(st)) st<-round(median(xcmsObj@peaks[,"rtmax"]-xcmsObj@peaks[,"rtmin"])/2, digits=2)
    if(is.null(sr)) sr<-0.5
    if(is.null(maxt)) maxt<-st*20
    
    sampnames<-row.names(xcmsObj@phenoData)
    data12<-groupval(xcmsObj, value="into")
    
    if(taglocation=="filepaths" & !is.null(MStag)) 
    { msfiles<-grep(MStag, xcmsObj@filepaths, ignore.case=TRUE)
    msmsfiles<-grep(idMSMStag, xcmsObj@filepaths, ignore.case=TRUE)
    if(length(intersect(msfiles, msmsfiles)>0)) 
    {stop("your MS and idMSMStag values do not generate unique file lists")}
    if(length(msfiles)!=length(msmsfiles)) 
    {stop("the number of MS files must equal the number of MSMS files")}
    data1<-t(data12[,msfiles])
    row.names(data1)<-sampnames[msfiles]
    data2<-t(data12[,msmsfiles])
    row.names(data2)<-sampnames[msmsfiles]  ##this may need to be changed to dimnames..
    times<-round(xcmsObj@groups[,"rtmed"], digits=3)
    if(any(is.na(times))) {
      do<-which(is.na(times))
      for(x in 1:length(do)) {
        times[do[x]]<-  as.numeric((xcmsObj@groups[do[x],"rtmin"]+ xcmsObj@groups[do[x],"rtmax"])/2)
      }
    }
    # if(any(is.na(times))) {stop("na values still present")} else {print("NAs removed")}
    mzs<-round(xcmsObj@groups[,"mzmed"], digits=4)
    } else {
      data1<-t(data12)
      msfiles<-1:nrow(data1)
      data2<-t(data12)
      times<-round(xcmsObj@groups[,"rtmed"], digits=3)
      mzs<-round(xcmsObj@groups[,"mzmed"], digits=4) 
      if(any(is.na(times))) {
        do<-which(is.na(times))
        for(x in 1:length(do)) {
          times[do[x]]<-  as.numeric((xcmsObj@groups[do[x],"rtmin"]+ xcmsObj@groups[do[x],"rtmax"])/2)
        }
      }
      #if(any(is.na(times))) {stop("na values still present")} else {print("NAs removed")}
    }
  }
  
  
  ########
  # ensure that we have all numeric values in the dataset. 
  # uses a noise addition 'jitter' around minimum values with missing data points.
  # this is mostly necessary for csv input, where other programs may not have used a 'fillPeaks' like step
  rpl1<-unique(c(which(is.na(data1)), which(is.nan(data1)), which(is.infinite(data1)), which(data1==0)))
  rpl2<-unique(c(which(is.na(data2)), which(is.nan(data2)), which(is.infinite(data2)), which(data2==0)))
  if(length(rpl1)>0) {data1[rpl1]<-abs(jitter(rep(min(data1, na.rm=TRUE), length(rpl1) ), amount=min(data1/100, na.rm=TRUE)))}
  if(length(rpl2)>0) {data2[rpl2]<-abs(jitter(rep(min(data2, na.rm=TRUE), length(rpl2) ), amount=min(data2/100, na.rm=TRUE)))}
  #if(length(rpl1)>0) {data1[rpl1]<-abs(jitter(rep(min(data1, na.rm=TRUE), length(rpl1) ), amount=min(data1/100, na.rm=TRUE)))}
  #if(length(rpl2)>0) {data2[rpl2]<-abs(jitter(rep(min(data2, na.rm=TRUE), length(rpl2) ), amount=min(data2/100, na.rm=TRUE)))}
  data1[which(data1<0)]<-abs(data1[which(data1<0)])
  data2[which(data2<0)]<-abs(data2[which(data2<0)])
  
  ########
  # Optional normalization of data, either Total ion signal or quantile
  
  if(normalize=="TIC") {
    data1<-(data1/rowSums(data1))*mean(rowSums(data1), na.rm=TRUE)
    data2<-(data2/rowSums(data2))*mean(rowSums(data2), na.rm=TRUE)
  }
  if(normalize=="quantile") {
    data1<-t(preprocessCore::normalize.quantiles(t(data1)))
    data2<-t(preprocessCore::normalize.quantiles(t(data2)))	
  }
  
  ########
  # data organization and parsing 
  # sort rt vector and data by retention time
  xcmsOrd<-order(times)
  data1<-data1[,order(times)]
  data2<-data2[,order(times)]
  mzs<-mzs[order(times)]
  times<-times[order(times)]
  featnames<-paste(mzs, "_", times, sep="")
  dimnames(data1)[[2]]<-featnames
  dimnames(data2)[[2]]<-featnames
  
  ########
  # establish some constants for downstream processing
  n<-ncol(data1)
  vlength<-(n*(n-1))/2
  nblocks<-floor(n/blocksize)
  
  ########
  # set off ff matrix system for holding data. 
  # manages RAM demands a bit.  
  ffmat<-ff(vmode="double", dim=c(n, n), initdata = 0) ##reset to 1 if necessary
  gc()
  #Sys.sleep((n^2)/10000000)
  #gc()
  
  ########
  # make list of all row and column blocks to evaluate
  eval1<-expand.grid(0:nblocks, 0:nblocks)
  names(eval1)<-c("j", "k") #j for cols, k for rows
  eval1<-eval1[which(eval1[,"j"]<=eval1[,"k"]),] #upper triangle only
  bl<-nrow(eval1)
  cat('\n', paste("calculating ramclustR similarity: nblocks = ", bl))
  cat('\n', "finished:")
  
  
  ########
  # Define the RCsim function used to calculate feature similarities on selected blocks of data
  RCsim<-function(bl)  {
    cat(bl,' ')
    j<-eval1[bl,"j"]  #columns
    k<-eval1[bl,"k"]  #rows
    startc<-min((1+(j*blocksize)), n)
    if ((j+1)*blocksize > n) {
      stopc<-n} else {
        stopc<-(j+1)*blocksize}
    startr<-min((1+(k*blocksize)), n)
    if ((k+1)*blocksize > n) {
      stopr<-n} else {
        stopr<-(k+1)*blocksize}
    if(startc<=startr) { 
      mint<-min(abs(outer(range(times[startr:stopr]), range(times[startc:stopc]), FUN="-")))
      if(mint<=maxt) {
        temp1<-round(exp(-(( (abs(outer(times[startr:stopr], times[startc:stopc], FUN="-"))))^2)/(2*(st^2))), 
                     
                     digits=20 )

        temp2<-round (exp(-((1-(pmax(  cor(data1[,startr:stopr], data1[,startc:stopc], method=cor.method),
                                       cor(data1[,startr:stopr], data2[,startc:stopc], method=cor.method),
                                       cor(data2[,startr:stopr], data2[,startc:stopc], method=cor.method)  )))^2)/(2*(sr^2))), 
                      
                      digits=20 )		
        #ffcor[startr:stopr, startc:stopc]<-temp
        temp<- 1-(temp1*temp2)
        temp[which(is.nan(temp))]<-1
        temp[which(is.na(temp))]<-1
        temp[which(is.infinite(temp))]<-1
        ffmat[startr:stopr, startc:stopc]<-temp
        rm(temp1); rm(temp2); rm(temp)
        gc()} 
      if(mint>maxt) {ffmat[startr:stopr, startc:stopc]<- 1}
    }
    gc()}
  
  ########
  # Call the similarity scoring function
  system.time(sapply(1:bl, RCsim))
  
  b<-Sys.time()
  
  ########
  # Report progress and timing
  cat('\n','\n' )
  cat(paste("RAMClust feature similarity matrix calculated and stored:", 
            round(difftime(b, a, units="mins"), digits=1), "minutes"))
  gc() 
  
  
  ########
  # extract lower diagonal of ffmat as vector
  blocksize<-mult*round(blocksize^2/n)
  nblocks<-floor(n/blocksize)
  remaind<-n-(nblocks*blocksize)
  
  ########
  # create vector for storing dissimilarities
  ramclustObj<-vector(mode="integer", length=vlength)
  
  ########
  # fill vector with dissimilarities
  for(k in 0:(nblocks)){
    startc<-1+(k*blocksize)
    if ((k+1)*blocksize > n) {
      stopc<-n} else {
        stopc<-(k+1)*blocksize}
    temp<-ffmat[startc:nrow(ffmat),startc:stopc]
    temp<-temp[which(row(temp)-col(temp)>0)]
    if(exists("startv")==FALSE) startv<-1
    stopv<-startv+length(temp)-1
    ramclustObj[startv:stopv]<-temp
    gc()
    startv<-stopv+1
    rm(temp)
    gc()
  }    
  rm(startv)
  gc()
  
  ########
  # convert vector to distance formatted object
  ramclustObj<-structure(ramclustObj, Size=(n), Diag=FALSE, Upper=FALSE, method="RAMClustR", Labels=featnames, class="dist")
  gc()
  
  c<-Sys.time()    
  cat('\n', '\n')
  cat(paste("RAMClust distances converted to distance object:", 
            round(difftime(c, b, units="mins"), digits=1), "minutes"))
  
  ########
  # cleanup
  delete.ff(ffmat)
  rm(ffmat)
  gc()
  
  
  ########
  # cluster using fastcluster package,
  system.time(ramclustObj<-hclust(ramclustObj, method=linkage))
  gc()
  d<-Sys.time()    
  cat('\n', '\n')    
  cat(paste("fastcluster based clustering complete:", 
            round(difftime(d, c, units="mins"), digits=1), "minutes"))
  if(minModuleSize==1) {
    clus<-cutreeDynamicTree(ramclustObj, maxTreeHeight=hmax, deepSplit=deepSplit, minModuleSize=2)
    sing<-which(clus==0)
    clus[sing]<-max(clus)+1:length(sing)
  }
  if(minModuleSize>1) {
    clus<-cutreeDynamicTree(ramclustObj, maxTreeHeight=hmax, deepSplit=deepSplit, minModuleSize=minModuleSize)
  }
  gc()
  
  ########
  # build results into ramclustObj
  ramclustObj$featclus<-clus
  ramclustObj$frt<-times
  ramclustObj$fmz<-mzs
  ramclustObj$xcmsOrd<-xcmsOrd
  msint<-rep(0, length(ramclustObj$fmz))
  for(i in 1:ncol(data1)){
    msint[i]<-weighted.mean(data1[,i], data1[,i])
  }
  ramclustObj$msint<-msint
  
  if(mslev==2) {
    msmsint<-rep(0, length(ramclustObj$fmz))	
    for(i in 1:ncol(data1)){	
      msmsint[i]<-weighted.mean(data2[,i], data2[,i])
    }
    ramclustObj$msmsint<-msmsint
  }
  
  clrt<-aggregate(ramclustObj$frt, by=list(ramclustObj$featclus), FUN="mean")
  ramclustObj$clrt<-clrt[which(clrt[,1]!=0),2]
  clrtsd<-aggregate(ramclustObj$frt, by=list(ramclustObj$featclus), FUN="sd")
  ramclustObj$clrtsd<-clrtsd[which(clrtsd[,1]!=0),2]
  
  ramclustObj$nfeat<-as.vector(table(ramclustObj$featclus)[2:max(ramclustObj$featclus)])
  ramclustObj$nsing<-length(which(ramclustObj$featclus==0))
  
  e<-Sys.time() 
  cat('\n', '\n')
  cat(paste("dynamicTreeCut based pruning complete:", 
            round(difftime(e, d, units="mins"), digits=1), "minutes"))
  
  f<-Sys.time()
  cat('\n', '\n')
  cat(paste("RAMClust has condensed", n, "features into",  max(clus), "spectra in", round(difftime(f, a, 
                                                                                                   
                                                                                                   units="mins"), digits=1), "minutes", '\n'))
  
  ramclustObj$ExpDes<-ExpDes
  strl<-nchar(max(ramclustObj$featclus)) - 1
  ramclustObj$cmpd<-paste("C", formatC(1:length(ramclustObj$clrt), digits = strl, flag = 0 ) , sep="")
  # cat(ramclustObj$cmpd[1:10], '\n')
  ramclustObj$ann<-ramclustObj$cmpd
  ramclustObj$annconf<-rep("", length(ramclustObj$clrt))
  ramclustObj$annnotes<-rep("", length(ramclustObj$clrt))
  ramclustObj$MSdata<-data1
  if(mslev==2) ramclustObj$MSMSdata<-data2  
  
  ########
  # collapse feature dataset into spectrum dataset
  if(collapse=="TRUE") {
    cat('\n', '\n', "... collapsing features into spectra")
    wts<-colSums(data1[])
    ramclustObj$SpecAbund<-matrix(nrow=nrow(data1), ncol=max(clus))
    for (ro in 1:nrow(ramclustObj$SpecAbund)) { 
      for (co in 1:ncol(ramclustObj$SpecAbund)) {
        ramclustObj$SpecAbund[ro,co]<- weighted.mean(data1[ro,which(ramclustObj$featclus==co)], wts[which(ramclustObj$featclus==co)])
      }
    }
    dimnames(ramclustObj$SpecAbund)[[2]]<-ramclustObj$cmpd
    if(!usePheno | is.null(xcmsObj)) {dimnames(ramclustObj$SpecAbund)[[1]]<-dimnames(ramclustObj$MSdata)[[1]]} 
    if(usePheno & !is.null(xcmsObj)) {dimnames(ramclustObj$SpecAbund)[[1]]<-as.vector(xcmsObj@phenoData[,1])[msfiles]}
    g<-Sys.time()
    cat('\n', '\n')
    cat(paste("RAMClustR has collapsed feature quantities
             into spectral quantities:", round(difftime(g, f, units="mins"), digits=1), "minutes", '\n'))
  }
  
  ########
  # further aggregate by sample names for 'SpecAbundAve' dataset
  
  rm(data1)
  rm(data2)
  if(!is.null(ramclustObj$SpecAbund)) {
    if(length(dimnames(ramclustObj$SpecAbund)[[1]])> length(unique(dimnames(ramclustObj$SpecAbund)[[1]]))) {
      ramclustObj$SpecAbundAve<-aggregate(ramclustObj$SpecAbund[,1:ncol(ramclustObj$SpecAbund)], 
                                 by=list(dimnames(ramclustObj$SpecAbund)[[1]]), 
                                 FUN="mean", simplify=TRUE)
      dimnames(ramclustObj$SpecAbundAve)[[1]]<-ramclustObj$SpecAbundAve[,1]
      ramclustObj$SpecAbundAve<-as.matrix(ramclustObj$SpecAbundAve[,2:ncol(ramclustObj$SpecAbundAve)])
      dimnames(ramclustObj$SpecAbundAve)[[2]]<-dimnames(ramclustObj$SpecAbund)[[2]]
      gc()
    }
  }
  gc()
  
  ########
  # write msp formatted spectra
  if(mspout==TRUE){ 
    cat(paste("writing msp formatted spectra...", '\n'))
    dir.create("spectra")
    libName<-paste("spectra/", ExpDes[[1]]["Experiment", 1], ".mspLib", sep="")
    file.create(file=libName)
    for (m in 1:as.numeric(mslev)){
      for (j in 1:max(ramclustObj$featclus)) {
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
    cat(paste('\n', "msp file complete", '\n')) 
  }  
  return(ramclustObj)
}
