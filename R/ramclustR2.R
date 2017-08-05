#' ramclustR
#'
#' Main clustering function 
#'
#' This is the Details section
#'
#' @param filename character Filename of the nmrML to check
#' @param ms MS1 intensities =MSdata, 
#' @param idmsms =ms,  
#' @param idMSMStag character e.g. "02.cdf"
#' @param featdelim character e.g. ="_"
#' @param timepos numeric 2
#' @param st numeric no clue e.g. = 5, 
#' @param sr numeric also no clue yet = 5, 
#' @param maxt numeric again no clue =20, 
#' @param deepSplit boolean e.g. =FALSE, 
#' @param blocksize integer number of features (scans?) processed in one block  =1000,
#' @param mult numeric =10
#'
#' @return A vector with the numeric values of the processed data
#' @author Corey Broeckling
#' @export

ramclustR2<- function(  xcmsObj=NULL,
                        ms=NULL, 
                        idmsms=NULL,
                        taglocation="filepaths",
                        MStag=NULL,
                        idMSMStag=NULL, 
                        featdelim="_", 
                        timepos=2, 
                        pkshape=FALSE,
                        sr=NULL, 
                        ss=NULL,
                        st=NULL, 
                        maxt=NULL, 
                        deepSplit=FALSE, 
                        blocksize=2000,
                        mult=5,
                        hmax=NULL,
                        sampNameCol=NULL,
                        collapse=TRUE,
                        usePheno=TRUE,
                        mspout=TRUE, 
                        mslev=1,
                        ExpDes=NULL,
                        normalize="TIC",
                        minModuleSize=2,
                        linkage="average",
                        mzdec=4,
                        #ffEIC=NULL,
                        #ffR=NULL,
                        #ffS=NULL,
                        #ffT=NULL,
                        saveFF=TRUE,
                        plots=TRUE,
                        cor.use="everything",
                        method="spearman",
                        min.scans=NULL) {
  
  require(xcms, quietly=TRUE)
  require(ff, quietly=TRUE)
  require(fastcluster, quietly=TRUE)
  require(dynamicTreeCut, quietly=TRUE)
  require(fastmatch, quietly=TRUE)
  if(plots) {require(ape, quietly=TRUE)}
  
  if(is.null(xcmsObj) & is.null(ms))  {
    stop("you must select either 
          1: an MS dataset with features as columns 
             (__)one column may contain sample names, defined by sampNameCol) 
          2: an xcmsObj. If you choose an xcms object, set taglocation: 'filepaths' by default
            and both MStag and idMSMStag")}
  
  if(!is.null(xcmsObj) & mslev==2 & any(is.null(MStag), is.null(idMSMStag), is.null(taglocation)))
  {stop("you must specify the the MStag, idMSMStag, and the taglocations")}
  
  if(is.null(ExpDes) & mspout==TRUE){
    cat("please enter experiment description (see popup window), or set 'mspout=FALSE'")
    ExpDes<-defineExperiment()
  }
  
  if(is.null(ExpDes) & mspout==FALSE){
    warning("using undefined instrumental settings")
    ExpDes<-structure(list(design = structure(list(ExpVars = c("", "", "", 
                                                               "", "", "-", "sn", "diet", "time", "tissue", "pos", "", "", "", 
                                                               "", ""), VarDesc = c("experiment name, no spaces", "lab and or user name", 
                                                                                    "species name (binomial latin name)", "sample type", "newLC newGC", 
                                                                                    "factor delimitor in your sample names", "Assign a name for your factors", 
                                                                                    "", "", "", "", "", "", "", "", "")), .Names = c("ExpVars", "VarDesc"
                                                                                    ), row.names = c("Experiment", "Contributor", "Species", "Sample", 
                                                                                                     "platform", "delim", "fact1name", "fact2name", "fact3name", "fact4name", 
                                                                                                     "fact5name", "fact6name", "fact7name", "fact8name", "fact9name", 
                                                                                                     "fact10name"), class = "data.frame"), instrument = structure(list(
                                                                                                       c.chrominst....Waters.UPLC..ACN.Gradient...msinst....Waters.Xevo.G2.QTOF... = structure(c(15L, 
                                                                                                                                                                                                 16L, 14L, 13L, 8L, 3L, 7L, 2L, 12L, 1L, 11L, 10L, 4L, 9L, 
                                                                                                                                                                                                 6L, 5L), .Label = c("0.05", "15-30", "2", "2200", "30", "50-2000", 
                                                                                                                                                                                                                     "6", "ACN, 0.1% formic acid", "Ar", "ESI", "P", "TOF", "Water", 
                                                                                                                                                                                                                     "Waters, 1 x 100 mm, 1.7 uM", "Waters UPLC: ACN Gradient", 
                                                                                                                                                                                                                     "Waters Xevo G2 QTOF"), class = "factor")), .Names = "c.chrominst....Waters.UPLC..ACN.Gradient...msinst....Waters.Xevo.G2.QTOF...", row.names = c("chrominst", 
                                                                                                                                                                                                                                                                                                                                                                       "msinst", "column", "solvA", "solvB", "MSlevs", "CE1", "CE2", 
                                                                                                                                                                                                                                                                                                                                                                       "mstype", "mzdifftof", "msmode", "ionization", "ESIvoltage", 
                                                                                                                                                                                                                                                                                                                                                                       "colgas", "msscanrange", "conevolt"), class = "data.frame")), .Names = c("design", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                "instrument"))
    
    
    
  }
  
  dir.create("ramclustR_out")
  if( normalize!="none"  & normalize!="TIC" & normalize!="quantile") {
    stop("please selected either 'none', 'TIC', or 'quantile' for the normalize setting")}
  
  a<-Sys.time()   
  
  if(is.null(hmax)) {hmax<-0.3}
  ##in using non-xcms data as input
  ##remove MSdata sets and save data matrix alone
  if(!is.null(ms)){
    if(is.null(st)) stop("please specify st: 
      a recommended starting point is half the value of 
      your average chromatographic peak width at half max (seconds)")
    if(is.null(sr)) sr<-0.5
    pkshape<-FALSE
    if(is.null(ss)) ss<-0.5
    if(is.null(maxt)) maxt<-60
    MSdata<-read.csv(ms, header=TRUE, check.names=FALSE)
    if(!is.null(idmsms)){
      MSMSdata<-read.csv(idmsms, header=TRUE, check.names=FALSE)}
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
  
  ##if xcms object is selected instead of an R dataframe/matrix
  if(!is.null(xcmsObj)){
    if(!class(xcmsObj)=="xcmsSet")
    {stop("xcmsObj must reference an object generated by XCMS, of class 'xcmsSet'")}
    
    if(is.null(st)) st<-round(median(xcmsObj@peaks[,"rtmax"]-xcmsObj@peaks[,"rtmin"])/2, digits=2)
    if(is.null(sr)) sr<-0.5
    if(is.null(maxt)) maxt<-st*20
    
    sampnames<-row.names(xcmsObj@phenoData)
    data12<-groupval(xcmsObj, value="into")
    
    if(taglocation=="filepaths" & !is.null(MStag) & mslev!=1) 
    { msfiles<-grep(MStag, xcmsObj@filepaths, ignore.case=TRUE)
    msmsfiles<-grep(idMSMStag, xcmsObj@filepaths, ignore.case=TRUE)
    if(length(intersect(msfiles, msmsfiles)>0)) 
    {stop("your MS and idMSMStag values do not generate unique file lists")}
    if(length(msfiles)!=length(msmsfiles)) 
    {stop("the number of MS files must equal the number of MSMS files")}
    if((length(msfiles)<1) | (length(msmsfiles)<1)) 
    {stop("your MStag and/or idMSMStag values return no files with matching strings")}
    data1<-t(data12[,msfiles])
    row.names(data1)<-sampnames[msfiles]
    data2<-t(data12[,msmsfiles])
    row.names(data2)<-sampnames[msmsfiles]  ##this may need to be changed to dimnames..
    times<-round(xcmsObj@groups[,"rtmed"], digits=3)
    mzs<-round(xcmsObj@groups[,"mzmed"], digits=4)
    } else {
      data1<-t(data12)
      msfiles<-1:nrow(data1)
      data2<-t(data12)
      times<-round(xcmsObj@groups[,"rtmed"], digits=3)
      mzs<-round(xcmsObj@groups[,"mzmed"], digits=4)      
    }
  }
  
  if((nrow(data1)/blocksize)-round(nrow(data1)/blocksize)<0.1) {blocksize=round(0.9*blocksize, digits=1)}
  ##replace na, inf, 0, and NaN with jittered min dataset value
  rpl1<-unique(c(which(is.na(data1)), which(is.nan(data1)), which(is.infinite(data1)), which(data1==0)))
  rpl2<-unique(c(which(is.na(data2)), which(is.nan(data2)), which(is.infinite(data2)), which(data2==0)))
  if(length(rpl1)>0) {data1[rpl1]<-jitter(rep(min(data1, na.rm=TRUE), length(rpl1) ), amount=min(data1/100, na.rm=TRUE))}
  if(length(rpl2)>0) {data2[rpl2]<-jitter(rep(min(data2, na.rm=TRUE), length(rpl2) ), amount=min(data2/100, na.rm=TRUE))}
  
  
  ##Optional normalization of data, either Total ion signal or quantile
  
  if(normalize=="TIC") {
    data1<-(data1/rowSums(data1))*mean(rowSums(data1), na.rm=TRUE)
    data2<-(data2/rowSums(data2))*mean(rowSums(data2), na.rm=TRUE)
  }
  if(normalize=="quantile") {
    library(preprocessCore)
    data1<-t(preprocessCore::normalize.quantiles(t(data1)))
    data2<-t(preprocessCore::normalize.quantiles(t(data2)))	
  }
  
  if(any(data1<0)) {
    warning("replacing negative MS intensity values with 0")
    data1[which(data1<0)]<-0
  }
  if(any(data2<0)) {
    warning("replacing negative MSMS intensity values with 0")
    data1[which(data2<0)]<-0
  }
  
  
  ##retention times and mzs vectors
  
  ##sort rt vector and data by retention time
  rtOrd<-order(times)
  data1<-data1[,rtOrd]
  data2<-data2[,rtOrd]
  mzs<-mzs[rtOrd]
  times<-times[rtOrd]
  
  ##generate feature names 
  featnames<-paste(mzs, "_", times, sep="")
  dimnames(data1)[[2]]<-featnames
  dimnames(data2)[[2]]<-featnames
  
  ##establish some constants for downstream processing
  n<-ncol(data1)
  vlength<-(n*(n-1))/2
  nblocks<-floor(n/blocksize)
  
  ##create three or four empty matrices, one each for the correlation (ffR) matrix, 
  ##the rt (ffT) matrix,  the peak shape matrix(ffS), and the product matrix (ffP)
  ##if the files already exist, simply load them instead (except product matrix)
  
  if(!file.exists("ramclustR_out/ffr.Rdata")) {ffr<-ff(vmode="double", dim=c(n, n), init=0)} else {
    cat("loading existing ffr")
    cat('/n')
    suppressWarnings(ffload(file="ramclustR_out/ffr", overwrite=TRUE))
    if(dim(xcmsObj@groups)[1]!=nrow(ffr)) {stop("dimensions of feature groups not equal to dimensions of ffr")}
  }
  if(!file.exists("ramclustR_out/fft.Rdata")) {fft<-ff(vmode="double", dim=c(n, n), init=0)} else {
    cat("loading existing fft")
    cat('/n')
    suppressWarnings(ffload(file="ramclustR_out/fft", overwrite=TRUE))
    if(dim(xcmsObj@groups)[1]!=nrow(fft)) {stop("dimensions of feature groups not equal to dimensions of fft")} 
  }
  if(pkshape & !file.exists("ramclustR_out/ffs.Rdata")) {ffs<-ff(vmode="double", dim=c(n, n), init=0)} else {
    if(file.exists("ramclustR_out/ffs.Rdata")) {
      cat("loading existing ffs")
      cat('/n')
      suppressWarnings(ffload(file="ramclustR_out/ffs", overwrite=TRUE))
      if(dim(xcmsObj@groups)[1]!=nrow(ffs)) {stop("dimensions of feature groups not equal to dimensions of ffs")}
    }
  }
  ffp<-ff(vmode="double", dim=c(n, n), init=0)
  gc()
  #Sys.sleep((n^2)/10000000)
  #gc()
  
  ##make list of all row and column blocks to evaluate
  eval1<-expand.grid(0:nblocks, 0:nblocks)
  names(eval1)<-c("j", "k") #j for cols, k for rows
  eval1<-eval1[which(eval1[,"j"]<=eval1[,"k"]),] #upper triangle only
  bl<-nrow(eval1)
  row.names(eval1)<-c(1:bl)
  
  
  ##PEAK SHAPE START
  if(pkshape) {
    if(file.exists("ramclustR_out/ffeic.Rdata")) {
      cat("loading existing ffeic")
      cat('/n')
      suppressWarnings(ffload(file="ramclustR_out/ffeic", overwrite=TRUE))
      if(dim(xcmsObj@groups)[1]!=nrow(ffeic)) {stop("dimensions of feature groups not equal to dimensions of ffEIC, delete existing ff objects from output directory")}
    } else {
      cat('\n', paste("Generating EICs for peak shape similarity scoring"), '\n')
      pks<-xcmsObj@peaks
      pks<-pks[-xcmsObj@filled,]
      grpidx<-unlist(xcmsObj@groupidx)
      grpidfun<-function(x) {
        rep(x, length(xcmsObj@groupidx[[x]]))
      }
      grpid<-unlist(sapply(1:length(xcmsObj@groupidx), FUN=grpidfun))
      
      matfun<- function (x) {
        fmatch(x, grpidx)
      }
      
      grp<-grpid[sapply(1:nrow(pks), FUN=matfun)]
      pks<-data.frame(pks, grp)
      #  pks[which(grp==1647),]
      
      rem<-which(is.na(grp))
      if(length(rem)>0) {
        pks<-pks[-rem,]
        grp<-grp[-rem]
      }
      
      
      rmz<-as.matrix(pks[,c("mzmin", "mzmax")])
      rrt<-as.matrix(pks[,c("rtmin", "rtmax")])
      
      rtrange<-c(min(xcmsObj@rt$raw[[1]]), max(xcmsObj@rt$raw[[1]]))
      rtlen<-length(xcmsObj@rt$raw[[1]])
      rtstep<-(rtrange[2]-rtrange[1])/rtlen
      
      ind<-round((pks[, "rt"]-rtrange[1])/rtstep)
      samp<-pks[,"sample"]
      rtcorfun<-function(x) {
        rtadjn<-	(xcmsObj@rt$raw[[samp[x]]][ind[x]])- 
          (xcmsObj@rt$corrected[[samp[x]]][ind[x]])
        return(rtadjn)
      } #  0.08 seconds for 10,000 features
      
      ## rtadj IS TO BE USED BELOW: The EIC values will be merged to create a collective EIC chromatogram matrix - with the rt slots adjusted using this value for each peak!
      rtadj<-unlist(sapply(1:length(grp), FUN=rtcorfun))
      rtadj[which(is.na(rtadj))]<-0
      
      eicrtrange<-c(floor(min(pks[,"rtmin"])), ceiling(max(pks[,"rtmax"])))
      npts<-5*eicrtrange[2]-eicrtrange[1]
      eictimes<-seq(from=eicrtrange[1], to=eicrtrange[2], by=((eicrtrange[2]-eicrtrange[1])/(npts-1)))
      fford<-order(rtOrd)
      
      ffeic <- ff(vmode = "double", dim = c(dim(xcmsObj@groups)[1], npts), init = 0)
      for(j in 1:max(samp)) {
        do<-which(samp==j)
        #if(any(ls()=="xseic")){t1<-xseic[[j]]; cat('\n', paste(j)) } else {
        t1<-getEIC(xcmsObj, as.matrix(pks[do,c("mzmin", "mzmax")]), as.matrix(pks[do,c("rtmin", "rtmax")]), sampleidx=j) #
        #}
        rtadjdo<-rtadj[do]
        grpdo<-grp[do]  #grpdo<-fford[grp[do]]
        for(i in 1:length(do)) {
          art<-t1@eic[[1]][[i]][,"rt"]+rtadjdo[i]
          aint<-t1@eic[[1]][[i]][,"intensity"]+rtadjdo[i]
          aint<-aint-min(aint)
          from<-which(eictimes-min(art)> 0 )[1]-1
          if(!is.na(from)){
            #cat("from =", from, '\n')
            #if(from<1){cat("from ", j," ", i, "\n")}
            to<- which(eictimes-max(art)> 0 )[1]
            if(to>ncol(ffeic)){cat("to ", j," ", i, "\n")}
            subtime<-eictimes[from:to]
            is<-ffeic[fford[grpdo[i]],from:to] 
            is[which(is.na(is))]<-0
            spl<-(round(spline(art, aint, n=length(subtime), method="fmm")$y))
            spl[which(is.na(spl))]<-0
            spl[which(spl<0)]<-0
            #spl<-spl-min(spl)
            ffeic[fford[grpdo[i]] ,from:to]<- is+spl
          }	
        }
        ffeic[which(ffeic[]<0)]<-0
      }
      gc()
      
      
    }
    
  } ##PEAK SHAPE END
  
  cat('\n', paste("calculating ramclustR similarity: nblocks = ", bl))
  cat('\n', "finished:")  
  
  makeFF<-function(block)  {
    cat(block,' ')
    j<-eval1[block,"j"]  #columns
    k<-eval1[block,"k"]  #rows
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
        tt<-(abs(outer(times[startr:stopr], times[startc:stopc], FUN="-")))
        fft[startr:stopr, startc:stopc]<-tt
        tr<-pmax(cor(data1[,startr:stopr], data1[,startc:stopc]),
                 cor(data1[,startr:stopr], data2[,startc:stopc]),
                 cor(data2[,startr:stopr], data2[,startc:stopc]))
        tr[which(tr<0)]<-0
        ffr[startr:stopr, startc:stopc]<-tr        
        if(pkshape) {
          tempeicr<-t(ffeic[startr:stopr,])
          tempeicc<-t(ffeic[startc:stopc,])
          if(length(intersect(which(rowSums(tempeicr, na.rm=TRUE)>0), which(rowSums(tempeicc, na.rm=TRUE)>0)) )>2){
            #dor<-which(rowSums(tempeicr, na.rm=TRUE)>0)
            #doc<-which(rowSums(tempeicc, na.rm=TRUE)>0)
            ts <- cor(tempeicr, tempeicc, use=cor.use, method=method)
            #ts <- cor(tempeicr[dor,], tempeicc[doc,], use=cor.use)
            ts[which(ts<0)]<-0
            ts[which(is.na(ts))]<-0
            if(!is.null(min.scans)) {
              out<-matrix(ncol=ncol(tempeicc), nrow=ncol(tempeicr))
              #system.time(
              for(i in 1:ncol(tempeicc)) {
                out[i:ncol(tempeicr),i]<-mapply(FUN=olf, i, c(i:ncol(tempeicr)))
              }
              #)
              out[which(out<min.scans)]<-0
              out[which(is.na(out))]<-0
              out[which(out>=min.scans)]<-1
              out[upper.tri(out)] <- t(out)[upper.tri(out)] 
              ts<-ts*out
            }
            
          } else {ts<-matrix(ncol=ncol(tr), nrow=nrow(tr)); ts[]<-0}
          
          ffs[startr:stopr, startc:stopc]<-ts
        }
        #c<-Sys.time()
        # b-a
        #c-b
        #if(dim(cr)!=dim(temp1)) {stop("dimensions incompatible")}
      }
      gc()} 
  }
  
  if(!file.exists("ramclustR_out/ffr.Rdata") & !file.exists("ramclustR_out/fft.Rdata")) {system.time(sapply(1:bl, makeFF))}
  #save ff objects, if not already used
  if(!file.exists("ramclustR_out/ffeic.Rdata") & saveFF) {
    ffsave(ffeic, file="ramclustR_out/ffeic")
  }
  
  RCsim<-function(block)  {
    cat(block,' ')
    j<-eval1[block,"j"]  #columns
    k<-eval1[block,"k"]  #rows
    startc<-min((1+(j*blocksize)), n)
    if ((j+1)*blocksize > n) {
      stopc<-n} else {
        stopc<-(j+1)*blocksize}
    startr<-min((1+(k*blocksize)), n)
    if ((k+1)*blocksize > n) {
      stopr<-n} else {
        stopr<-(k+1)*blocksize}
    if(startc<=startr) { 
      tt<-fft[startr:stopr, startc:stopc]
      mint<-min(tt)
      tr<-ffr[startr:stopr, startc:stopc]
      if(mint<=maxt) {
        temp1<-round(exp(-(( tt )^2)/(2*(st^2))), digits=8 )
        if(is.null(sr)) {temp2<-tr} else {temp2<-round(exp(-((1-tr)^2)/(2*(sr^2))), digits=8 )}
        #temp2<-round(exp(-((1-(tr))^2)/(2*(sr^2))), digits = 8 )  
        
        #ffcor[startr:stopr, startc:stopc]<-temp
        if(pkshape) {
          ts<-ffs[startr:stopr, startc:stopc]
          if(is.null(ss)) {temp3<-ts} else {temp3<-round(exp(-((1-ts)^2)/(2*(ss^2))), digits=8 )}
          cat("\n")
        }
        
        if(pkshape) {temp<-1-((temp1*temp2*temp3)^(1/3))} else {temp<- 1-((temp1*temp2)^(1/2))}
        
        temp[which(is.nan(temp))]<-1
        temp[which(is.na(temp))]<-1
        temp[which(is.infinite(temp))]<-1
        ffp[startr:stopr, startc:stopc]<-temp
        rm(temp1); rm(temp2); rm(temp); rm(temp3)
        gc()} 
      if(mint>maxt) {ffp[startr:stopr, startc:stopc]<- 1}
    }
    gc()} 
  # ffmat[995:1002,995:1002]
  
  
  ##Call the similarity scoring function
  system.time(sapply(1:bl, RCsim))
  #RCsim(bl=1:bl)
  
  b<-Sys.time()
  
  if(!file.exists("ramclustR_out/ffr.Rdata") & saveFF) {
    ffsave(ffr, file="ramclustR_out/ffr")
  }
  if(!file.exists("ramclustR_out/ffs.Rdata") &  saveFF) {
    ffsave(ffs, file="ramclustR_out/ffs")
  }
  if(!file.exists("ramclustR_out/fft.Rdata") & saveFF) {
    ffsave(fft, file="ramclustR_out/fft")
  }
  cat('\n','\n' )
  cat(paste("RAMClust feature similarity matrix calculated and stored:", 
            round(difftime(b, a, units="mins"), digits=1), "minutes"))
  
  
  gc() 
  
  
  ##extract lower diagonal of ffmat as vector
  blocksize<-mult*round(blocksize^2/n)
  nblocks<-floor(n/blocksize)
  remaind<-n-(nblocks*blocksize)
  
  ##create vector for storing dissimilarities
  RCd<-vector(mode="integer", length=vlength)
  
  for(k in 0:(nblocks)){
    startc<-1+(k*blocksize)
    if ((k+1)*blocksize > n) {
      stopc<-n} else {
        stopc<-(k+1)*blocksize}
    temp<-ffp[startc:nrow(ffp),startc:stopc]
    temp<-temp[which(row(temp)-col(temp)>0)]
    if(exists("startv")==FALSE) startv<-1
    stopv<-startv+length(temp)-1
    RCd[startv:stopv]<-temp
    gc()
    startv<-stopv+1
    rm(temp)
    gc()
  }    
  rm(startv)
  gc()
  
  
  ##convert vector to distance formatted object
  RCd<-structure(RCd, Size=(n), Diag=FALSE, Upper=FALSE, method="RAMClustR", Labels=featnames, class="dist")
  gc()
  
  c<-Sys.time()    
  cat('\n', '\n')
  cat(paste("RAMClust distances converted to distance object:", 
            round(difftime(c, b, units="mins"), digits=1), "minutes"))
  
  
  gc()
  
  
  ##cluster using fastcluster package,
  system.time(RC<-hclust(RCd, method=linkage))
  gc()
  d<-Sys.time()    
  cat('\n', '\n')    
  cat(paste("fastcluster based clustering complete:", 
            round(difftime(d, c, units="mins"), digits=1), "minutes"))
  clus<-as.vector(cutreeHybrid(dendro=RC, distM=as.matrix(RCd), minClusterSize = max(minModuleSize, 2), deepSplit = deepSplit, cutHeight=hmax)$labels)
  #clus<-cutreeDynamicTree(RC, maxTreeHeight=hmax, deepSplit=deepSplit, minModuleSize=max(minModuleSize, 2))
  if(minModuleSize==1) {
    sing<-which(clus==0)
    clus[sing]<-max(clus)+1:length(sing)
  }
  gc()
  
  
  RC$featclus<-clus
  RC$frt<-times
  RC$fmz<-mzs
  RC$rtOrd<-rtOrd
  msint<-rep(0, length(RC$fmz))
  for(i in 1:ncol(data1)){
    msint[i]<-weighted.mean(data1[,i], data1[,i])
  }
  RC$msint<-msint
  
  if(mslev==2) {
    msmsint<-rep(0, length(RC$fmz))	
    for(i in 1:ncol(data1)){	
      msmsint[i]<-weighted.mean(data2[,i], data2[,i])
    }
    RC$msmsint<-msmsint
  }
  clrt<-aggregate(RC$frt, by=list(RC$featclus), FUN="mean")
  RC$clrt<-clrt[which(clrt[,1]!=0),2]
  clrtsd<-aggregate(RC$frt, by=list(RC$featclus), FUN="sd")
  RC$clrtsd<-clrtsd[which(clrtsd[,1]!=0),2]
  
  RC$nfeat<-as.vector(table(RC$featclus)[2:max(RC$featclus)])
  RC$nsing<-length(which(RC$featclus==0))
  
  e<-Sys.time() 
  cat('\n', '\n')
  cat(paste("dynamicTreeCut based pruning complete:", 
            round(difftime(e, d, units="mins"), digits=1), "minutes"))
  
  f<-Sys.time()
  cat('\n', '\n')
  cat(paste("RAMClust has condensed", n, "features into",  max(clus), "spectra in", round(difftime(f, a, units="mins"), digits=1), "minutes", '\n'))
  
  RC$ExpDes<-ExpDes
  RC$cmpd<-paste("C", 1:length(RC$clrt), sep="")
  RC$ann<-RC$cmpd
  RC$annconf<-rep("", length(RC$clrt))
  RC$annnotes<-rep("", length(RC$clrt))
  RC$MSdata<-data1
  if(mslev==2) RC$MSMSdata<-data2  
  
  if(collapse=="TRUE") {
    cat('\n', '\n', "... collapsing features into spectra")
    wts<-colSums(data1[])
    RC$SpecAbund<-matrix(nrow=nrow(data1), ncol=max(clus))
    for (ro in 1:nrow(RC$SpecAbund)) { 
      for (co in 1:ncol(RC$SpecAbund)) {
        RC$SpecAbund[ro,co]<- weighted.mean(data1[ro,which(RC$featclus==co)], wts[which(RC$featclus==co)])
      }
    }
    dimnames(RC$SpecAbund)[[2]]<-paste("C", 1:ncol(RC$SpecAbund), sep="")
    if(!usePheno) {dimnames(RC$SpecAbund)[[1]]<-dimnames(RC$MSdata)[[1]]} 
    if(usePheno) {dimnames(RC$SpecAbund)[[1]]<-xcmsObj@phenoData[,1][msfiles]}
    g<-Sys.time()
    cat('\n', '\n')
    cat(paste("RAMClustR has collapsed feature quantities
             into spectral quantities:", round(difftime(g, f, units="mins"), digits=1), "minutes", '\n'))
  }
  
  rm(data1)
  rm(data2)
  if(!is.null(RC$SpecAbund)) {
    if(length(dimnames(RC$SpecAbund)[[1]])> length(unique(dimnames(RC$SpecAbund)[[1]]))) {
      RC$SpecAbundAve<-aggregate(RC$SpecAbund[,1:ncol(RC$SpecAbund)], 
                                 by=list(dimnames(RC$SpecAbund)[[1]]), 
                                 FUN="mean", simplify=TRUE)
      dimnames(RC$SpecAbundAve)[[1]]<-RC$SpecAbundAve[,1]
      RC$SpecAbundAve<-as.matrix(RC$SpecAbundAve[,2:ncol(RC$SpecAbundAve)])
      dimnames(RC$SpecAbundAve)[[2]]<-dimnames(RC$SpecAbund)[[2]]
      gc()
    }
  }
  gc()
  
  if(mspout==TRUE){ 
    cat(paste("writing msp formatted spectra...", '\n'))
    dir.create("spectra")
    libName<-paste("spectra/", ExpDes[[1]]["Experiment", 1], ".mspLib", sep="")
    file.create(file=libName)
    for (m in 1:as.numeric(mslev)){
      for (j in 1:max(RC$featclus)) {
        #print(paste(j,"_", sep=""))
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
        mz<-round(RC$fmz[sl][order(wm, decreasing=TRUE)], digits=mzdec)
        rt<-RC$frt[sl][order(wm, decreasing=TRUE)]
        wm<-round(wm[order(wm, decreasing=TRUE)])
        mrt<-mean(rt)
        npeaks<-length(mz)
        for (l in 1:length(mz)) {
          ion<- paste(mz[l], wm[l])
          if(l==1) {specdat<-ion} 
          if(l>1)  {specdat<-c(specdat, " ", ion)}
        }
        cat(
          paste("Name: C", j, sep=""), '\n',
          paste("SYNON: $:00in-source", sep=""), '\n',
          paste("SYNON: $:04", sep=""), '\n', 
          paste("SYNON: $:05", if(m==1) {ExpDes[[2]]["CE1", 1]} else {ExpDes$instrument["CE2", "InstVals"]}, sep=""), '\n',
          paste("SYNON: $:06", ExpDes[[2]]["mstype", 1], sep=""), '\n',           #mstype
          paste("SYNON: $:07", ExpDes[[2]]["msinst", 1], sep=""), '\n',           #msinst
          paste("SYNON: $:09", ExpDes[[2]]["chrominst", 1], sep=""), '\n',        #chrominst
          paste("SYNON: $:10", ExpDes[[2]]["ionization", 1], sep=""),  '\n',      #ionization method
          paste("SYNON: $:11", ExpDes[[2]]["msmode", 1], sep=""), '\n',           #msmode
          paste("SYNON: $:12", ExpDes[[2]]["colgas", 1], sep=""), '\n',           #collision gas
          paste("SYNON: $:14", ExpDes[[2]]["msscanrange", 1], sep=""), '\n',      #ms scanrange
          paste("SYNON: $:16", ExpDes[[2]]["conevolt", 1], sep=""), '\n',         #conevoltage
          paste("Comment: Rt=", round(mrt, digits=2), 
                "  Contributor=", ExpDes[[1]]["Contributor", 1], 
                "  Study=", ExpDes[[1]]["Experiment", 1], 
                sep=""), '\n',
          paste("Num Peaks:", npeaks), '\n',
          paste(specdat), '\n', '\n', sep="", file=libName, append= TRUE)
      }
    }
    cat(paste('\n', "msp file complete", '\n')) 
  } 
  if(plots) {
    
    filename<-paste(getwd(),"/ramclustR_out/ramculstRplots_", strftime(Sys.time(), "%Y%m%d%H%M%S"), ".pdf", sep="")
    pdf(file=filename, useDingbats=FALSE, width=10, height=10)
    th<-hist(RCd, breaks=1000)
    ax.fact<-ceiling(max(th$counts)/max(th$counts[which(th$breaks<0.2)]))
    plot(th, ylim=c(0, max(th$counts[which(th$breaks<0.2)])), xlab="RC dissimilarities", main=paste("dissimilarity histogram: y-axis scaled by", ax.fact, "fold"))
    points(th$mids, round(th$counts/ax.fact), type="h", col="red",lwd=2.5)
    
    
    RCp<-as.phylo(RC)
    RCp$tip.label<-paste(RC$featclus, round(RC$fmz, digits=2), round(RC$frt, digits=1), sep="_")
    col<-RC$featclus+1
    plot(RCp, type="fan", show.tip.label=TRUE, cex=0.05, tip.color=col, edge.width=0.25)
    for(x in 1:max(RC$featclus)) {
      tmpeic<-as.matrix(ffeic[which(RC$featclus==x),])
      pres<-which((colSums(tmpeic, na.rm=TRUE))>0)
      startt<-max((min(pres)-10),1)
      stopt<-min((max(pres+10)), length(eictimes))
      maxy<-max(tmpeic[,startt:stopt], na.rm=TRUE)
      expon<-0.2
      #cat(paste("made it ", x, startt, stopt))
      #cat(tmpeic[1,])
      #       plot(eictimes[startt:stopt], tmpeic[1,startt:stopt]/max(tmpeic[1,startt:stopt], na.rm=TRUE), type="l",  
      #            main=paste("cl", x, "  rt =", round(RC$clrt[x], digits=2)), xlab="rt", ylab="intensity", ylim=c(0,1.4))
      #       text((eictimes[startt:stopt][which.max(tmpeic[1,startt:stopt])]), 1.2, round(RC$fmz[which(RC$featclus==x)][1], digits=3), cex=0.5 )
      plot(eictimes[startt:stopt], (tmpeic[1,startt:stopt])^expon, type="l",  
           main=paste("cl", x, "  rt =", round(RC$clrt[x], digits=2)), xlab="rt", ylab="intensity", ylim=c(0, (maxy*1.4)^expon ))
      text((eictimes[startt:stopt][which.max(tmpeic[1,startt:stopt])]), (maxy*1.2)^expon, round(RC$fmz[which(RC$featclus==x)][1], digits=3), cex=0.5 )
      if(nrow(tmpeic)>1) {
        for(y in 2:nrow(tmpeic)) {
          points(eictimes[startt:stopt], (tmpeic[y,startt:stopt])^expon, type="l", col=y)
          text((eictimes[startt:stopt][which.max(tmpeic[y,startt:stopt])]), jitter( (maxy*1.2)^expon), round(RC$fmz[which(RC$featclus==x)][y], digits=3), cex=0.5 )
        }
      }
    }
    dev.off()
  }
  
  #cleanup
  delete.ff(ffr)
  rm(ffr)
  delete.ff(fft)
  rm(fft) 
  if(pkshape) {
    delete.ff(ffeic)
    rm(ffeic)
    delete.ff(ffs)
    rm(ffs)
  }
  return(RC)
}
