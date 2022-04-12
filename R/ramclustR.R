#' ramclustR
#'
#' Main clustering function for grouping features based on their analytical behavior.  
#'
#' @param xcmsObj xcmsObject: containing grouped feature data for clustering by ramclustR
#' @param ms filepath: optional csv input. Features as columns, rows as samples. Column header mz_rt
#' @param idmsms filepath: optional idMSMS / MSe csv data.  same dim and names as ms required
#' @param MStag character: character string in 'taglocation' to designat MS / MSe files e.g. "01.cdf"
#' @param idMSMStag character: character string in 'taglocation' to designat idMSMS / MSe files e.g. "02.cdf"
#' @param taglocation character: "filepaths" by default, "phenoData[,1]" is another option. refers to xcms slot
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
#' @param usePheno logical: transfer phenotype data from XCMS object to SpecAbund dataset?
#' @param mspout logical: write msp formatted spectra to file?
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param normalize character: either "none", "TIC", "quantile", or "batch.qc" normalization of feature intensities.  see batch.qc overview in details. 
#' @param qc.inj.range integer: how many injections around each injection are to be scanned for presence of QC samples when using batch.qc normalization?  A good rule of thumb is between 1 and 3 times the typical injection span between QC injections.  i.e. if you inject QC ever 7 samples, set this to between 7 and 21.  smaller values provide more local precision but make normalization sensitive to individual poor outliers (though these are first removed using the boxplot function outlier detection), while wider values provide less local precision in normalization but better stability to individual peak areas.
#' @param batch integer vector with length equal to number of injections in xset or csv file
#' @param order integer vector with length equal to number of injections in xset or csv file
#' @param qc logical vector with length equal to number of injections in xset or csv file.  
#' @param minModuleSize integer: how many features must be part of a cluster to be returned? default = 2
#' @param linkage character: heirarchical clustering linkage method - see ?hclust
#' @param mzdec integer: number of decimal places used in printing m/z values
#' @param cor.method character: which correlational method used to calculate 'r' - see ?cor
#' @param rt.only.low.n logical: default = TRUE  At low injection numbers, correlational relationships of peak intensities may be unreliable.  by defualt ramclustR will simply ignore the correlational r value and cluster on retention time alone.  if you wish to use correlation with at n < 5, set this value to FALSE.
#' @param fftempdir valid path: if there are file size limitations on the default ff package temp directory  - getOptions('fftempdir') - you can change the directory used as the fftempdir with this option.
#' @param replace.zeros logical: TRUE by default.  NA, NaN, and Inf values are replaced with zero, and zero values are sometimes returned from peak peaking.  When TRUE, zero values will be replaced with a small amount of noise, with noise level set based on the detected signal intensities for that feature. 
#' @details Main clustering function output - see citation for algorithm description or vignette('RAMClustR') for a walk through.  batch.qc. normalization requires input of three vectors (1) batch (2) order (3) qc.   This is a feature centric normalization approach which adjusts signal intensities first by comparing batch median intensity of each feature (one feature at a time) QC signal intensity to full dataset median to correct for systematic batch effects and then secondly to apply a local QC median vs global median sample correction to correct for run order effects.
#' @return   $featclus: integer vector of cluster membership for each feature
#' @return   $frt: feature retention time, in whatever units were fed in (xcms uses seconds, by default)
#' @return   $fmz: feature retention time, reported in number of decimal points selected in ramclustR function
#' @return   $xcmsOrd: the original XCMS (or csv) feature order for cross referencing, if need be
#' @return   $clrt: cluster retention time
#' @return   $clrtsd: retention time standard deviation of all the features that comprise that cluster
#' @return   $nfeat: number of features in the cluster
#' @return   $nsing: number of 'singletons' - that is the number of features which clustered with no other feature
#' @return   $ExpDes: the experimental design object used when running ramclustR.  List of two dataframes. 
#' @return   $cmpd: compound name.  C#### are assigned in order of output by dynamicTreeCut.  Compound with the most features is classified as C0001...
#' @return   $ann: annotation.  By default, annotation names are identical to 'cmpd' names.  This slot is a placeholder for when annotations are provided
#' @return   $MSdata:  the MSdataset provided by either xcms or csv input
#' @return   $MSMSdata: the (optional) MSe/idMSMS dataset provided be either xcms or csv input
#' @return   $SpecAbund: the cluster intensities after collapsing features to clusters
#' @return   $SpecAbundAve: the cluster intensities after averaging all samples with identical sample names
#' @return   - 'spectra' directory is created in the working directory.  In this directory a .msp is (optionally) created, which contains the spectra for all compounds in the dataset following clustering.  if MSe/idMSMS data are provided, they are listed width he same compound name as the MS spectrum, with the collision energy provided in the ExpDes object provided to distinguish low from high CE spectra. 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis boxplot hist "legend" "par" "plot" "points" "title"
#' @importFrom stats aggregate cor fitted lm loess median predict quantile sd weighted.mean
#' @importFrom utils edit read.csv read.delim read.delim2 write.csv packageVersion
#' @importFrom ff ff
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom e1071 skewness
#' @importFrom gplots heatmap.2
#' @importFrom pcaMethods pca
#' @importFrom jsonlite fromJSON
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom utils citation packageVersion
#' 
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
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
#' @examples
#' ## Choose input file with feature column names `mz_rt` (expected by default).
#' ## Column with sample name is expected to be first (by default).
#' ## These can be adjusted with the `featdelim` and `sampNameCol` parameters.
#' wd <- getwd()
#' filename <- system.file("extdata", "peaks.csv", package = "RAMClustR", mustWork = TRUE)
#' print(filename)
#' head(data.frame(read.csv(filename)), c(6L, 5L))
#' 
#' ## If the file contains features from MS1, assign those to the `ms` parameter.
#' ## If the file contains features from MS2, assign those to the `idmsms` parameter.
#' ## If you ran `xcms` for the feature detection, the assign the output to the `xcmsObj` parameter.
#' ## In this example we use a MS1 feature table stored in a `csv` file.
#' setwd(tempdir())
#' ramclustobj <- ramclustR(ms = filename, st = 5, maxt = 1, blocksize = 1000)
#' 
#' ## Investigate the deconvoluted features in the `spectra` folder in MSP format
#' ## or inspect the `ramclustobj` for feature retention times, annotations etc.
#' print(ramclustobj$ann)
#' print(ramclustobj$nfeat)
#' print(ramclustobj$SpecAbund[,1:6])
#' setwd(wd)

ramclustR  <- function(xcmsObj=NULL,
                       ms=NULL, 
                       idmsms=NULL,
                       taglocation="filepaths",
                       MStag=NULL,
                       idMSMStag=NULL, 
                       featdelim="_", 
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
                       qc.inj.range = 20,
                       order = NULL,
                       batch = NULL,
                       qc = NULL,
                       minModuleSize=2,
                       linkage="average",
                       mzdec=3,
                       cor.method="pearson",
                       rt.only.low.n = TRUE,
                       fftempdir = NULL,
                       replace.zeros = TRUE
) {
  
  ########
  # If experimental design is NULL: 
  if(is.null(ExpDes)) {
    ExpDes <- defineExperiment(force.skip = TRUE)
    warning('\n', "you failed to define your experimental descriptor using 'defineExperiment()'", '\n',
            "RAMClustR must now guess at what you are trying to do ", '\n',
            "and your exported spectra will be labelled incorrectly")
    if(!is.null(idmsms)) {
      ExpDes[[2]][which(row.names(ExpDes[[2]]) == "MSlevs"),1]<-2
    }
  }
  
  if(normalize == "batch.qc") {
    if(is.null(order) | is.null(batch) |is.null(qc)) {
      stop('to use batch.qc normalization you must provide vectors for batch, order (run order) and qc information as vectors.  see help ?ramclustR')
    }
  }
  
  if(!is.null(fftempdir)) {
    origffdir<-getOption("fftempdir")
    options("fftempdir" = fftempdir)
  }
  
  ########
  # define ms levels, used several times below
  mslev <- as.integer(as.numeric(as.character(ExpDes[[2]][which(row.names(ExpDes[[2]]) == "MSlevs"),1])))
  
  history <- paste("Raw mass spectrometry data were processed using an R based workflow for feature detection, retention time alignment, feature grouping, peak filling, feature clustering.")

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
  
  if( normalize!="none"  & normalize!="TIC" & normalize!="quantile" & normalize!="batch.qc") {
    stop("please selected either 'none', 'TIC', or 'quantile' for the normalize setting")}
  
  
  a<-Sys.time()   
  
  ########
  # set default parameters of not defined in function call
  
  if(is.null(hmax)) {hmax<-0.3}
  
  cat(paste("  organizing dataset", '\n'))
  
  ########
  # if csv input of MS data, do this: 
  if(!is.null(ms)){
    if(is.null(st)) stop("please specify st: 
      a recommended starting point is half the value of 
      your average chromatographic peak width at half max (seconds)")
    if(is.null(sr)) sr<-0.5
    if(is.null(maxt)) maxt<-60
    
    history <- paste(history,
                     " Feature data was input as .csv files with ramclustR parameter settings of ",
                     " st = ", st,
                     " sr = ", sr,
                     " and maxt = ", maxt, ".", sep = ""
    )
    
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
    if(any(is.na(times))) {
      stop("column(s) ", which(is.na(times)), " rt cannot be made numeric. ")
    }
    if(any(is.na(mzs))) {
      stop("column(s) ", which(is.na(times)), " mz cannot be made numeric. ")
    }
    rm(rtmz)
  }
  
  newXCMS <- FALSE
  ########
  # if xcms object as input, do this:
  if(!is.null(xcmsObj)){
    
    if (!requireNamespace("xcms", quietly = TRUE)) {
      stop("The use of this function requires package 'xcms'. Please ",
           "install with 'Biobase::install(\"xcms\")'")
    }
    
    OK <- FALSE
    
    if(inherits(xcmsObj, "xcmsSet")) {OK <- TRUE} 
    if(inherits(xcmsObj, "XCMSnExp")) {
      OK <- TRUE
      newXCMS <- TRUE
    }
    if(!OK) {stop("xcmsObj must reference an object generated by XCMS, of class 'xcmsSet'")}
    
    if(is.null(st) & !newXCMS) st<-round(median(xcmsObj@peaks[,"rtmax"]-xcmsObj@peaks[,"rtmin"])/2, digits=2)
    if(is.null(st) & newXCMS) st<-round(median(xcmsObj@msFeatureData$chromPeaks[,"rtmax"]-xcmsObj@msFeatureData$chromPeaks[,"rtmin"])/2, digits=2)
    if(is.null(sr)) sr<-0.5
    if(is.null(maxt)) maxt<-st*20
    
    history <- paste(history,
                     " XCMS",  paste0("(v.", packageVersion('xcms'), ")"), "was used for feature detection and retention time alighment. ",
                     "Processing was performed using R", paste0("(v.", citation()$author, " ",  citation()$year, ")."),
                     " Feature data was input as an xcms object with ramclustR parameter settings of ",
                     " st = ", st,
                     " sr = ", sr,
                     " and maxt = ", maxt, ".", sep = ""
    )
    
    if(!newXCMS) {sampnames<-row.names(xcmsObj@phenoData)}
    if(newXCMS)  {sampnames <- xcmsObj@phenoData[[1]]}
    if(!newXCMS) {data12 <- xcms::groupval(xcmsObj, value="into")}
    if(newXCMS)  {data12 <- xcms::featureValues(xcmsObj, value = "into")}
    
    if(taglocation=="filepaths" & !is.null(MStag)){ 
      if(!newXCMS) {
        msfiles<-grep(MStag, xcmsObj@filepaths, ignore.case=TRUE)
        msmsfiles<-grep(idMSMStag, xcmsObj@filepaths, ignore.case=TRUE)
      }
      if(newXCMS) {
        msfiles<-grep(MStag, xcmsObj@processingData@files, ignore.case=TRUE)
        msmsfiles<-grep(idMSMStag, xcmsObj@processingData@files, ignore.case=TRUE)
      }
      if(length(msfiles) == 0) {stop("no msfiles recognized")}
      if(length(msfiles) == 0) {stop("no msmsfiles recognized")}
      if(length(intersect(msfiles, msmsfiles)>0)) {stop("your MS and idMSMStag values do not generate unique file lists")}
      if(length(msfiles)!=length(msmsfiles))  {stop("the number of MS files must equal the number of MSMS files")}
      data1<-t(data12[,msfiles])
      row.names(data1)<-sampnames[msfiles]
      data2<-t(data12[,msmsfiles])
      row.names(data2)<-sampnames[msmsfiles]  ##this may need to be changed to dimnames..
      
      if(normalize == "batch.qc") {
        qc2     <- qc[msmsfiles]
        qc      <- qc[msfiles]
        order2  <- order[msmsfiles]
        order   <- order[msfiles]
        batch2  <- batch[msmsfiles]
        batch   <- batch[msfiles]
      }
      
      if(!newXCMS) {times<-round(xcmsObj@groups[,"rtmed"], digits=3)}
      if(newXCMS) {times<-round(xcmsObj@msFeatureData$featureDefinitions$rtmed, digits = 3)}
      
      if(any(is.na(times))) {
        do<-which(is.na(times))
        for(x in 1:length(do)) {
          if(!newXCMS) {
            times[do[x]]<-  as.numeric((xcmsObj@groups[do[x],"rtmin"]+ xcmsObj@groups[do[x],"rtmax"])/2)
          }
          if(newXCMS) {
            times[do[x]] <- as.numeric((xcmsObj@msFeatureData$featureDefinitions$rtmin[do[x]] + xcmsObj@msFeatureData$featureDefinitions$rtmax[do[x]])/2)
          }
        }
      }
      # if(any(is.na(times))) {stop("na values still present")} else {print("NAs removed")}
      if(!newXCMS) {mzs<-round(xcmsObj@groups[,"mzmed"], digits=4)}
      if(newXCMS) {mzs<-round(xcmsObj@msFeatureData$featureDefinitions$mzmed, digits = 4)}
      
    } else {
      data1<-t(data12)
      msfiles<-1:nrow(data1)
      data2<-t(data12)
      
      row.names(data1)<-sampnames[msfiles]
      row.names(data2)<-sampnames[msfiles]
      
      if(!newXCMS) {times<-round(xcmsObj@groups[,"rtmed"], digits=3)}
      if(newXCMS) {times<-round(xcmsObj@msFeatureData$featureDefinitions$rtmed, digits = 3)}
      
      if(any(is.na(times))) {
        do<-which(is.na(times))
        for(x in 1:length(do)) {
          if(!newXCMS) {
            times[do[x]]<-  as.numeric((xcmsObj@groups[do[x],"rtmin"]+ xcmsObj@groups[do[x],"rtmax"])/2)
          }
          if(newXCMS) {
            times[do[x]] <- as.numeric((xcmsObj@msFeatureData$featureDefinitions$rtmin[do[x]] + xcmsObj@msFeatureData$featureDefinitions$rtmax[do[x]])/2)
          }
        }
      }
      # if(any(is.na(times))) {stop("na values still present")} else {print("NAs removed")}
      if(!newXCMS) {mzs<-round(xcmsObj@groups[,"mzmed"], digits=4)}
      if(newXCMS) {mzs<-round(xcmsObj@msFeatureData$featureDefinitions$mzmed, digits = 4)}
      
    }
    
    if (normalize == "batch.qc") {
      qc2 <- qc[msfiles]
      qc <- qc[msfiles]
      order2 <- order[msfiles]
      order <- order[msfiles]
      batch2 <- batch[msfiles]
      batch <- batch[msfiles]
    }
    
    #if(any(is.na(times))) {stop("na values still present")} else {print("NAs removed")}
    #}
  }
  
  ## if using batch.qc check for proper information
  if(normalize == "batch.qc") {
    if(!all.equal(length(batch), length(qc), length(order), nrow(data1))) {
      stop("all lengths must be identical and are not: ", '\n',
           "  length(batch) = ", length(batch), '\n', 
           "  length(order) = ", length(order), '\n',
           "  length(qc) = ", length(qc), '\n',
           "  number of injections = ", nrow(data1), '\n')
    }
  }
  
  
  history <- paste(history, "RAMClustR (version ",
                   packageVersion("RAMClustR"),
                   ") was utilized to cluster features into spectra (Broeckling 2014).",
                   sep = "")
  
  if(mslev == 2) {
    history <- paste(history,
                     "Feature data included both MS and indiscriminant MS/MS data (Broeckling 2012)."
    )
  }
  ########
  # ensure that we have all numeric values, 
  # then optionally ensure we have all non-zero values in the dataset. 
  # uses a noise addition 'jitter' around minimum values with missing data points.
  # this is mostly necessary for csv input, where other programs may not have used a 'fillPeaks' like step
  rpl1<-unique(c(which(is.na(data1)), which(is.nan(data1)), which(is.infinite(data1))))
  rpl2<-unique(c(which(is.na(data2)), which(is.nan(data2)), which(is.infinite(data2))))
  
  if(length(rpl1)>0) {data1[rpl1]<- 0 }
  if(length(rpl2)>0) {data2[rpl2]<- 0 }
  
  data1[which(data1<0)]<-abs(data1[which(data1<0)])
  data2[which(data2<0)]<-abs(data2[which(data2<0)])

  if(replace.zeros == TRUE) {
    rpl1 <- which(data1==0)
    rpl2 <- which(data2==0)
    if(length(rpl1)>0) {data1[rpl1]<-abs(jitter(rep(min(data1, na.rm=TRUE), length(rpl1) ), amount=min(data1/100, na.rm=TRUE)))}
    if(length(rpl2)>0) {data2[rpl2]<-abs(jitter(rep(min(data2, na.rm=TRUE), length(rpl2) ), amount=min(data2/100, na.rm=TRUE)))}
  }
  
  ########
  # Optional normalization of data, either Total ion signal or quantile
  
  if(!normalize == "none") cat(paste("  normalizing dataset", '\n'))
  
  ## for batch.qc method, we will order the datasets first by batch and run order
  # backup <- data1
  if(normalize == "batch.qc") {
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
    
    ndf <- data.frame(batch = batch2, order = order2, qc = qc2)
    new.ord <- order(ndf$order)
    ndf <- ndf[new.ord,]
    data2 <- data2[new.ord,]
    new.ord <- order(ndf$batch)
    ndf <- ndf[new.ord,]
    data2 <- data2[new.ord,]
    batch2 <- ndf[,"batch"]
    qc2 <- ndf[,"qc"]
    order2 <- ndf[,"order"]
    
  }
  
  data1raw <- data1
  if(mslev > 1) {data2raw <- data2}
  
  
  if(normalize=="TIC") {
    data1<-(data1/rowSums(data1))*mean(rowSums(data1), na.rm=TRUE)
    data2<-(data2/rowSums(data2))*mean(rowSums(data2), na.rm=TRUE)
    
    history <- paste(history,
                     " Features were normalized to total signal using 'tic' normalization."
    )
  }
  
  if(normalize=="quantile") {
    tmpnames1<-dimnames(data1)
    tmpnames2<-dimnames(data2)
    data1<-t(preprocessCore::normalize.quantiles(t(data1)))
    data2<-t(preprocessCore::normalize.quantiles(t(data2)))
    dimnames(data1)<-tmpnames1
    dimnames(data2)<-tmpnames2
    history <- paste(history,
                     " Features were normalized using 'quantile' normalization."
    )
  }
  
  #  data1 <- data1raw
  
  if(normalize == "batch.qc") {
    
    history <- paste(history,
                     " Features were normalized to nearby QC samples on a feature-by-feature basis using the 'batch.qc' option ",
                     "with qc.inj.range = ", qc.inj.range, ".", 
                     " QC normalization was applied to ", length(order), " injections in ", length(unique(batch)), 
                     " bactches, and normization was based on ", length(which(qc)), " recognized QC samples.", sep = "") 
    
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
        # names(doqc) <- names(tmp[doqc])
        
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
        # cat("i: ", i, '\n')
        # cat("f: ", f, '\n')
        # cat("batchmed: ", batchmed, '\n')
        # cat("featmed: ", featmed, '\n')
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
      plot(tmp, col = batch, cex = (qc + 1)/2, ylim = c(0.9,1.11)*range(tmp), 
           main = paste("all:", round(sd(tmp)/mean(tmp), digits = 2), '\n',
                        "qc:", round(sd(tmp[qc])/mean(tmp[qc]), digits = 2)))
      plot(tmpn, col = batch, pch = 19, cex = (qc + 1)/2, , ylim = c(0.9,1.11)*range(tmp), 
           main = paste("all:", round(sd(tmpn)/mean(tmpn), digits = 2), '\n',
                        "qc:", round(sd(tmpn[qc])/mean(tmpn[qc]), digits = 2)))
      rm(tmp); rm(tmpn); rm(tmpnqc); gc()
      
    }
    
    if(mslev == 2) {
      qc    <- qc2
      batch <- batch2
      order <- order2
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
        plot(tmp, col = batch, cex = (qc + 1)/2, ylim = c(0.9,1.11)*range(tmp), 
             main = paste("all:", round(sd(tmp)/mean(tmp), digits = 2), '\n',
                          "qc:", round(sd(tmp[qc])/mean(tmp[qc]), digits = 2)))
        plot(tmpn, col = batch, pch = 19, cex = (qc + 1)/2, , ylim = c(0.9,1.11)*range(tmp), 
             main = paste("all:", round(sd(tmpn)/mean(tmpn), digits = 2), '\n',
                          "qc:", round(sd(tmpn[qc])/mean(tmpn[qc]), digits = 2)))
        rm(tmp); rm(tmpn); rm(tmpnqc); gc()
        
      }
      
    }
    
    dev.off()
  }
  
  
  ########
  # data organization and parsing 
  # sort rt vector and data by retention time
  xcmsOrd<-order(times)
  # cat(times[1:100], '\n')
  # cat('dim data1:', dim(data1), '\n')
  data1<-data1[,order(times)]
  data2<-data2[,order(times)]
  mzs<-mzs[order(times)]
  times<-times[order(times)]
  featnames<-paste(mzs, "_", times, sep="")
  dimnames(data1)[[2]]<-featnames
  dimnames(data2)[[2]]<-featnames
  n<-ncol(data1)
  
  ########
  # calculate similarity matrix
  ramclustObj <- calculate.similarity(n, data1, data2, times, blocksize, mult, maxt, st, sr, rt.only.low.n, cor.method)
  
  ########
  # convert vector to distance formatted object
  ramclustObj <-
    structure(
      ramclustObj,
      Size = (n),
      Diag = FALSE,
      Upper = FALSE,
      method = "RAMClustR",
      Labels = featnames,
      class = "dist"
    )
  gc()
  
  c <- Sys.time()
  cat("RAMClust distances converted to distance object", '\n')
  
  ########
  # cluster using fastcluster package,
  system.time(ramclustObj<-fastcluster::hclust(ramclustObj, method=linkage))
  history <- paste(history,
                   "The feature similarity matrix was clustered using fastcluster package heirarchical clustering method using the",
                   linkage, "method."
  )
  
  gc()
  d<-Sys.time()    
  cat("fastcluster based clustering complete", '\n')
  if(minModuleSize==1) {
    clus<-dynamicTreeCut::cutreeDynamicTree(ramclustObj, maxTreeHeight=hmax, deepSplit=deepSplit, minModuleSize=2)
    sing<-which(clus==0)
    clus[sing]<-max(clus)+1:length(sing)
  }
  if(minModuleSize>1) {
    clus<-dynamicTreeCut::cutreeDynamicTree(ramclustObj, maxTreeHeight=hmax, deepSplit=deepSplit, minModuleSize=minModuleSize)
  }
  gc()
  
  history <- paste(history,
                   " The dendrogram was cut using the cutreeDynamicTree function from the dynamicTreeCut package. ",
                   " Cutting parameters were set to ",
                   "minModuleSize = ", minModuleSize,
                   ", hmax = ", hmax, 
                   ", and deepSplit = ", deepSplit, ".", sep = "")
  history <- paste(history, '\n', '\n')
  
  
  
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
  cat("dynamicTreeCut based pruning complete", '\n')
  
  f<-Sys.time()
  cat(paste("RAMClust has condensed", n, "features into",  max(clus), "spectra", '\n'))
  
  history <- paste(history, n, "features were collapsed into",  max(clus), "spectra.")
  
  ramclustObj$ExpDes<-ExpDes
  strl<-nchar(max(ramclustObj$featclus)) - 1
  ramclustObj$cmpd<-paste("C", formatC(1:length(ramclustObj$clrt), digits = strl, flag = 0 ) , sep="")
  # cat(ramclustObj$cmpd[1:10], '\n')
  ramclustObj$ann<-ramclustObj$cmpd
  ramclustObj$annconf<-rep(4, length(ramclustObj$clrt))
  ramclustObj$annnotes<-rep("", length(ramclustObj$clrt))
  ramclustObj$MSdata_unnormalized <- data1raw
  if(mslev == 2) {
    ramclustObj$MSMSdata_unnormalized <- data2raw
    ramclustObj$MSMSdata<-data2 
  }
  ramclustObj$MSdata<-data1
  
  ########
  # collapse feature dataset into spectrum dataset
  if(collapse=="TRUE") {
    cat("collapsing feature into spectral signal intensities", '\n')
    wts<-colSums(data1[])
    ramclustObj$SpecAbund<-matrix(nrow=nrow(data1), ncol=max(clus))
    for (ro in 1:nrow(ramclustObj$SpecAbund)) { 
      for (co in 1:ncol(ramclustObj$SpecAbund)) {
        ramclustObj$SpecAbund[ro,co]<- weighted.mean(data1[ro,which(ramclustObj$featclus==co)], wts[which(ramclustObj$featclus==co)])
      }
    }
    dimnames(ramclustObj$SpecAbund)[[2]]<-ramclustObj$cmpd
    if(!is.null(ms) & normalize=="quantile") {dimnames(ramclustObj$SpecAbund)[[1]]<-tmpnames1[[1]]}
    if(!usePheno | is.null(xcmsObj)) {dimnames(ramclustObj$SpecAbund)[[1]]<-dimnames(ramclustObj$MSdata)[[1]]} 
    if(usePheno & !is.null(xcmsObj)) {
      if(!newXCMS) {dimnames(ramclustObj$SpecAbund)[[1]]<-as.vector(xcmsObj@phenoData[,1])[msfiles]}
      if(newXCMS) {dimnames(ramclustObj$SpecAbund)[[1]]<-as.vector(xcmsObj@phenoData[[1]])[msfiles]}
    }
    
    g<-Sys.time()
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
    cat(paste("writing msp formatted spectra", '\n'))
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
    cat(paste("msp file complete", '\n')) 
  }  
  ramclustObj$history <- history
  if(nrow(ramclustObj$MSdata) < 5 & rt.only.low.n) {
    warning('\n', "too few samples to use correlational similarity, clustering by retention time only", '\n')
    ramclustObj$history <- paste(ramclustObj$history,
                                 "Since there were fewer than five injections, clustering was performed only using retention time simiilarity.")
  }
  return(ramclustObj)
}

