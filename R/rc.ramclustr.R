#' rc.ramclustr
#'
#' Main clustering function for grouping features based on their analytical behavior.  
#'
#' @param ramclustObj ramclustR object: containing ungrouped features.  constructed by rc.get.xcms.data, for example
#' @param st numeric: sigma t - time similarity decay value 
#' @param sr numeric: sigma r - correlational similarity decay value
#' @param maxt numeric: maximum time difference to calculate retention similarity for - all values beyond this are assigned similarity of zero
#' @param deepSplit logical: controls how agressively the HCA tree is cut - see ?cutreeDynamicTree
#' @param blocksize integer: number of features (scans?) processed in one block  =1000,
#' @param mult numeric: internal value, can be used to influence processing speed/ram usage
#' @param hmax numeric: precut the tree at this height, default 0.3 - see ?cutreeDynamicTree
#' @param minModuleSize integer: how many features must be part of a cluster to be returned? default = 2
#' @param linkage character: heirarchical clustering linkage method - see ?hclust
#' @param cor.method character: which correlational method used to calculate 'r' - see ?cor 'method' option. default = "pearson"
#' @param cor.use character: which data points to use to calculate 'r' - see ?cor 'use' option. default = "pairwise.complete.obs"
#' @param rt.only.low.n logical: default = TRUE  At low injection numbers, correlational relationships of peak intensities may be unreliable.  by defualt ramclustR will simply ignore the correlational r value and cluster on retention time alone.  if you wish to use correlation with at n < 5, set this value to FALSE.
#' @param collapse logical: if true (default), feature quantitative values are collapsed into spectra quantitative values. 
#' @details Main clustering function output - see citation for algorithm description or vignette('RAMClustR') for a walk through.  batch.qc. normalization requires input of three vectors (1) batch (2) order (3) qc.   This is a feature centric normalization approach which adjusts signal intensities first by comparing batch median intensity of each feature (one feature at a time) QC signal intensity to full dataset median to correct for systematic batch effects and then secondly to apply a local QC median vs global median sample correction to correct for run order effects.
#' @return   $featclus: integer vector of cluster membership for each feature
#' @return   $clrt: cluster retention time
#' @return   $clrtsd: retention time standard deviation of all the features that comprise that cluster
#' @return   $nfeat: number of features in the cluster
#' @return   $nsing: number of 'singletons' - that is the number of features which clustered with no other feature
#' @return   $cmpd: compound name.  C#### are assigned in order of output by dynamicTreeCut.  Compound with the most features is classified as C0001...
#' @return   $ann: annotation.  By default, annotation names are identical to 'cmpd' names.  This slot is a placeholder for when annotations are provided
#' @return   $SpecAbund: the cluster intensities after collapsing features to clusters
#' @return   $SpecAbundAve: the cluster intensities after averaging all samples with identical sample names
#' @importFrom graphics abline axis boxplot hist "legend" "par" "plot" "points" "title"
#' @importFrom stats aggregate cor fitted lm loess median predict quantile sd weighted.mean
#' @importFrom utils edit read.csv read.delim read.delim2 write.csv packageVersion
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom e1071 skewness
#' @importFrom gplots heatmap.2
#' @importFrom pcaMethods pca
#' @importFrom jsonlite fromJSON
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
rc.ramclustr  <- function(
  ramclustObj=NULL,
  st=NULL, 
  sr=NULL, 
  maxt=NULL, 
  deepSplit=FALSE, 
  blocksize=2000,
  mult=5,
  hmax=NULL,
  collapse=TRUE,
  minModuleSize=2,
  linkage="average",
  cor.method="pearson",
  cor.use = "pairwise.complete.obs",
  rt.only.low.n = TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("please supply a ramclustR object as input. see rc.get.xcms.data, for example")
  }
    
  ########
  # define ms levels, used several times below
  if(is.null(ramclustObj$MSMSdata)) {
    mslev <- 1
  } else {
    mslev <- 2
  }
  
  a<-Sys.time()   
  
  ########
  # set default parameters of not defined in function call
  
  if(is.null(hmax)) {hmax<-0.3}
  if(is.null(st)) {
    if(!is.null(ramclustObj$st)) {
      st <- ramclustObj$st
    } else {
      # use max feature retention time / 500 as a reasonable approximation
      # not really validated.... 
      st <- max(ramclustObj$frt)/500
    }
  }
  if(is.null(maxt)) {maxt <- 100*st}
  if(is.null(sr)) {sr <- 0.7}
  
  
  params <- c(
    "st"=st, 
    "sr"=sr, 
    "maxt"=maxt, 
    "deepSplit"=deepSplit, 
    "blocksize"=blocksize,
    "mult"=mult,
    "hmax"=hmax,
    "collapse"=collapse,
    "minModuleSize"=minModuleSize,
    "linkage"=linkage,
    "cor.method"=cor.method,
    "rt.only.low.n" = rt.only.low.n
  )

  ########
  # establish some constants for downstream processing
  data1 <- ramclustObj$MSdata
  if(mslev == 2) {
    data2 <- ramclustObj$MSMSdata
  } else {
    data2 <- data1
  }
  
  n<-ncol(data1)
  vlength<-(n*(n-1))/2
  nblocks<-floor(n/blocksize)
  
  times <- ramclustObj$frt
  mzs <- ramclustObj$fmz
  featnames <- paste(ramclustObj$frt, ramclustObj$fmz, sep = "_")
    
  ########
  # make list of all row and column blocks to evaluate
  eval1<-expand.grid(0:nblocks, 0:nblocks)
  names(eval1)<-c("j", "k") #j for cols, k for rows
  eval1<-eval1[which(eval1[,"j"]<=eval1[,"k"]),] #upper triangle only
  bl<-nrow(eval1)
   
  
  tmp.ramclustObj <- calculate.similarity(
    numcols = n, 
    data1 = data1, data2 = data2,
    times = times, blocksize = blocksize, 
    mult = mult, maxt = maxt, 
    st = st, sr = sr, 
    rt.only.low.n = rt.only.low.n, 
    cor.method = cor.method, 
    cor.use = cor.use)
  
  ########
  # convert vector to distance formatted object
  tmp.ramclustObj<-structure(
    tmp.ramclustObj, 
    Size=(n), 
    Diag=FALSE, 
    Upper=FALSE, 
    method="RAMClustR", 
    Labels=featnames, 
    class="dist")
  
  
  tmp.ramclustObj[which(is.na(tmp.ramclustObj))] <- 1
  
  gc()
  
  c<-Sys.time()    
  # cat("RAMClust distances converted to distance object", '\n')
    
  ########
  # cluster using fastcluster package,
  system.time(tmp.ramclustObj<-fastcluster::hclust(
    tmp.ramclustObj, 
    method=linkage))
  
  gc()
  d<-Sys.time()    
  cat("fastcluster based clustering complete", '\n')
  if(minModuleSize==1) {
    clus<-dynamicTreeCut::cutreeDynamicTree(tmp.ramclustObj, maxTreeHeight=hmax, deepSplit=deepSplit, minModuleSize=2)
    sing<-which(clus==0)
    clus[sing]<-max(clus)+1:length(sing)
  }
  if(minModuleSize>1) {
    clus<-dynamicTreeCut::cutreeDynamicTree(
      tmp.ramclustObj, 
      maxTreeHeight=hmax, 
      deepSplit=deepSplit, 
      minModuleSize=minModuleSize)
  }
  gc()
  
  ########
  # build results into ramclustObj
  ramclustObj$merger <- tmp.ramclustObj$merger
  ramclustObj$height <- tmp.ramclustObj$height
  ramclustObj$order <- tmp.ramclustObj$order
  ramclustObj$labels <- tmp.ramclustObj$labels
  ramclustObj$method <- tmp.ramclustObj$method
  ramclustObj$call <- tmp.ramclustObj$call
  ramclustObj$dist.method <- tmp.ramclustObj$dist.method
  ramclustObj$featclus<-clus
  
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
  
  strl<-nchar(max(ramclustObj$featclus)) - 1
  ramclustObj$cmpd<-paste("C", formatC(1:length(ramclustObj$clrt), digits = strl, flag = 0 ) , sep="")
  # cat(ramclustObj$cmpd[1:10], '\n')
  ramclustObj$ann<-ramclustObj$cmpd
  ramclustObj$annconf<-rep(4, length(ramclustObj$clrt))
  ramclustObj$annnotes<-rep("", length(ramclustObj$clrt))
  
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
    g<-Sys.time()
  }
  
  dimnames(ramclustObj$SpecAbund)[[1]] <- ramclustObj$sample_names
  
  if(is.null(dimnames(ramclustObj$SpecAbund)[[1]])) {
    stop('this appears to be an older format ramclustR object and does not have a "phenoData" slot with sample names')
  }
  
  
  if(nrow(ramclustObj$MSdata) < 5 & rt.only.low.n) {
    warning('\n', "too few samples to use correlational similarity, clustering by retention time only", '\n')
  }
  
  ramclustObj$history$ramclustr <- paste0(
    "Features were clustered using the ramclustR algorithm (Broeckling 2014). Parameter settings were as follows: ",
    "st = ", st,
    ", sr = ", sr, 
    ", maxt = ", maxt,
    ", deepSplit = ", deepSplit,
    ", hmax = ", hmax,
    ", minModuleSize = ", minModuleSize,
    ", and cor.method = ", cor.method, ". ",
    if(nrow(ramclustObj$MSdata) < 5 & rt.only.low.n) {
      "There were too few samples to use correlational similarity, thus clustering was performed using retention time only. "
    } else {""}
  )
  
  
  if(is.null(ramclustObj$params)) {ramclustObj$params <- list()}
  ramclustObj$params$rc.ramclustr <- params
  
  return(ramclustObj)
}
