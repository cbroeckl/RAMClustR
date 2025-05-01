#' rc.merge.split.clusters
#'
#' Cluster refinement - scanning instruments (quadrupole, as in GC-MS) can display cluster splitting, possibily due to slight differences in measured peak retentiont time as a function of mass due to scan dynamics.  this function enables a second pass clustering designed to merge two clusters if the second cluster is within a small retention time window and shows a sufficiently strong correlation.  
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param merge.threshold numeric. value between -1 and 1 indicating the correlational r threshold above which two clusters will be merged
#' @param cor.method character: which correlational method used to calculate 'r' - see ?cor 'method' option. default = "pearson"
#' @param cor.use character: which data points to use to calculate 'r' - see ?cor 'use' option. default = "pairwise.complete.obs"
#' @param rt.sd.factor numeric.  default = 3.  clusters within rt.sd.factor * ramclustObj$rtsd (cluster retention time standard deviation) are considered for merging.
#' @param sample.name.column character. column name from ramclustObj$phenoData which should be used as row.names in new ramclustObj$SpecAbund dataset.  
#' @details exports files to a directory called 'spectra'.  If one.file = FALSE, a new directory 'spectra/msp' is created to hold the individual msp files. if do.findman has been run, spectra are written as ms2 spectra, else as ms1. 
#' @return new ramclustR object, with (generally) fewer clusters than the input ramclustR object.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept xcms
#' @author Corey Broeckling
#' @importFrom methods is
#' @export 
#' 
#' 
rc.merge.split.clusters <- function(
    ramclustObj = NULL,
    merge.threshold = 0.7,
    cor.method = 'spearman',
    rt.sd.factor = 3,
    cor.use = "everything",
    sample.name.column = "sample.ID"
    
) {
  
  if(is.null(ramclustObj)) stop('please provide a ramclustObj', '\n')
  if(!is.numeric(merge.threshold)) merge.threshold <- as.numeric(merge.threshold)
  if(merge.threshold < 0 | merge.threshold > 1) stop("'merge.threshold' must be between zero and one", '\n')
  
  
  ## for all clusters, see if there are clusters at nearby retentiontimes which highly correlate
  ## if there are, join them, renumber featclus, and regenerate SpecAbund file.  
  orig.cl.n <- max(ramclustObj$featclus)
  cls <- max(ramclustObj$featclus):2
  # ramclustObj <- RC
  for(i in cls) {
    if(max(diff(sort(unique(ramclustObj$featclus)))) > 1) stop("error 1")
    potential.merges <- which(
      (abs(ramclustObj$clrt[1:i] - ramclustObj$clrt[i])) < (rt.sd.factor*ramclustObj$clrtsd[i])
    )
    potential.merges <- potential.merges[!potential.merges == i]
    if(length(potential.merges) == 0) next
    
    rval <- cor(ramclustObj$SpecAbund[,i], ramclustObj$SpecAbund[,potential.merges], method = cor.method, use = cor.use)
    merges <- potential.merges[which(rval >= merge.threshold)]
    if(length(merges) == 0) next
    
    merges <- max(merges)
    
    ramclustObj$featclus[which(ramclustObj$featclus == i)] <- merges
    # if(max(diff(sort(unique(ramclustObj$featclus)))) > 1) cat(i, "max diff = ", max(diff(sort(unique(ramclustObj$featclus)))), '\n')
    old.featclus <- ramclustObj$featclus
    new.featclus <- old.featclus
    decend.by.one <- which(old.featclus > i)
    new.featclus[decend.by.one] <- (old.featclus[decend.by.one])-1
    if(max(diff(sort(unique(new.featclus)))) > 1) stop("error 2")
    ramclustObj$featclus <- new.featclus
    
  }
  # sort(unique(ramclustObj$featclus))
  # sort(unique(RC$featclus))
  
  ## store SpecAbund sample names
  sa.rn <- dimnames(ramclustObj$SpecAbund)[[1]]
  
  # 0.99888
  
  # collapse feature dataset into spectrum dataset
  data1 <- ramclustObj$MSdata
  wts<-colMeans(data1[], na.rm = TRUE)
  ramclustObj$SpecAbund<-matrix(nrow=nrow(data1), ncol=max(ramclustObj$featclus))
  for (ro in 1:nrow(ramclustObj$SpecAbund)) { 
    for (co in 1:ncol(ramclustObj$SpecAbund)) {
      ramclustObj$SpecAbund[ro,co]<- weighted.mean(data1[ro,which(ramclustObj$featclus==co)], wts[which(ramclustObj$featclus==co)], na.rm = TRUE)
    }
  }
  dimnames(ramclustObj$SpecAbund)[[1]] <- sa.rn
  
  strl <- nchar(max(ramclustObj$featclus)) - 1
  ramclustObj$cmpd <- paste("C", formatC(1:max(ramclustObj$featclus), digits = strl, flag = 0 ) , sep="")
  ramclustObj$ann <- ramclustObj$cmpd
  
  clrt<-aggregate(ramclustObj$frt, by=list(ramclustObj$featclus), FUN="mean")
  ramclustObj$clrt<-clrt[which(clrt[,1]!=0),2]
  clrtsd<-aggregate(ramclustObj$frt, by=list(ramclustObj$featclus), FUN="sd")
  ramclustObj$clrtsd<-clrtsd[which(clrtsd[,1]!=0),2]
  ramclustObj$nfeat<-as.vector(table(ramclustObj$featclus)[2:max(ramclustObj$featclus)])
  ramclustObj$nsing<-length(which(ramclustObj$featclus==0))
  ramclustObj$annconf<-rep(4, length(ramclustObj$clrt))
  ramclustObj$annnotes<-rep("", length(ramclustObj$clrt))
  dimnames(ramclustObj$SpecAbund)[[1]] <- ramclustObj$phenoData[,sample.name.column]
  
  new.cl.n <- max(ramclustObj$featclus)
  cat(paste("Original cluster number =", orig.cl.n, '\n', "New cluster number =", new.cl.n, '\n'))
  
  return(ramclustObj)
  
}

