#' RCQC
#'
#' filter RC object and summarize quality control sample variation  
#'
#' @param ramclustObj ramclustR object to analyze
#' @param qctag "QC" by default - rowname tag to identify QC samples
#' @param npc number of Principle components to calcuate and plot
#' @param scale "pareto" by default: PCA scaling method used
#' @param which.data which dataset to use.  "SpecAbund" by default
#' @param outfile name of output pdf file. 
#'
#' @details plots a ramclustR summary plot.  first page represents the correlation of each cluster to all other clusters, sorted by retention time.  large blocks of yellow along the diaganol indicate either poor clustering or a group of coregulated metabolites with similar retention time.  It is an imperfect diagnostic, particularly with lipids on reverse phase LC or sugars on HILIC LC systems.  Page 2: histogram of r values from page 1 - only r values one position from the diagonal are used.  Pages 3:5 - PCA results, with QC samples colored red.  relative standard deviation calculated as sd(QC PC scores) / sd(all PC scores).  Page 6: histogram of CV values for each compound int he dataset, QC samples only.  
#' @return   new RC object, with QC samples moved to new slot.  prints output summary plots to pdf.
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

RCQC<-function(ramclustObj=NULL,
               qctag="QC",
               npc=4,
               scale="pareto",
               which.data="SpecAbund",
               outfile="ramclustQC.pdf"
               
){
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC1", '\n')
  }
  
  dir.create("QC")
  pdf(file=paste("QC/", "ramclustQC2.pdf", sep=""), useDingbats=FALSE, width=8, height=8)  
  #visualize clustering
  ## if clustering was perfect, we should see a normal distribution of 
  ## correlational r values 1 step from the diagonal
  ## imperfect clustering introduces right skew
  ## load("datasets/RCobject.Rdata")
  if(!is.null(ramclustObj$clrt)) {
    o<-order(ramclustObj$clrt)
    c<-cor(ramclustObj$SpecAbund[,o])
    d<-diag(as.matrix((c[2:(nrow(c)), 1:ncol(c)-1])))
    hist(d, breaks=50, main="")
    title(main="histogram of pearson's r for each cluster to its adjacent cluster (by time)", cex.main=0.8,
          sub=paste("skew =", round(skewness(d), digits=3), " :values near zero are better"), cex.sub=0.6)
    
    # ideally heatmap will have a bright yellow diagonal with no yellow squares near the diagonal
    # this is slow for larger numbers of clusters
    gplots::heatmap.2(c^2, trace="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, main="pearsons r^2, clusters sorted by rt", cex.main=0.5,
              cexRow=0.02 + 1/log10(length(o)), cexCol=0.02 + 1/log10(length(o)))
  }
  ## PCA of QC samples
  while(any(search()=="tf")) {detach(tf)}
  td<-ramclustObj[[which.data]]   ##move as.matrix to ramclustR function
  qc<-grep("QC", dimnames(td)[[1]])
  if(length(qc)>1) { 
    cols<-rep(8, nrow(td))
    cols[qc]<-2
    PCA<-pca(td, scale=scale, nPcs=npc, center=TRUE)
    sc<-PCA@scores
    write.csv(sc, file="QC/RCpcascores.csv")
    ld<-PCA@loadings
    for(i in 1:(ncol(sc)-1)) {
      plot(sc[,i], sc[,i+1], col=cols, pch=19, main="PCA analysis, QC samples vs full set",
           xlab=paste("PC", i, "::  r^2 =", round(PCA@R2[i], digits=2), "::  QC(rel sd) = ", 
                      round(sd(sc[qc,i])/sd(sc[,i]), digits=2) ),
           ylab=paste("PC", i+1, "::  r^2 =", round(PCA@R2[i+1], digits=2), "::  QC(rel sd) = ", 
                      round(sd(sc[qc,i+1])/sd(sc[,i+1]), digits=2) )
      )
      legend(qctag, text.col=2, x="topright", bty="n")
    }
    ## histogram of QC relative standard deviations for all compounds/clusters
    sds<-apply(td[qc,], 2, FUN="sd", na.rm=TRUE)
    #cat(sds, '\n')
    means<-apply(td[qc,], 2, FUN="mean", na.rm=TRUE)
    cvs<-sds/means
    qs<-quantile(cvs, probs=seq(0,1,0.2), na.rm=TRUE)
    hist(cvs, breaks=50, main="", na.rm=TRUE)
    title("histogram of cluster CVs of QC samples", line=2.7)
    title("20% quantiles in red on top axis", col.main =2, cex.main=0.7, line=2)
    axis(side=3, col=2, col.ticks=2, col.axis=2, round(qs, digits=3), labels=TRUE, las=2, cex.axis=0.4)
    nonqc<-ramclustObj[["SpecAbund"]][-grep(qctag, dimnames(ramclustObj[[i]])[[1]]),]
    ramclustObj$clcvqc<-cvs
  } else {nonqc<-ramclustObj[["SpecAbund"]]}
  
  ## histogram of replicate injection relative standard deviations
  keep<-table(row.names(nonqc))
  keep<-names(keep[which(keep>=2)])
  if(length(keep)>0){
    class<-as.factor(keep)
    levs<-levels(class)
    mean1<-matrix(nrow=length(levs), ncol=ncol(nonqc))
    sds<-matrix(nrow=length(levs), ncol=ncol(nonqc))
    for (i in 1:length(levs)){
      mean1[i,]<-apply(nonqc[which(as.character(row.names(nonqc))==levs[i]),], 2, "mean")
      sds[i,]<-apply(nonqc[which(as.character(row.names(nonqc))==levs[i]),], 2, "sd")
    }
    cvs<-apply(sds/mean1, 2, FUN="median", na.rm=TRUE)
    means<-apply(mean1, 2,  FUN="median", na.rm=TRUE)
    write.csv(data.frame(means, cvs), file="QC/cvs.csv")
    #ordmeans<-sort(means, decreasing=TRUE)
    #fivecut<-ordmeans[round(length(means)*0.05)] 
    #up25<-which(means>quantile(means)[4])
    #up5<-which(means>fivecut)
    qs<-quantile(cvs, probs=seq(0,1,0.2), na.rm=TRUE)
    hist(cvs, breaks=50, main="")
    title("histogram of cluster median CVs for replicate injections", line=2.7)
    title("20% quantiles in red on top axis", col.main =2, cex.main=0.7, line=2)
    axis(side=3, col=2, col.ticks=2, col.axis=2, round(qs, digits=3), labels=TRUE, las=2, cex.axis=0.4)
    ramclustObj$clcvrepinj<-cvs
  }
  dev.off()
  for(i in c("SpecAbund", "SpecAbundAve")) {
    if(!is.null(ramclustObj[[i]])) {
      qc<-grep(qctag, dimnames(ramclustObj[[i]])[[1]])
      if(length(qc)>0) {
        ramclustObj[[paste("QC_", i, sep="")]]<-  ramclustObj[[i]][qc,]
      } else {
        ramclustObj[[paste("QC_", i, sep="")]]<-NA
      }
    }
  }
  
  for(i in c("SpecAbund", "SpecAbundAve")) {
    if(!is.null(ramclustObj[[i]])) {
      qc<-grep(qctag, dimnames(ramclustObj[[i]])[[1]])
      if(length(qc)>0) {
        ramclustObj[[i]]<-  ramclustObj[[i]][-qc,]
      } 
    }
  }
  
  ramclustObj$history <- paste(ramclustObj$history, 
                               "Variance in quality control samples was described using the",
                               "RCQC function within ramclustR. Summary statistics are provided",
                               "including the relative standard deviation of QC samples to all",
                               "samples in PCA space, as well as the relative standard deviation",
                               "of each feature in QC samples, plotted as a histogram.")
  
  return(ramclustObj)
}


