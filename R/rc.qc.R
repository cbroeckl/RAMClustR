#' rc.qc
#'
#' summarize quality control for clustering and for quality control sample variation based on compound ($SpecAbund) and feature ($MSdata and $MSMSdata, if present)
#'
#' @param ramclustObj ramclustR object to analyze
#' @param qc.tag qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.  
#' @param remove.qc logical - if TRUE (default) QC injections will be removed from the returned ramclustObj (applies to $MSdata, $MSMSdata, $SpecAbund, $phenoData, as appropriate). If FALSE, QC samples remain.
#' @param npc number of Principle components to calcuate and plot
#' @param scale "pareto" by default: PCA scaling method used
#' @param outfile.basename base name of output files. Extensions added internally. default = "ramclustQC"
#' @param view.hist logical.  should histograms be plotted? 
#' @details plots a ramclustR summary plot.  first page represents the correlation of each cluster to all other clusters, sorted by retention time.  large blocks of yellow along the diaganol indicate either poor clustering or a group of coregulated metabolites with similar retention time.  It is an imperfect diagnostic, particularly with lipids on reverse phase LC or sugars on HILIC LC systems.  Page 2: histogram of r values from page 1 - only r values one position from the diagonal are used.  Pages 3:5 - PCA results, with QC samples colored red.  relative standard deviation calculated as sd(QC PC scores) / sd(all PC scores).  Page 6: histogram of CV values for each compound int he dataset, QC samples only.  
#' @return   new RC object. Saves output summary plots to pdf and .csv summary tables to new 'QC' directory. If remove.qc = TRUE, moves QC samples to new $QC slot from original position. 
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

rc.qc<-function(ramclustObj=NULL,
                qc.tag="QC",
                remove.qc = FALSE,
                npc=4,
                scale="pareto",
                outfile.basename ="ramclustQC",
                view.hist = TRUE
                
){
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(is.null(qc.tag)) {
    stop("qc.tag = NULL; qc.tag must be defined to enable QC variance examination.", '\n')
  }
  
  if(is.null(outfile.basename)) {
    outfile.basename <- "ramclustQC"
  }
  
  do.sets <- c("MSdata", "SpecAbund")
  
  if(is.null(ramclustObj$SpecAbund)) {
    do.sets <- do.sets[!(do.sets %in% "SpecAbund")]
  } 
  
  do.sets.rows <- sapply(
    c(do.sets, "phenoData"), 
    FUN = function(x) {
      nrow(ramclustObj[[x]])
    })
  
  if(!sd(do.sets.rows) == 0) {
    stop("number of rows in MSdata, SpecAbund, and phenoData sets are not identical.")
  }
  
  ## define QC samples in each set
  if(length(qc.tag) == 1) {
    qc <- grepl(qc.tag[1], ramclustObj$phenoData$sample.names)
  } 
  if(length(qc.tag) == 2) {
    qc <- grepl(qc.tag[1], ramclustObj$phenoData[[qc.tag[2]]])
  }
  
  if(length(which(qc)) == 0) {
    stop("no QC samples found using the qc.tag ", "'", qc.tag, "'", '\n')
  }
  
  ## create directory
  dir.create("QC")
  
  ## if cv threshold for compounds has been applied, use it, else create 'all cmpds' vector
  if(!is.null(ramclustObj$SpecAbund)) {
    if(!is.null(ramclustObj$cmpd.use)) {
      cmpd.use <- ramclustObj$cmpd.use
    } else {
      cmpd.use <- rep(TRUE, length(ramclustObj$ann))
    }
  }
  
  #visualize clustering
  ## if clustering was perfect, we should see a normal distribution of 
  ## correlational r values 1 step from the diagonal
  ## imperfect clustering introduces right skew
  ## load("datasets/RCobject.Rdata")
  if(!is.null(ramclustObj$clrt)) {
    
    ## create file to collect figures. 
    pdf(file=paste("QC/", "ramclust_clustering_diagnostic.pdf", sep=""), 
        useDingbats=FALSE, width=8, height=8)  
    o<-order(ramclustObj$clrt[cmpd.use])
    c<-cor(ramclustObj$SpecAbund[,cmpd.use][,o])
    d<-diag(as.matrix((c[2:(nrow(c)), 1:ncol(c)-1])))
    hist(d, breaks=50, main="")
    title(main="histogram of pearson's r for each cluster to its adjacent cluster (by time)", cex.main=0.8,
          sub=paste("skew =", round(e1071::skewness(d), digits=3), " :values near zero are better", '\n', 
                    'WARNING:metabolic relationships will confound interpretation of this plot'), cex.sub=0.6)
    
    # ideally heatmap will have a bright yellow diagonal with no yellow squares near the diagonal
    # this is slow for larger numbers of clusters
    gplots::heatmap.2(c^2, trace="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, main="pearsons r^2, clusters sorted by rt", cex.main=0.5,
                      cexRow=0.02 + 1/log10(length(o)), cexCol=0.02 + 1/log10(length(o)))
    dev.off()
  }
  
  
  ## PCA of QC samples
  ## histogram of feature and/or compound CVs for QC samples
  
  qc <- which(qc)
  
  cols<-rep(8, nrow(ramclustObj$phenoData))
  cols[qc]<-2
  
  for(x in do.sets) {
    
    ## PCA plot
    if(x == "SpecAbund") {
      td <- ramclustObj[[x]][,cmpd.use]
    } else {
      td <- ramclustObj[[x]]
    }
    
    # if(!is.null(ramclustObj$MSMSdata) & x == "MSdata") {
    #   td <- td + ramclustObj$MSMSdata
    # }
    
    if(min(dim(td)) < npc) {npc <- min(dim(td))}
    PCA<-pcaMethods::pca(td, scale=scale, nPcs=npc, center=TRUE)
    sc<-PCA@scores
    write.csv(sc, file = paste0("QC/", outfile.basename, "_", x, "_pcascores.csv"))
    pdf(file = paste0("QC/", outfile.basename, "_", x, "_qc_diagnostic.pdf"), useDingbats=FALSE, width=8, height=8)  
    
    ld<-PCA@loadings
    for(i in 1:(ncol(sc)-1)) {
      plot(sc[,i], sc[,i+1], col=cols, pch=19, main=paste(
        "PCA analysis:", x, if(x == "SpecAbund") {
          "(compounds)"
        } else {"(features)"}
      ),
      xlab=paste("PC", i, "::  r^2 =", round(PCA@R2[i], digits=2), "::  QC(rel sd) = ", 
                 round(sd(sc[qc,i])/sd(sc[,i]), digits=2) ),
      ylab=paste("PC", i+1, "::  r^2 =", round(PCA@R2[i+1], digits=2), "::  QC(rel sd) = ", 
                 round(sd(sc[qc,i+1])/sd(sc[,i+1]), digits=2) )
      )
      legend(qc.tag, text.col=2, x="topright", bty="n")
    }
    
    ## histogram of QC relative standard deviations for all compounds/clusters
    sds<-apply(td[qc,], 2, FUN="sd", na.rm=TRUE)
    #cat(sds, '\n')
    means<-apply(td[qc,], 2, FUN="mean", na.rm=TRUE)
    cvs<-sds/means
    
    if(x == "MSdata") {
      ramclustObj$qc.cv.feature <- cvs
      ramclustObj$qc.cv.feature.msdata <- cvs
      if(!is.null(ramclustObj$MSMSdata)) {
        sds<-apply(ramclustObj$MSMSdata[qc,], 2, FUN="sd", na.rm=TRUE)
        #cat(sds, '\n')
        means<-apply(ramclustObj$MSMSdata[qc,], 2, FUN="mean", na.rm=TRUE)
        msms.cvs<-sds/means
        ramclustObj$qc.cv.feature.msmsdata <- msms.cvs
        cvs <- pmin(ramclustObj$qc.cv.feature.msdata, msms.cvs)
        ramclustObj$qc.cv.feature <- cvs
      }
      
    } else {
      ramclustObj$qc.cv.cmpd <- cvs
    }
    qs<-quantile(cvs, probs=seq(0,1,0.2), na.rm=TRUE)
    hist(cvs, breaks=50, main="")
    title(paste("histogram of", x,  "CVs from QC samples"), line=2.7)
    title("20% quantiles in red on top axis", col.main =2, cex.main=0.7, line=2)
    axis(side=3, col=2, col.ticks=2, col.axis=2, round(qs, digits=3), labels=TRUE, las=2, cex.axis=0.4)
    dev.off()
    
    if(view.hist) {
      hist(cvs, breaks=50, main="")
      title(paste("histogram of", x,  "CVs from QC samples"), line=2.7)
      title("20% quantiles in red on top axis", col.main =2, cex.main=0.7, line=2)
      axis(side=3, col=2, col.ticks=2, col.axis=2, round(qs, digits=3), labels=TRUE, las=2, cex.axis=0.4)
      
    }
    
    
    if(x == "SpecAbund") {
      out <- data.frame(
        "cmpd" = ramclustObj$cmpd[cmpd.use],
        "annotation" = ramclustObj$ann[cmpd.use],
        "rt" = ramclustObj$clrt[cmpd.use], 
        "mean.int" = means,
        "cv" = cvs
      )
    } else {
      if(is.null(ramclustObj$fmz)) next
      out <- data.frame(
        "mz" = ramclustObj$fmz,
        "rt" = ramclustObj$frt,
        "mean.int" = means,
        "cv" = cvs
      )
      if(length(ramclustObj$labels) > 0) {
        out <- data.frame(
          out, 
          "feature" = ramclustObj$labels,
          "cluster" = ramclustObj$featclus
        )
      }
    }
    write.csv(out, file = paste0("QC/", outfile.basename, "_", x, "_cv_summary.csv"))
  }
  
  if(remove.qc) {
    ramclustObj$qc <- list()
    for(x in c("phenoData", do.sets)) {
      ramclustObj$qc[[x]] <- ramclustObj[[x]][qc,]
      ramclustObj[[x]] <- ramclustObj[[x]][-qc,]
    }
  }
  
  ramclustObj$history$qc.summary <- paste(
    
    "Variance in quality control samples was described using the",
    "rc.qc function within ramclustR. Summary statistics are provided",
    "including the relative standard deviation of QC samples to all",
    "samples in PCA space, as well as the relative standard deviation",
    "of each feature/compound in QC samples, plotted as a histogram.", 
    if(!is.null(ramclustObj$cmpd.use)) {" Only compounds which passed the CV filter are reported."}
  )
  
  return(ramclustObj)
}


