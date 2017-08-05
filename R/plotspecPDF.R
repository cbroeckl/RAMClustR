
#' ramclustR
#'
#' Main clustering function 
#'
#' This is the Details section
#'
#' @param ramclustObj=RC,
#' @param spec=NULL, 
#' @param mslev=as.numeric(as.character(RC$ExpDes[[2]]["MSlevs", 1])),  
#' @param mzdiff=2.5
#'
#' @return NULL
#' @author Corey Broeckling

plotspecPDF<-function(
	ramclustObj=RC,
	spec=NULL, 
	mslev=as.numeric(as.character(RC$ExpDes[[2]]["MSlevs", 1])),
	mzdiff=2.5
	) {
	dir.create("pdf_spectra")
	if(grepl(ramclustObj$ExpDes[[2]]["mstype",1], 'quad', ignore.case=TRUE)) digits<-0
	if(grepl(ramclustObj$ExpDes[[2]]["mstype",1], 'tof', ignore.case=TRUE)) digits<-3
	if(grepl(ramclustObj$ExpDes[[2]]["mstype",1], 'orbe', ignore.case=TRUE)) digits<-4
	
	if(is.null(spec)) spec<-c(1, length(ramclustObj$cmpd))
	digs<-nchar(max(spec))
	for (j in spec) {
		cmpd<-j
		pdf(file=paste("pdf_spectra/", formatC(cmpd, width=digs, format="d", flag = 0), ".pdf", sep=""), useDingbats=FALSE, height=5*mslev, width=7)
		par(mfrow=c(mslev,1))
			for (m in 1:length(mslev)){
	      		sl<-which(ramclustObj$featclus==cmpd)
     				mz<-ramclustObj$fmz[sl]
      			rt<-ramclustObj$frt[sl]
      			mrt<-mean(rt)
      			npeaks<-length(mz)
				if(m==1) y<-ramclustObj$msint[sl]
				if(m==2) y<-ramclustObj$msmsint[sl]
				mz<-mz[order(y, decreasing=TRUE)]
				y<-y[order(y, decreasing=TRUE)]
				mzlab<-mz
				ylab<-y
				mzdif<-abs(outer(mzlab, mzlab, "-"))
				rem<-vector()
				for(i in 1:(ncol(mzdif)-1)) {
					do<-mzdif[i:nrow(mzdif), i]
					rem<-union(rem, (i+which(do!=0 & do <mzdiff)-1))
					}
				if(length(rem)>0) {
					mzlab<-mzlab[-rem]
					ylab<-ylab[-rem]
					}
				plot(mz, y, type="h", xlab="m/z", ylab="intensity", main=paste("Compound ", j, "; rt = ", round(mrt), "sec", sep=""), ylim=c(0, max(y)*1.2), yaxs="i", xlim=c(0.7*min(mz), 1.3*max(mz)))
				text(mzlab, ylab, round(mzlab, digits=digits), cex=0.6, pos=4, srt=30, offset=0.25)	
			}
		dev.off()
		}				
	}
