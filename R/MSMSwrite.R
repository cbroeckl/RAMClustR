#' MSMSwrite
#'
#' Write out MSMS data 
#'
#' This is the Details section
#'
#' @param j numeric a parameter 
#' @return nothing really
#' @author Corey Broeckling
#' @export

MSMSwrite <-
function (j) {
sl<-feat[which(clus==j)]
mz<-as.numeric(unlist(strsplit(sl, "_"))[seq(from=1, to=((2*length(sl))), by=2)])
rt<-as.numeric(unlist(strsplit(sl, "_"))[seq(from=2, to=((2*length(sl))), by=2)])
wts<-rowSums(MSMSdata[,sl])
wm<-vector()
for (k in 1:length(sl)) {
	wm<-c(wm, weighted.mean(MSMSdata[,sl[k]], wts))
	}
mz<-mz[order(wm, decreasing=TRUE)]
rt<-rt[order(wm, decreasing=TRUE)]
wm<-wm[order(wm, decreasing=TRUE)]
mrt<-mean(rt)
write(paste("Name: C", j, "_idMSMS", "_Rt", round(mrt,digits=1), sep=""), libName, append=TRUE)
npeaks<-length(mz)
write(paste("SYNON: $:00in-source", sep=""), libName, append=TRUE)
write(paste("SYNON: $:04", sep=""), libName, append=TRUE)
write(paste("SYNON: $:0515-30", sep=""), libName, append=TRUE)
write(paste("SYNON: $:06Q-TOF", sep=""), libName, append=TRUE)
write(paste("SYNON: $:07Waters Xevo G2 QTOF", sep=""), libName, append=TRUE)
write(paste("SYNON: $:09UPLC: C18 ACN Gradient", sep=""), libName, append=TRUE)
write(paste("SYNON: $:10ESI", sep=""), libName, append=TRUE)
write(paste("SYNON: $:11P", sep=""), libName, append=TRUE)
write(paste("SYNON: $:12Ar", sep=""), libName, append=TRUE)
write(paste("SYNON: $:1450-1200", sep=""), libName, append=TRUE)
write(paste("SYNON: $:1630", sep=""), libName, append=TRUE)
write(paste("Comment: Rt=", round(mrt, digits=2), 
	"  bpmz=", round(mz[which(wm==max(wm))], digits=4), 
	"  Contributor=\"Colorado State University Proteomics and Metabolomics Facility\"", 
	"  Study=", Experiment, 
	sep=""), 
	file=libName, append= TRUE)
write(paste("Num Peaks:", npeaks), file=libName, append= TRUE)
}
