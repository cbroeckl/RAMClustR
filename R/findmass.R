findmass<-function(
  ramclustObj = RC,
  mz = NULL,
  mztol = 0.02
) {
  if(is.null(mz)) {stop("must set 'mz'", '\n')}
  tar<-which(abs(ramclustObj$mz - mz) <= mztol)
}