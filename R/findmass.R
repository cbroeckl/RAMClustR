#' findmass
#'
#' see if any features match a given mass, and whether they are plausibly M0
#' @param ramclustObj - the ramclustR object to explore
#' @param mz - numeric mz value to search for
#' @param mztol - absolute mass tolerance around mz
#' @details  a convenience function to perform a targeted search of all feaures for a mass of interest.  Also performs a crude plausibility check as to whether the matched feature could be M0, based on the assumption of approximately 1 carbon per 17 m/z units and natural isottopic abundance of 1.1% 13C
#' @param rttol - when examining isotope patterns, feaure retention time tolerance around features matching mz +- mztol
#' @param m.check - logical - check whether the matching masses are plausibly M0.  Looks for ions above and below the target m/z at the same time that have intensities consistent with target ion being a non-M0 isotope.
#' @return returns a table to the console listing masses which match, their retention time and intensity, and whether it appears to be plausible as M0
#'#' @export

findmass<-function(
  ramclustObj = RC,
  mz = NULL,
  mztol = 0.02,
  rttol = 2,
  m.check = TRUE
) {
  if(is.null(mz)) {stop("must set 'mz'", '\n')}
  if(is.null(mztol)) {stop("must set 'mztol'", '\n')}
  tar<-which(abs(ramclustObj$fmz - mz) <= mztol)
  if(length (tar)==0) {
    stop("no masses found within", mztol, "mz units of", mz, '\n')
  } else {
    out<-data.frame("featn" = tar, 
                    "featclus" = ramclustObj$featclus[tar],
                    "mz" = ramclustObj$fmz[tar],
                    "rt" = ramclustObj$frt[tar],
                    "int" = ramclustObj$msint[tar],
                    "M0_plausible" = rep(NA, length(tar))
                    )
    if(m.check) {
      for(i in 1:length(tar)) {
        ## check to see if there is a signal at mz - proton mass with intensity inconsistent with mz  M0 isotope
        check1<-which((abs(ramclustObj$fmz - mz + 1.0078) <= mztol) & (abs(ramclustObj$frt - ramclustObj$frt[tar[i]]) <= rttol ))
        if(length(check1)>0) {
          negrange <- c(0.5,2)* (ramclustObj$msint[tar[i]]  / ((ramclustObj$fmz[tar[i]] / 17 ) * 0.011))
          res1<-any(ramclustObj$msint[check1] > negrange[1] & ramclustObj$msint[check1] < negrange[2] )
        } else {res1 <- FALSE}
        
        ## check to see if there is a signal at mz + proton mass with intensity inconsistent with mz  M0 isotope
        check2<-which((abs(ramclustObj$fmz - mz - 1.0078) <= mztol) & (abs(ramclustObj$frt - ramclustObj$frt[tar[i]]) <= rttol ))
        if(length(check1)>0) {
          posrange <- c(0.5,2)* (ramclustObj$msint[tar[i]]  * ((ramclustObj$fmz[tar[i]] / 17 ) * 0.011))
          res2<-any(ramclustObj$msint[check2] > negrange[1] & ramclustObj$msint[check2] < negrange[2] )
        } else {res2 <- FALSE}
        out[i, "M0_plausible"]<-!any(c(res1, res2))
      }
    }
  }
  print(out)
}
