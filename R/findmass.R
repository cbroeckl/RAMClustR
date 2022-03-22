#' findmass
#'
#' see if any features match a given mass, and whether they are plausibly M0
#' @details  a convenience function to perform a targeted search of all feaures for a mass of interest.  Also performs a crude plausibility check as to whether the matched feature could be M0, based on the assumption of approximately 1 carbon per 17 m/z units and natural isotopic abundance of 1.1% 13C.  Note that this function returns the cluster to which the feature is assigned, but that the M0_plausibility is independent of cluster membership.
#' 
#' @param ramclustObj R object: the ramclustR object to explore
#' @param mz numeric: mz value to search for
#' @param mztol numeric: absolute mass tolerance around mz
#' @param rttol numeric: when examining isotope patterns, feature retention time tolerance around features matching mz +- mztol
#' @param zmax integer: maximum charge state to consider.  default is 6.  
#' @param m.check logical:  check whether the matching masses are plausibly M0.  That is, we look for ions 1 proton mass (from charge state 1:zmax) below the target m/z at the same time that have intensities consistent with target ion being a non-M0 isotope.
#' @return returns a table to the console listing masses which match, their retention time and intensity, and whether it appears to be plausible as M0
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

findmass <- function (ramclustObj = NULL, 
                      mz = NULL, 
                      mztol = 0.02, 
                      rttol = 2, 
                      zmax = 6, 
                      m.check = TRUE) 
{
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if (is.null(mz)) {
    stop("must set 'mz'", "\n")
  }
  if (is.null(mztol)) {
    stop("must set 'mztol'", "\n")
  }
  tar <- which(abs(ramclustObj$fmz - mz) <= mztol)
  if (length(tar) == 0) {
    out <- data.frame(featn = NA, featclus = NA, 
                      mz = NA, rt = NA, 
                      int = NA, M0_plausible = NA)
    out<-out[0,]
  }
  else {
    out <- data.frame(featn = tar, featclus = ramclustObj$featclus[tar], 
                      mz = ramclustObj$fmz[tar], rt = ramclustObj$frt[tar], 
                      int = ramclustObj$msint[tar], M0_plausible = rep(NA, 
                                                                       length(tar)))
    if (m.check) {
      for (i in 1:length(tar)) {
        check <- vector()
        for (j in 1:zmax) {
          check1 <- which((abs(ramclustObj$fmz - mz + 
                                 1.007276) <= mztol) & (abs(ramclustObj$frt - 
                                                              ramclustObj$frt[tar[i]]) <= rttol))
          check <- unique(c(check, check1))
        }
        if (length(check) > 0) {
          negrange <- c(0.5, 2) * (ramclustObj$msint[tar[i]]/((ramclustObj$fmz[tar[i]]/17) * 
                                                                0.011))
          out[i, "M0_plausible"] <- !any(ramclustObj$msint[check] > 
                                           negrange[1] & ramclustObj$msint[check] < 
                                           negrange[2])
        }
        else {
          out[i, "M0_plausible"] <- TRUE
        }
      }
    }
  }
  return(out)
}

