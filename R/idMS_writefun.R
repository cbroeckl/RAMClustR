#' idMS_writefun
#'
#' Write out MS data 
#'
#' This is the Details section
#'
#' @param l numeric a parameter 
#' @return nothing really
#' @author Corey Broeckling
#' @export


idMS_writefun <-
function (l) {
	ion<- paste(round(mz[l], digits=4), " ", round(wm[l]), ";", sep="")
	write(ion, file=libName, append= TRUE) }
