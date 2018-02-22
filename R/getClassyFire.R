#' getClassyFire
#'
#' use classyfire web API to look up full ClassyFire heirarchy for each inchikey
#' @details The $inchikey slot is used to look up the 
#' 
#' @param ramclustObj ramclustR object to ClassyFy
#' @return returns a ramclustR object.  new dataframe in $classyfire slot with rows equal to number of compounds.  
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'classyFire'
#' @author Corey Broeckling
#' @references Djoumbou Feunang Y, Eisner R, Knox C, Chepelev L, Hastings J, Owen G, Fahy E, Steinbeck C, Subramanian S, Bolton E, Greiner R, and Wishart DS. ClassyFire: Automated Chemical Classification With A Comprehensive, Computable Taxonomy. Journal of Cheminformatics, 2016, 8:61. DOI: 10.1186/s13321-016-0174-y

#' @export 


getClassyFire <- function(
  ramclustObj = RC
) {
  
  library(jsonlite)
  
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikey slot found, please 'annotate' first", '\n')
  }
  
  ramclustObj$classyfire <- data.frame(
    "inchikey" = rep(NA, length(ramclustObj$inchikey)),
    "kingdom" = rep(NA, length(ramclustObj$inchikey)),
    "superclass" = rep(NA, length(ramclustObj$inchikey)),
    "class" = rep(NA, length(ramclustObj$inchikey)),
    "subclass" = rep(NA, length(ramclustObj$inchikey)),
    "parent" = rep(NA, length(ramclustObj$inchikey)),
    "description" = rep(NA, length(ramclustObj$inchikey))
  )
  
  url = 'http://classyfire.wishartlab.com'
  
  for(i in 1:length(ramclustObj$inchikey)) {
    if(is.na(ramclustObj$inchikey[i])) {next} 
    out<- tryCatch(fromJSON(paste0(url, "/entities/", ramclustObj$inchikey[i], ".json")), 
                   error = function(x) {return(NA)})
    if(length(out)<=1) {
      rm(out)
      next
    }
    a <- ramclustObj$inchikey[i]
    b <- out$kingdom$name; if(is.null(b)) b<-NA
    c <- out$superclass$name; if(is.null(c)) c<-NA
    d <- out$class$name; if(is.null(d)) d<-NA
    e <- out$subclass$name; if(is.null(e)) e<-NA
    f <- out$direct_parent$name; if(is.null(f)) f<-NA
    g <- out$description; if(is.null(g)) g<-NA
    
    ramclustObj$classyfire[i,]<- c(a,b,c,d,e,f,g)
    Sys.sleep(0.2)
  }
  return(ramclustObj)
}


