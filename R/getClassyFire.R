#' getClassyFire
#'
#' use classyfire web API to look up full ClassyFire heirarchy for each inchikey
#' @details The $inchikey slot is used to look up the 
#' 
#' @param ramclustObj ramclustR object to ClassyFy
#' @param get.all logical; if TRUE, when inchikey classyfire lookup fails, submits for classyfication.  Can be slow and cause errors. default = FALSE.
#' @return returns a ramclustR object.  new dataframe in $classyfire slot with rows equal to number of compounds.  
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'classyFire'
#' @author Corey Broeckling
#' @references Djoumbou Feunang Y, Eisner R, Knox C, Chepelev L, Hastings J, Owen G, Fahy E, Steinbeck C, Subramanian S, Bolton E, Greiner R, and Wishart DS. ClassyFire: Automated Chemical Classification With A Comprehensive, Computable Taxonomy. Journal of Cheminformatics, 2016, 8:61. DOI: 10.1186/s13321-016-0174-y

#' @export 


getClassyFire <- function(
  ramclustObj = RC,
  get.all = FALSE
) {
  
  library(jsonlite)
  
  if(get.all) {
    warning("'get.all' can sometimes take a bit of time", '\n')
  }
  
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
      ## insert a call to 
      if(get.all) {
        if(!is.na(ramclustObj$smiles[i])) {
          params <- list(label = ramclustObj$inchikey[i],
                         query_input = ramclustObj$smiles[i],
                         query_type = "STRUCTURE")
          
          submit <- httr::POST(
            "http://classyfire.wishartlab.com/queries",
            body = params,
            encode = "json",
            httr::accept_json(),
            httr::add_headers("Content-Type" = "application/json")
          )
          
          query_id <- jsonlite::fromJSON(httr::content(submit, 'text')) # %>% unlist() %>% as.list()
          
          out<-list()
          out$classification_status <- "not done"
          Sys.sleep(1)
          while(out$classification_status != "Done") {
            Sys.sleep(1)
            out<- fromJSON(
              paste0(
                "http://classyfire.wishartlab.com/queries/",
                query = query_id$id,
                ".json"
              )
            )
          }
          newinchikey<-out$entities$inchikey
          newinchikey<-gsub("InChIKey=", "", newinchikey)
          out<-fromJSON(paste0(url, "/entities/", newinchikey, ".json"))
        }
      }
    }
    a <- ramclustObj$inchikey[i]
    b <- out$kingdom$name; if(is.null(b)) b<-NA
    c <- out$superclass$name; if(is.null(c)) c<-NA
    d <- out$class$name; if(is.null(d)) d<-NA
    e <- out$subclass$name; if(is.null(e)) e<-NA
    f <- out$direct_parent$name; if(is.null(f)) f<-NA
    g <- out$description; if(is.null(g)) g<-NA
    
    ramclustObj$classyfire[i,]<- c(a,b,c,d,e,f,g)
    rm(out)
  }
  return(ramclustObj)
}


