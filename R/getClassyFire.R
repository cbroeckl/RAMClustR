#' getClassyFire
#'
#' use classyfire web API to look up full ClassyFire heirarchy for each inchikey
#' @details The $inchikey slot is used to look up the 
#' 
#' @param ramclustObj ramclustR object to ClassyFy
#' @param get.all logical; if TRUE, when inchikey classyfire lookup fails, submits for classyfication.  Can be slow. max.wait (below) sets max time to spend on each compound before moving on. default = FALSE.
#' @param max.wait  numeric; maximum time to wait per compound when 'get.all' = TRUE.   
#' @param posts.per.minute  integer; a limit set when 'get.all' is true.  ClassyFire server accepts no more than 5 posts per minute when calculating new ClassyFire results.  Slows down submission process to keep server from denying access.  
#' @return returns a ramclustR object.  new dataframe in $classyfire slot with rows equal to number of compounds.  
#' @importFrom jsonlite fromJSON
#' @importFrom RCurl url.exists
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @concept classyFire
#' @author Corey Broeckling
#' @references Djoumbou Feunang Y, Eisner R, Knox C, Chepelev L, Hastings J, Owen G, Fahy E, Steinbeck C, Subramanian S, Bolton E, Greiner R, and Wishart DS. ClassyFire: Automated Chemical Classification With A Comprehensive, Computable Taxonomy. Journal of Cheminformatics, 2016, 8:61. DOI: 10.1186/s13321-016-0174-y

#' @export 

getClassyFire <- function (ramclustObj = NULL, get.all = TRUE, max.wait = 10, posts.per.minute = 5) 
{
  if (is.null(ramclustObj$inchikey)) {
    stop("no inchikey slot found, please 'annotate' first", 
         "\n")
  }
  if (get.all & is.null(ramclustObj$smiles)) {
    stop("obtaining new classyfication (get.all option) requires a smiles notation", 
         "\n", "and no smiles slot found in ramclustObj.  Please first run 'getSmilesInchi()'", 
         "\n")
  }
  if (any(names(ramclustObj) == "classyfire")) {
    redo <- TRUE} else {redo <- FALSE}
  if (!redo) {
    ramclustObj$classyfire <- data.frame(inchikey = rep(NA, 
                                                        length(ramclustObj$inchikey)), kingdom = rep(NA, 
                                                                                                     length(ramclustObj$inchikey)), superclass = rep(NA, 
                                                                                                                                                     length(ramclustObj$inchikey)), class = rep(NA, length(ramclustObj$inchikey)), 
                                         subclass = rep(NA, length(ramclustObj$inchikey)), 
                                         parent = rep(NA, length(ramclustObj$inchikey)), 
                                         description = rep(NA, length(ramclustObj$inchikey)))
  }
  url = "http://classyfire.wishartlab.com"
  for (i in 1:length(ramclustObj$inchikey)) {
    if (is.na(ramclustObj$inchikey[i])) {
      next
    }
    out <- tryCatch(jsonlite::fromJSON(paste0(url, "/entities/", ramclustObj$inchikey[i], 
                                    ".json")), error = function(y) {
                                      return(NA)
                                    })
    if (length(out) > 1) {
      a <- ramclustObj$inchikey[i]
      b <- out$kingdom$name
      if (is.null(b)) 
        b <- NA
      c <- out$superclass$name
      if (is.null(c)) 
        c <- NA
      d <- out$class$name
      if (is.null(d)) 
        d <- NA
      e <- out$subclass$name
      if (is.null(e)) 
        e <- NA
      f <- out$direct_parent$name
      if (is.null(f)) 
        f <- NA
      g <- out$description
      if (is.null(g)) 
        g <- NA
      ramclustObj$classyfire[i, ] <- c(a, b, c, d, e, 
                                       f, g)
      rm(out)
    }
  }
  if (get.all) {
    get.full <- which(!is.na(ramclustObj$inchikey) & is.na(ramclustObj$classyfire[, 
                                                                                  2]))
    cat("this will take some time - maximum of", posts.per.minute, "posts per minute", '\n')
    for (i in get.full) {

      cat(i)
      if (!is.na(ramclustObj$smiles[i])) {
        params <- list(label = "ramclustR", query_input = ramclustObj$smiles[i], 
                       query_type = "STRUCTURE")
        submit <- httr::POST("http://classyfire.wishartlab.com/queries", 
                             body = params, encode = "json", httr::accept_json(), 
                             httr::add_headers(`Content-Type` = "application/json"))
        Sys.sleep(1)
        query_id <- jsonlite::fromJSON(httr::content(submit, 
                                                     "text"))
        
        if(any(names(query_id) == "error")) {
          if(query_id$error == "Limit exceeded") {
            Sys.sleep(60)
            params <- list(label = "ramclustR", query_input = ramclustObj$smiles[i], 
                           query_type = "STRUCTURE")
            submit <- httr::POST("http://classyfire.wishartlab.com/queries", 
                                 body = params, encode = "json", httr::accept_json(), 
                                 httr::add_headers(`Content-Type` = "application/json"))
            query_id <- jsonlite::fromJSON(httr::content(submit, 
                                                         "text"))
          }
        }
        
        
        if(any(names(query_id) == "status")) {
          if(query_id$status == "500") {
            cat(" failed", '\n')
            next
          }
        }
        
        out <- list()
        out$classification_status <- "not done"
        out$number_of_elements <- 0
        Sys.sleep(1)
        skiptonext <- FALSE
        time.a <- Sys.time()
        while (out$number_of_elements == 0) {
          Sys.sleep(1)
          
          if(!RCurl::url.exists(paste0("http://classyfire.wishartlab.com/queries/",  query = query_id$id, ".json"))) {
            cat(" not done", '\n')
            break
          }
          out <- tryCatch( 
            {
            out <- jsonlite::fromJSON(paste0("http://classyfire.wishartlab.com/queries/", query = query_id$id, ".json"))
            }, 
            error = function(y) {
              out <- list()
              out$classification_status <- "not done"
              out$number_of_elements <- 0
              out
              }
            )
          
          if (round(as.numeric(difftime(Sys.time(), 
                                        time.a, units = "secs")), 3) >= max.wait) {
            cat("timed out", "\n")
            skiptonext <- TRUE
            break
          }
        }
        if (out$number_of_elements == 0) {
          ramclustObj$classyfire[i, ] <- c(ramclustObj$inchikey[i], 
                                           rep(NA, 6))
          rm(out)} else {
            a <- out$entities$inchikey
            if(is.null(a)) {
              a <- ramclustObj$inchikey[i]
            } else {
              a <- gsub("InChIKey=", "", a)
            }
            b <- out$kingdom$name
            if (is.null(b)) {
              b <- NA
              c <- NA
              d <- NA 
              e <- NA
              f <- NA
              g <- NA
            } else {
              c <- out$superclass$name
              if (is.null(c)) {c <- NA}
              d <- out$class$name
              if (is.null(d)) {d <- NA}
              e <- out$subclass$name
              if (is.null(e)) {e <- NA}
              f <- out$direct_parent$name
              if (is.null(f)) {f <- NA}
              g <- out$description
              if (is.null(g)) {g <- NA}
            }
            cat(" done", '\n')
            ramclustObj$classyfire[i, ] <- c(a, b, c, 
                                             d, e, f, g)
            rm(out)
          }
      }
      Sys.sleep(ceiling(60/posts.per.minute))
    }
  }
  
  ramclustObj$history <- paste(ramclustObj$history, 
                               "Compounds were assigned to chemical ontogenies using the ClassyFire API (Djoumbou 2016).")
  
  
  return(ramclustObj)
}

