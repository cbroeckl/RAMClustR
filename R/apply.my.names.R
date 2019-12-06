#' apply.my.names
#'
#' replace DB names with names of your chosing, based on inchikey matching.    
#' @details database chemical names can be a bit odd at times, and sometimes are just database reference numbers. This function takes output from 'get.inchikey.from.names' and uses it to rename compounds based on inchikey or pubchem CID matches. 
#' @param ramclustObj a ramclustR object to name compounds
#' @param cmpds data.frame.  looking for column names 'cmpd.name', 'CID', and 'inchikey'.  one of 'CID' and 'inchikey' must be present. 
#' @return returns a ramclustObj with new annotation names based on your input. 
#' @author Corey Broeckling
#' @export 
apply.my.names <- function(
  ramclustObj = NULL,
  cmpds = NULL
) {
  
  if(is.null(cmpds)) stop("cmpd.names must be a dataframe with column names including 'chem.name' and one of 'CID' or 'inchikey", '\n')
  
  if(!is.data.frame(cmpds)) stop("cmpd.names must be a dataframe with column names including 'chem.name' and one of 'CID' or 'inchikey", '\n')
  
  if(!any(dimnames(cmpds)[[2]] == "cmpd.name")) stop("cmpd.names must be a dataframe with column names including 'cmpd.name' and one of 'CID' or 'inchikey", '\n')
  
  if(!any(dimnames(cmpds)[[2]] == "CID") | !any(dimnames(cmpds)[[2]] == "inchikey")) stop("cmpd.names must be a dataframe with column names including 'chem.name' and one of 'CID' or 'inchikey", '\n')
  
  if(any(dimnames(cmpds)[[2]] == "inchikey")) {
    use = "inchikey"
  } else {use = "CID"}
  
  if(use == "inchikey") {
    for(i in 1:nrow(cmpds)) {
      m <- grep(cmpds[i,use], ramclustObj$inchikey)
      if(length(m)==0) next
      # cat(i, '\n')
      ramclustObj$ann[m] <- make.unique(rep(as.character(cmpds[i,"cmpd.name"]), length(m)))
    }
  }

  if(use == "CID") {
    if(is.null(ramclustObj$pubchem)) stop('please run getSmilesInchi() first to enable pubchem CID matching', '\n')
    for(i in 1:nrow(cmpds)) {
      m <- which(ramclustObj$pubchem[,"CID"] == cmpds[i,"CID"])
      if(length(m)==0) next
      # cat(i, '\n')
      ramclustObj$ann[m] <- make.unique(rep(as.character(cmpds[i,"cmpd.name"]), length(m)))
    }
  }
  
  return(ramclustObj)
}
