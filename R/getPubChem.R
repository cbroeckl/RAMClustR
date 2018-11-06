#' getSmilesInchi
#'
#' use PubChem API to look up full smiles and inchi notation for each inchikey
#' @details The $inchikey slot is used to look up parameters from pubchem. PubChem CID, a pubchem URL, smiles (canonical) and inchi are returned.  if smiles and inchi slots are alread present (from MSFinder, for example) pubchem smiles and inchi are used to fill in missing values only, not replace. 
#' 
#' @param ramclustObj ramclustR object to ClassyFy
#' @return returns a ramclustR object.  new dataframe in $classyfire slot with rows equal to number of compounds.  
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'PubChem', smiles, inchi, inchikey
#' @author Corey Broeckling
#' @references Kim S, Thiessen PA, Bolton EE, Bryant SH. PUG-SOAP and PUG-REST: web services for programmatic access to chemical information in PubChem. Nucleic Acids Res. 2015;43(W1):W605-11.

#' @export 


getSmilesInchi <- function(
  ramclustObj = NULL
) {
  
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikey slot found, please 'annotate' first", '\n')
  }
  
  
  url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
  
  params <- c("CID", "inchi", "smiles")
  
  pubchem <- data.frame(matrix(nrow = length(ramclustObj$cmpd), ncol = length(params)))
  dimnames(pubchem)[[1]] <- ramclustObj$cmpd
  dimnames(pubchem)[[2]] <- params
  
  for(i in 1:length(ramclustObj$inchikey)) {
    if(is.na(ramclustObj$inchikey[i])) {next} 
    out<- tryCatch(fromJSON(paste0(url, ramclustObj$inchikey[i], "/JSON")), 
                   error = function(x) {return(NA)})
    
    if(is.list(out)) {
      
      pubchem[i,"CID"] <- out$PC_Compounds$id[1,"id"]
      
      prop.names  <-out$PC_Compounds$props[[1]][[1]]
      prop.values <- out$PC_Compounds$props[[1]][[2]]
      
      sm <- grep("smiles", prop.names[,"label"], ignore.case = TRUE)
      if(length(sm) > 1) {
        can <- grep("canonical", prop.names[,"name"], ignore.case = TRUE)
        if(length(can) == 0) {sm <- sm[1]} else {
          if(length(intersect(sm, can)) == 1)  {
            sm <-  intersect(sm, can)
          } else {sm <- sm[1]}
        }
        pubchem[i,"smiles"] <- prop.values[sm,"sval"]
      }
      
      
      
      inchi <- which(
        grepl("inchi", prop.names[,"label"], ignore.case = TRUE) & 
          !grepl("inchikey", prop.names[,"label"], ignore.case = TRUE)
      )
      if(length(inchi == 1)) {
        pubchem[i,"inchi"] <- prop.values[inchi,"sval"]
      }
    }
    
  }
  
  ramclustObj$pubchem <- pubchem
  ramclustObj$pubchem.url <- rep(NA, nrow(pubchem))
  do <- which(!is.na(pubchem[,"CID"]))
  ramclustObj$pubchem.url[do] <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem[do,"CID"])

  if(is.null(ramclustObj$smiles)) {ramclustObj$smiles <- pubchem$smiles} else {
    fix <- which(is.na(ramclustObj$smiles) & !is.na(pubchem$smiles))
    if(length(fix) > 0) {
      ramclustObj$smiles[fix] <- pubchem$smiles[fix]
    }
  }
  if(is.null(ramclustObj$inchi)) {ramclustObj$inchi <- pubchem$inchi} else {
    fix <- which(is.na(ramclustObj$inchi) & !is.na(pubchem$inchi))
    if(length(fix) > 0) {
      ramclustObj$inchi[fix] <- pubchem$inchi[fix]
    }
  }
  
  return(ramclustObj)
}


