#' getSmilesInchi
#'
#' use PubChem API to look up full smiles and inchi notation for each inchikey
#' @details The $inchikey slot is used to look up parameters from pubchem. PubChem CID, a pubchem URL, smiles (canonical) and inchi are returned.  if smiles and inchi slots are alread present (from MSFinder, for example) pubchem smiles and inchi are used to fill in missing values only, not replace. 
#' 
#' @param ramclustObj ramclustR object to look up smiles and inchi for each inchikey (without a smiles/inchi). Must provide one of ramclustObj or inchikey.
#' @param inchikey character vector of inchikey strings.  Must provide one of ramclustObj or inchikey.
#' @param ignore.stereo logical.  default = TRUE. If the Pubchem databases does not have the full inchikey string, should we search by the first (non-stereo) block of the inchikey?  When true, returns the first pubchem match to the inchikey block one string.  If the full inchikey is present, that is used preferentially.
#' @return returns a ramclustR object.  new vector of $smiles and $inchi with length equal to number of compounds.  
#' @importFrom jsonlite fromJSON
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept smiles
#' @concept inchi
#' @concept inchikey
#' @author Corey Broeckling
#' @references Kim S, Thiessen PA, Bolton EE, Bryant SH. PUG-SOAP and PUG-REST: web services for programmatic access to chemical information in PubChem. Nucleic Acids Res. 2015;43(W1):W605-11.

#' @export 


rc.cmpd.get.smiles.inchi <- function(
  ramclustObj = NULL,
  inchikey = NULL,
  ignore.stereo = TRUE
) {
  
  if(is.null(ramclustObj) & is.null(inchikey)) {
    stop("must supply ramclustObj or inchikey vector as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(is.null(ramclustObj)) {
    ramclustObj <- list()
    ramclustObj[['inchikey']] <- inchikey
    ramclustObj$cmpd <- paste0("C", 1:length(inchikey))
  }
  
  params <- c(
    "ramclustObj" = ramclustObj,
    "inchikey" = inchikey,
    "ignore.stereo" = ignore.stereo
  )
  

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
    
    ## start time of loop
    time.start <- Sys.time()
    out<- tryCatch(jsonlite::fromJSON(paste0(url, ramclustObj$inchikey[i], "/JSON")), 
                   error = function(x) {return(NA)})
    
    if(!is.list(out)) {
      if(ignore.stereo) {
        tmp.inchi <- as.character(strsplit(ramclustObj$inchikey[i], "-")[[1]][1])
        out<- tryCatch(jsonlite::fromJSON(paste0(url, tmp.inchi, "/JSON")), 
                       error = function(x) {return(NA)})
      }
    } 
    
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
    # time.start <- Sys.time()
    time.stop <- Sys.time()
    proc.time <- formatC(time.stop - time.start)
    if(attr(proc.time, "units") == "secs") {
      proc.time <- as.numeric(proc.time[1])
      if(proc.time < 0.2) {
        Sys.sleep(0.2 - proc.time)
      }
    }
  }
  
  
  if(is.null(ramclustObj$pubchem.url)) {
    ramclustObj$pubchem.url <- rep(NA, nrow(pubchem))
    do <- which(!is.na(pubchem[,"CID"]))
    ramclustObj$pubchem.url[do] <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem[do,"CID"])
  } else  {
    fix <- which(is.na(ramclustObj$pubchem.url) & !is.na(pubchem$pubchem.url))
    if(length(fix) > 0) {
      ramclustObj$pubchem.url[fix] <- pubchem$pubchem.url[fix]
    }
  }
  
  
  if(is.null(ramclustObj$smiles)) {ramclustObj$smiles <- pubchem$smiles} else {
    fix <- which(is.na(ramclustObj$smiles) & !is.na(pubchem$smiles))
    if(length(fix) > 0) {
      ramclustObj$smiles[fix] <- pubchem$smiles[fix]
    }
  }
  if(is.null(ramclustObj[['inchi']])) {ramclustObj[['inchi']] <- pubchem$inchi} else {
    fix <- which(is.na(ramclustObj$inchi) & !is.na(pubchem$inchi))
    if(length(fix) > 0) {
      ramclustObj$inchi[fix] <- pubchem$inchi[fix]
    }
  }
  if(is.null(ramclustObj$cid)) {ramclustObj$cid <- pubchem$CID} else {
    fix <- which(is.na(ramclustObj$cid) & !is.na(pubchem$CID))
    if(length(fix) > 0) {
      ramclustObj$cid[fix] <- pubchem$CID[fix]
    }
  }
  
  if(!is.null(ramclustObj$inchi) & !is.null(ramclustObj$smiles)) {
    get.inchi <- which(!is.na(ramclustObj$smiles) & is.na(ramclustObj$inchi))
    if(length(get.inchi) > 0) {
      for(i in get.inchi) {
        tmp <- unlist(webchem::cs_convert(ramclustObj$smiles[i], from = "smiles", to = "inchi"))
        if(length(tmp) == 1) {
          if(grepl("InChI=", tmp)) {
            ramclustObj$inchi[i] <- tmp
          }
        } else {
          if(grepl("InChI=", tmp)) {
            ramclustObj$inchi[i] <- tmp[1]
          }
        }
        rm(tmp)
      }
    }
  }
  
  if(is.null) {ramclustObj$params <- list()}
  ramclustObj$params$rc.cmpd.get.smiles.inchi <- params
  
  ramclustObj$history.smiles.inchi <- paste(
    "Smiles structures were retrieved for each inchikey without a structure using the Pubchem API (Djoumbou 2016) called from RAMClustR using the getSmilesInchi function."
  )
  
  
  return(ramclustObj)
}


