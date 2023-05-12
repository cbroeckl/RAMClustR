#' get.taxon.cids
#'
#' use pubchem rest to retreive pubchem CIDS known to be found in a given species.  NCBI taxid should be used as input.  i.e. Homo sapiens subsp. 'Denisova' is taxid 741158 
#' @details this function enables return of a list of pubchem CIDs which can be used for prioritizing annotations.  If a genus level taxid is selected, setting the sub.taxa.n option > 0 will return metabolites associated with that taxid and all (assuming n is large enough) subtaxa.  i.e. seting taxid to 9605 (genus = 'Homo') will return metabolites associated with Homo sapiens, Homo heidelbergensis, Homo sapiens subsp. 'Denisova', etc. 
#' @param taxid integer NCBI taxid for the taxon to search.  
#' @param taxstring taxonomy string for the taxon of interest.
#' @param sub.taxa.n integer value for the number of subtaxa to consider.  Note that if the sub.taxa.n value is less the the availabe number of subtaxa, only the first sub.taxa.n values, as reported by rentrez, are returned.  If you require specific subtaxa, you should call those taxids explicitly to ensure those results are returned.
#' @param get.inchikey logical whether to get the InChIKeys as well (default TRUE).
#' @return returns a vector of integer pubchem cids (and optionally inchikeys if get.inchikey was set to TRUE)
#' @author Corey Broeckling
#' 
#' @export 
#' 

get.taxon.cids <- function(taxid = NULL, taxstring = NULL, sub.taxa.n = 1000, get.inchikey = TRUE) {
  
  ## confirm rentrez
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    stop("The use of this function requires package 'rentrez'.")
  }
  
  ## ensure proper input
  if(is.null(taxid) & is.null(taxstring)) {
    stop("you must submit a valid inter NCBI Taxonomy ID value (taxid) or taxonomy string (taxstring) for the taxon of interest.", '\n',
         " - i.e. taxid = 9606 for Homo sapiens ", '\n',
         " - or taxstring = 'Homo sapiens'", '\n',
         " - or taxstring = 'Homo'", '\n')
  }
  
  ## ensure no ambiguous request
  if(!is.null(taxid) & !is.null(taxstring)) {
    message("you must submit either a valid inter NCBI Taxonomy ID value (taxid) or taxonomy string (taxstring) for the taxon of interest.", '\n',
            " -  taxid will be used in the case that both are submitted.", '\n')
    taxstring <- NULL
  }
  
  ## convert taxstring to taxid
  if(is.null(taxid)) {
    taxid <- rentrez::entrez_search(db = "taxonomy", term = paste0(taxstring, "[All names]"))
    if(length(taxid$ids) == 0) warning("No taxon matched: returning empty dataframe", '\n')
    if(length(taxid$ids) > 1) warning("More than one taxon matched - only the smallest taxid will be used", '\n')
    taxid <- taxid$ids[which.min(taxid$ids)]
  }
  
  ## if sub.taxa.n >0, get subtaxa taxid values
  if(sub.taxa.n > 0) {
    sub.taxid <- rentrez::entrez_search(db = "taxonomy", term = paste0("txid", taxid, "[Subtree]"), retmax = sub.taxa.n)
    sub.taxid <- as.integer(sub.taxid$ids)
    if(length(sub.taxid)> 0) {
      taxid <- sort(unique(c(taxid, sub.taxid)))
    }
  }
  
  ## cleanup connections to ensure no connections remain open
  closePubchemConnections <- function (desc.rem = "pubchem") {
    d <- showConnections(all = TRUE)
    desc <- d[,"description"]
    desc <- desc[grepl(desc.rem, desc)]
    set <- as.integer(as.numeric(names(desc)))
    if(length(set) > 0) {
      for (i in seq_along(set)) close(getConnection(set[i]))
    }
    gc()
    invisible()
  }
  closePubchemConnections()
  
  
  ## collect CIDS for all taxid values
  if(length(taxid) == 0)  {
    out <- data.frame(
      'cid' = all.cids
    )
  } else {
    all.cids <- vector(length = 0, mode = "integer")
    
    cat("retrieving metabolites for taxid: ", '\n')
    
    for(i in 1:length(taxid)) {
      cat(taxid[i], " ")
      cids <- vector(length = 0, mode = "integer")
      ## metabolites
      url.csv1 <-  paste0(
        "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22consolidatedcompoundtaxonomy%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",
        taxid[i],
        "%22},{%22srccmpdkind%22:%22Metabolite%22}]},%22order%22:[%22cid,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",
        taxid[i],
        "_consolidatedcompoundtaxonomy%22}")
      d1 <- read.csv(url.csv1)
      if(any(colnames(d1) == "cid")) {
        cids <- c(cids, as.numeric(d1[which(!(d1[,"cid"] == "NULL")),"cid"]))
      }
      
      ## natural products
      url.csv2 <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22consolidatedcompoundtaxonomy%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",
        taxid[i],
        "%22},{%22srccmpdkind%22:%22Natural%20Product%22}]},%22order%22:[%22cid,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",
        taxid[i],
        "_consolidatedcompoundtaxonomy%22}"
      )
      d2 <- read.csv(url.csv2)
      if(any(colnames(d2) == "cid")) {
        cids <- c(cids, as.numeric(d2[which(!(d2[,"cid"] == "NULL")),"cid"]))
      }
      
      ## metabolic pathway metabolites
      url.csv3 <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?task=pathway_chemical&taxid=",
        taxid[i],
        "&start=1&limit=10000000&download=true&downloadfilename=TaxID_",
        taxid[i],
        "_pcget_pathway_chemical&infmt=json&outfmt=csv"
      )
      d3 <- read.csv(url.csv3)
      if(any(colnames(d3) == "cid")) {
        cids <- c(cids, as.numeric(d3[which(!(d3[,"cid"] == "NULL")),"cid"]))
      }
      
      ## glycans
      url.csv4 <- paste0( 
        "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22glycosmos_glycan%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",
        taxid[i],
        "%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",
        taxid[i],
        "_glycosmos_glycan%22}"
      )
      d4 <- read.csv(url.csv4)
      if(any(colnames(d4) == "cid")) {
        cids <- c(cids, as.numeric(d4[which(!(d4[,"cid"] == "NULL")),"cid"]))
      }
      
      cids <- unique(cids)
      all.cids <- sort(unique(c(all.cids, cids)))
    }
    
    out <- data.frame(
      'cid' = all.cids
    )
  }
  
  if(get.inchikey & !any(is.na(out))) {
    all.inchikeys <- rep(NA, length(all.cids))
    do.ind <- split(1:length(all.cids), ceiling(seq_along(1:length(all.cids))/100))
    do.l <- split(all.cids, ceiling(seq_along(all.cids)/100))
    for(i in 1:length(do.l)) {
      keep <- which(!all.cids[do.ind[[i]]]=="NA")
      # cat(do.l[[i]][keep], '\n')
      Sys.sleep(0.2)
      if(length(keep)==0) next
      # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularWeight/TXT
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(all.cids[do.ind[[i]]][keep], collapse = ","),
                     "/property/", "inchikey", "/JSON")
      out <- tryCatch(
        {
          jsonlite::read_json(html)
        },
        error=function(cond) {
          closePubchemConnections()
          return(NA)
        },
        warning=function(cond) {
          closePubchemConnections()
          return(NA)
        },
        finally={
          closePubchemConnections()
          cons <- suppressWarnings(showConnections(all = TRUE)); rm(cons)
        }
      )
      if(is.na(out[1])) next
      tmp <- lapply(1:length(out$PropertyTable$Properties),
                    FUN = function(x) {
                      unlist(out$PropertyTable$Properties[[x]])
                    })
      tmp <- data.frame(t(data.frame(tmp, stringsAsFactors = FALSE)), stringsAsFactors = FALSE)
      all.inchikeys[do.ind[[i]][keep]] <- tmp$InChIKey[keep]
    }
    
    out <- data.frame(
      'cid' = all.cids, 
      'inchikey' = all.inchikeys
    )
    
    
  }  else {
    out <- cbind(out, "inchikey" = NA)
  }
  
  message(paste("retreived", nrow(out), "structures"), '\n')
  return(out)
}


