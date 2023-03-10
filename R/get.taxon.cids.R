#' get.taxon.cids
#'
#' use pubchem rest to retreive pubchem CIDS known to be found in a given species.  NCBI taxid should be used as input.  i.e. Homo sapiens subsp. 'Denisova' is taxid 741158 
#' @details this function enables return of a list of pubchem CIDs which can be used for prioritizing annotations.  If a genus level taxid is selected, setting the sub.taxa.n option > 0 will return metabolites associated with that taxid and all (assuming n is large enough) subtaxa.  i.e. seting taxid to 9605 (genus = 'Homo') will return metabolites associated with Homo sapiens, Homo heidelbergensis, Homo sapiens subsp. 'Denisova', etc. 
#' @param taxid integer NCBI taxid for the taxon to search.  
#' @param sub.taxa.n integer value for the number of subtaxa to consider.  Note that if the sub.taxa.n value is less the the availabe number of subtaxa, only the first sub.taxa.n values, as reported by rentrez, are returned.  If you require specific subtaxa, you should call those taxids explicitly to ensure those results are returned.  
#' @return returns a vector of integer pubchem cids.  
#' @author Corey Broeckling
#' 
#' @export 
#' 

get.taxon.cids <- function(taxid = 4496, sub.taxa.n = 1000) {
  
  if(sub.taxa.n > 0) {
    sub.taxid <- rentrez::entrez_search(db = "taxonomy", term = paste0("txid", taxid, "[Subtree]"), retmax = sub.taxa.n)
    sub.taxid <- as.integer(sub.taxid$ids)
    if(length(sub.taxid)> 0) {
      taxid <- sort(unique(c(taxid, sub.taxid)))
    }
  }
  

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
  
  return(all.cids)
} 

# oat.cids <- get.taxon.cids(taxid = 4496)
# oat.inchikeys <- rc.cmpd.get.pubchem(cmpd.cid = oat.cids, get.bioassays = FALSE, get.synonyms = FALSE, get.vendors = FALSE)
# oat.inchikeys <- oat.inchikeys$properties$InChIKey

