get.pathways <- function(cid = 5793) {
  url.pre <- "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22where%22:{%22ands%22:[{%22cid%22:%22"
  url.mid <- "%22},{%22core%22:%221%22}]},%22order%22:[%22taxname,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22CID_"
  url.post <- "_pathway%22}"
  
  d <- read.csv(paste0(url.pre, cid, url.mid, cid, url.post))
  
  if(!is.data.frame(d)) {cat("not a data.frame")}
  return(d)
}

d2 <- get.pathways()
d3 <- get.pathways(cid = 4496 )
head(d3)


## metabolites url https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22consolidatedcompoundtaxonomy%22,%22where%22:{%22ands%22:[{%22taxid%22:%224498%22},{%22srccmpdkind%22:%22Metabolite%22}]},%22order%22:[%22cid,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_4498_consolidatedcompoundtaxonomy%22}
get.taxon.cids <- function(taxid = 4498) {
  
  cids <- vector(length = 0, mode = "integer")
  
  ## metabolites
  url.csv1 <-  paste0(
    "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22consolidatedcompoundtaxonomy%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",
    taxid,
    "%22},{%22srccmpdkind%22:%22Metabolite%22}]},%22order%22:[%22cid,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",
    taxid,
    "_consolidatedcompoundtaxonomy%22}")
  d1 <- read.csv(url.csv1)
  if(any(colnames(d1) == "cid")) {
    cids <- c(cids, as.numeric(d1[which(!(d1[,"cid"] == "NULL")),"cid"]))
  }

  ## natural products
  url.csv2 <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22consolidatedcompoundtaxonomy%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",
    taxid,
    "%22},{%22srccmpdkind%22:%22Natural%20Product%22}]},%22order%22:[%22cid,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",
    taxid,
    "_consolidatedcompoundtaxonomy%22}"
  )
  d2 <- read.csv(url.csv2)
  if(any(colnames(d2) == "cid")) {
    cids <- c(cids, as.numeric(d2[which(!(d2[,"cid"] == "NULL")),"cid"]))
  }
  
  ## metabolic pathway metabolites
  url.csv3 <- paste0("https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?task=pathway_chemical&taxid=",
  taxid,
  "&start=1&limit=10000000&download=true&downloadfilename=TaxID_",
  taxid,
  "_pcget_pathway_chemical&infmt=json&outfmt=csv"
  )
  d3 <- read.csv(url.csv3)
  if(any(colnames(d3) == "cid")) {
    cids <- c(cids, as.numeric(d3[which(!(d3[,"cid"] == "NULL")),"cid"]))
  }
  
  ## glycans
  url.csv4 <- paste0( 
    "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22glycosmos_glycan%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",
    4932,
    "%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",
    4932,
    "_glycosmos_glycan%22}"
  )
  d4 <- read.csv(url.csv4)
  if(any(colnames(d4) == "cid")) {
    cids <- c(cids, as.numeric(d4[which(!(d4[,"cid"] == "NULL")),"cid"]))
  }
  
  return(cids)
} 
  
  