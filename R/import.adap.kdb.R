#' import.adap.kdb
#'
#' use pubchem rest and view APIs to retrieve structures, CIDs (if a name or inchikey is given), synonyms, and optionally vendor data, when available. 
#' @details useful for moving from chemical name to digital structure representation. greek letters are assumed to be 'UTF-8' encoded, and are converted to latin text before searching.   if you are reading in your compound name list, do so with 'encoding' set to 'UTF-8'. 
#' @param ramclustObj ramclustR object to be annotated.
#' @param annotations  file name/path to annotations .xlsx file.  generally 'simple_export.xlsx'
#' @param min.score 700 (out of 1000) by default
#' @param annotate logical.  TRUE by default.  for now please leave default
#' @param manual.name when looking up inchikey/names, should manual input be used to fill ambiguous names? generally recommend TRUE
#' @return returns a ramclustR structured object suitable for down stream processing steps. 
#' @author Corey Broeckling
#' 
#' @export 
#' 

import.adap.kdb <- function(
  ramclustObj = NULL,
  annotations = NULL,
  min.score = 700,
  annotate = TRUE,
  manual.name = TRUE
) {
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("The use of this function requires package 'readxl'.")
  }
  
  ## import annotations file
  if(is.null(annotations)) {
    annotations <- readxl::read_xlsx(
      "ADAP_BIG/simple_export.xlsx",  sheet = 2,
      skip = 1
    )
  } else {
    annotations <- readxl::read_xlsx(
      annotations,  sheet = 2,
      skip = 1)
  } 
  
  unannotated <- unique(c(
    which(is.na(annotations$Compound.Name)),
    grep("Unknown", annotations$Compound.Name)
  ))
  annotations <- annotations[-unannotated,]
  annotations <- annotations[which(annotations$Fragmentation.Score.by.matching.with.Experimental.spectra >= min.score),]
  
  cmpd <- ramclustObj$cmpd
  ann <- ramclustObj$ann
  adap.score <- rep(NA, length(cmpd))
  for(i in 1:length(cmpd)) {
    use <- annotations[which(annotations$Numerical.Signal.ID.assigned.by.ADAP.KDB == i),]
    if(length(use) == 0) next
    use <- use[order(use$Fragmentation.Score.by.matching.with.Experimental.spectra, decreasing = TRUE),]
    ann[i] <- use$Compound.Name[1]
    adap.score[i] <- use$Fragmentation.Score.by.matching.with.Experimental.spectra[1]
  }
  ann <- sapply(1:length(ann), FUN = function(x) trimws(unlist(strsplit( ann[x], "]"))[2]))
  
  flag <- rep(FALSE, length(ann))
  met.names <- rep(NA, length(ann))
  for(i in 1:length(ann)) { # 182 and 188
    if(is.na(ann[i])) next
    cmpd.name <- ann[i]
    cmpd.name <- trimws(unlist(strsplit(cmpd.name, ",", fixed = TRUE)))
    if(length(cmpd.name) == 1) flag[i] <- TRUE
    cmpd.name <- cmpd.name[!grepl("derivative", cmpd.name)]
    cmpd.name <- cmpd.name[!grepl("TMS", cmpd.name)]
    cmpd.name <- cmpd.name[!grepl("trimethylsilyl", cmpd.name)]
    cmpd.name <- cmpd.name[!grepl("methyloxime", cmpd.name)]
    cmpd.name <- cmpd.name[!grepl("methoxime", cmpd.name)]
    if(length(cmpd.name)>1) flag[i] <- TRUE
    nc <- nchar(cmpd.name)
    cmpd.name <- cmpd.name[which(nc>2)]
    cmpd.name <- paste(cmpd.name, collapse = ",")
    met.names[i] <- cmpd.name
    
    if(manual.name & flag[i]) {
      tmp <- readline(prompt=paste("Enter metabolite name (or 'enter' to accept current: ", met.names[i], "  "))
      if(nchar(tmp) > 0) {
        met.names[i] <- tmp
      }
    }
  }
  
  inchikey <- rep(NA, length(ann))
  inchi <- rep(NA, length(ann))
  smiles <- rep(NA, length(ann))
  formula <- rep(NA, length(ann))
  pubchem.cid <- rep(NA, length(ann))
  do.met.names <- which(!is.na(met.names))
  
  pc <- rc.cmpd.get.pubchem(cmpd.names = met.names[do.met.names], all.props = FALSE,
                                       get.synonyms = FALSE, get.bioassays = FALSE, get.properties = TRUE,
                            manual.entry = manual.name, write.csv = FALSE)
  
  inchikey[do.met.names]    <- pc$properties$InChIKey
  inchi[do.met.names]       <- pc$properties$InChI
  smiles[do.met.names]      <- pc$properties$CanonicalSMILES
  pubchem.cid[do.met.names] <- pc$pubchem$cid
  formula[do.met.names]     <- pc$properties$MolecularFormula
  
  ramclustObj$cmpd <- cmpd
  ramclustObj$ann  <- met.names
  ramclustObj$ann.derivitive <- ann
  ramclustObj$annconf <- "2a"
  ramclustObj$formula <- formula
  ramclustObj$inchikey <- inchikey
  ramclustObj$smiles <- smiles
  ramclustObj$pubchem.cid <- pubchem.cid
  
  return(ramclustObj)
  
}
