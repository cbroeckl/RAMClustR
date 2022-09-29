#' adap.to.rc
#'
#' use pubchem rest and view APIs to retreive structures, CIDs (if a name or inchikey is given), synonyms, and optionally vendor data, when available. 
#' @details useful for moving from chemical name to digital structure represtation. greek letters are assumed to be 'UTF-8' encoded, and are converted to latin text before searching.   if you are reading in your compound name list, do so with 'encoding' set to 'UTF-8'. 
#' @param seq file name/path to sequence file - expect filenames in column 1 and sample names in column 2.  filenames should match those in spec.abund
#' @param spec.abund file name/path to adap-big export of signal intensities. .csv file expected
#' @param msp file name/path to .msp file created by adap-big (optional)
#' @param annotations  file name/path to annotations .xlsx file.  generally 'simple_export.xlsx' (optional)
#' @param mzdec mz decimals to report for internal storage/reporting.  generally we want 0 for adap kdb low resolution GC-quad instruments.
#' @param min.score 700 (out of 1000) by default
#' @param manual.name when looking up inchikey/names, should manual input be used to fill ambiguous names? generally recommend TRUE, but requires user input. for unsupervised workflows set to FALSE. 
#' @param qc.tag a character string by which to recognize a sample as a qc sample.  i.e. 'QC' or 'qc'. 
#' @param blank.tag a character string by which to recognize a sample as a blank sample.  i.e. 'blank' or 'Blank'.
#' @return returns a ramclustR structured object suitable for down stream processing steps. 
#' @author Corey Broeckling
#' 
#' @export 
#' 

adap.to.rc <- function(
  seq = '', 
  spec.abund = '',  
  msp = '',  
  annotations = '', 
  mzdec = 0, 
  min.score = 700, 
  manual.name = TRUE,
  qc.tag = "qc",
  blank.tag = "blank"
) {
  
  require(MSnbase)
  ## read sequence file into R
  if(is.null(seq)) {
    seq <- read.csv("seq.csv", 
                    header = TRUE, 
                    check.names = FALSE, 
                    stringsAsFactors = FALSE)
  } else {
    seq <-  read.csv(seq, 
                     header = TRUE, 
                     check.names = FALSE, 
                     stringsAsFactors = FALSE)
  }
  
  ## try reading in experimental design
  if(file.exists("ExpDes.Rdata")) load("ExpDes.Rdata")
  if(file.exists("datasets/ExpDes.Rdata")) load("datasets/ExpDes.Rdata")
  if(any(ls()=="ExpDes")) {
    factor.names <- ExpDes[[1]][grep("fact", row.names(ExpDes[[1]])),1]
  }
  
  
  ## read spectra into R and parse
  if(is.null(msp)) {
    msp <- readLines("ADAP_BIG/spectra.msp")
  } else {
    msp <- readLines(msp)
  }
  
  names <- gsub("Name: ", "", msp[grep("Name: ", msp)])
  rts <- as.numeric(gsub("RT: ", "", msp[grep("RT: ", msp)]))
  npeaks <- as.numeric(gsub("Num Peaks: ", "", msp[grep("Num Peaks: ", msp)]))
  ms.start  <- 1 + grep("Num Peaks: ", msp)
  ms.stop <- c(((grep("Name: ", msp)[2:length(names)])-2), (length(msp)-1))
  
  if(length(names) != length(rts)) stop("unequal lengths")
  if(length(names) != length(npeaks)) stop('unequal lengths')
  
  spectra <- rep(NA, length(names))
  names(spectra) <- paste0(
    "C",
    formatC(1:length(spectra), width = nchar(length(spectra)), format = "d", flag = "0")
  )
  spectra <- as.list(spectra)
  
  for(i in 1:length(spectra)) {
    s <- msp[ms.start[i]:ms.stop[i]]
    mz <- round(as.numeric(sapply(1:length(s), FUN = function(x) {
      unlist(strsplit(s[x], " "))[1]
    })), digits = mzdec)
    int <- as.numeric(sapply(1:length(s), FUN = function(x) {
      unlist(strsplit(s[x], " "))[2]
    }))
    
    
    s <- new("Spectrum1",
             mz = mz, 
             polarity = 1L, 
             centroided = TRUE,
             intensity = int, 
             rt = round(60*rts[i], digits = 1))
    spectra[[i]] <- s
    
  }
  
  
  ## read intensities into R and parse
  if(is.null(spec.abund)) {
    spec.abund <- read.csv(
      "ADAP_BIG/specAbund.csv", 
      header = TRUE, 
      check.names = FALSE, 
      stringsAsFactors = FALSE)
  } else {
    spec.abund <- read.csv(
      spec.abund, 
      header = TRUE, 
      check.names = FALSE, 
      stringsAsFactors = FALSE)
  }
  ramclustObj$history$adap.big <- paste(
    paste0(
      "ADAP-BIG (Smirnov 2019) was used to detect peaks, align samples, and deconvolve spectra. ",
      "Results were exported from ADAP and imported into R using the adap.to.rc function. ")
  )
  spec.abund <- spec.abund[, endsWith(names(spec.abund), " area")]
  file.names <- gsub(" area", "", names(spec.abund))
  spec.abund <- t(spec.abund)
  row.names(spec.abund) <- file.names
  dimnames(spec.abund)[[2]] <- names(spectra)
  
  
  
  ## create pheno dataframe
  pheno <- data.frame(
    "filename" = file.names,
    "sample.names" = rep(NA, length(file.names))
  )
  
  for(i in 1:nrow(pheno)) {
    mtch <- grep(unlist(basename(tools::file_path_sans_ext(file.names[i]))), seq[,1], ignore.case = TRUE)
    if(length(mtch) !=1) stop(paste0("missing ", file.names[i], '\n'))
    pheno[i,2] <- seq[mtch,2]
  }
  
  
  ## import annotations file
  if(is.null(annotations)) {
    annotations <- NA
  } else {
    annotations <- xlsx::read.xlsx(
      annotations,  sheetIndex = 2,
      startRow = 2)
    unannotated <- unique(c(
      which(is.na(annotations$Compound.Name)),
      grep("Unknown", annotations$Compound.Name)
    ))
    annotations <- annotations[-unannotated,]
    annotations <- annotations[which(annotations$Fragmentation.Score.by.matching.with.Experimental.spectra > min.score),]
    
    ramclustObj$history$adap.kdb <- paste(
      paste0(
        "ADAP-KDB (Smirnov 2021) was for assigning annotations to library spectra and consensus spectra from Metabolomics Workbench. ",
        "These annotations were exported and imported into R. ")
    )
  } 
  

  cmpd <- names(spectra)
  ann <- rep(NA, length(cmpd))
  adap.score <- rep(NA, length(cmpd))
  if(is.data.frame(annotations)) {
    for(i in 1:length(spectra)) {
      use <- annotations[which(annotations$Numerical.Signal.ID.assigned.by.ADAP.KDB == i),]
      if(length(use) == 0) next
      use <- use[order(use$Fragmentation.Score.by.matching.with.Experimental.spectra, decreasing = TRUE),]
      ann[i] <- use$Compound.Name[1]
      adap.score[i] <- use$Fragmentation.Score.by.matching.with.Experimental.spectra[1]
    }
  }

  # ann <- sapply(1:length(ann), FUN = function(x) trimws(unlist(strsplit( ann[x], "]"))[2]))
  
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
  
  phosinchikey <- rep(NA, length(ann))
  inchi <- rep(NA, length(ann))
  smiles <- rep(NA, length(ann))
  formula <- rep(NA, length(ann))
  pubchem.cid <- rep(NA, length(ann))
  do.met.names <- which(!is.na(met.names))
  
  
  if(length(do.met.names)>0) {
    pc <- rc.cmpd.get.pubchem(cmpd.names = met.names[do.met.names], all.props = FALSE,
                              get.synonyms = FALSE, get.bioassays = FALSE, get.properties = TRUE,
                              manual.entry = manual.name, write.csv = FALSE)
    
    inchikey[do.met.names]    <- pc$properties$InChIKey
    inchi[do.met.names]       <- pc$properties$InChI
    smiles[do.met.names]      <- pc$properties$CanonicalSMILES
    pubchem.cid[do.met.names] <- pc$pubchem$cid
    formula[do.met.names]     <- pc$properties$MolecularFormula
  } else {
    inchikey        <- rep(NA, length(cmpd))
    inchi           <- rep(NA, length(cmpd))
    smiles          <- rep(NA, length(cmpd))
    pubchem.cid     <- rep(NA, length(cmpd))
    formula         <- rep(NA, length(cmpd))
  }

  
  
  is.qc <- grep(qc.tag, pheno[,2])
  is.blank <- grep(blank.tag, pheno[,2])
  spec.abund.samp <- spec.abund[-unique(c(is.qc, is.blank)),, drop = FALSE]
  spec.abund.qc <- spec.abund[is.qc,, drop = FALSE]
  spec.abund.blank <- spec.abund[is.blank, , drop = FALSE]
  phenoData.samp <- pheno[-unique(c(is.qc, is.blank)),, drop = FALSE]
  phenoData.qc <- pheno[is.qc,,drop = FALSE]
  phenoData.blank <- pheno[is.blank,,drop = FALSE]

  
  ramclustObj <- list()
  if(any(ls()=="factor.names")) {
    ramclustObj$factor.names <- factor.names
  }
  ramclustObj$phenoData <- pheno
  ramclustObj$SpecAbund <- apply(spec.abund, 2, FUN = "median", na.rm = TRUE)
  ramclustObj$qc.cv.feature.msdata <- {
    apply(spec.abund[is.qc,], 2, FUN = "sd", na.rm = TRUE)/apply(spec.abund[is.qc,], 2, FUN = "mean", na.rm = TRUE)
  }
  ramclustObj$cmpd <- cmpd
  ramclustObj$clrt <- rts
  met.names[which(is.na(ann))] <- cmpd[which(is.na(ann))]
  ann[which(is.na(ann))] <- cmpd[which(is.na(ann))]
  ramclustObj$ann  <- met.names
  ramclustObj$ann.derivitive <- ann
  ramclustObj$annconf <- "2a"
  ramclustObj$SpecAbund <- spec.abund
  
  # ramclustObj$qc.spec.abund <- spec.abund.qc
  # ramclustObj$qc.phenoData <- phenoData.qc
  # ramclustObj$blank.spec.abund <- spec.abund.blank
  # ramclustObj$blank.phenoData <- phenoData.blank
  ramclustObj$formula <- formula
  ramclustObj$inchikey <- inchikey
  ramclustObj$smiles <- smiles
  ramclustObj$pubchem.cid <- pubchem.cid
  


  
  return(ramclustObj)
  
}
