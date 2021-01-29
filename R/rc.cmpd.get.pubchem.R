#' rc.cmpd.get.pubchem
#'
#' use PubChem API to look up pubchem data for compounds assigned annotations.
#' @details Uses pubchem rest to retrieve standard pubchem properties for a compound.  The $inchikey slot and, if available, the MSFinder.structure CID value are used to look up descriptions in pubchem. Highest priority is for a listed pubchem CID from MSFinder.  Next the full inchikey is used.  If no match is returned from pubchem, the first block (bonds) is used.  Note that this discards stereochemistry. PubChem CID, a pubchem URL, smiles (canonical), inchi and other parameters are returned.  if smiles and inchi slots are alread present (from MSFinder, for example) pubchem smiles and inchi are used to fill in missing values only - not to replace MSFinder derived values. You can create an artificial ramclustObj by creating a list with $inchikey and $cmpd (compound name) slots. 
#' 
#' @param ramclustObj ramclustR object. must contain a $inchikey slot and/or or a $MSfinder.structure slot in which to look for a pubchem CID. 
#' @param all.properties logical - if TRUE, a long list of properties is returned.  if FALSE, a short list.  FALSE is faster, TRUE more comprehensive.
#' @param get.synonyms logical - if TRUE, pubchem synonyms are returned based on pubchem inchikey. 
#' @param get.descriptions logical - if TRUE, pubchem text description is provided, with URL, if available. 
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


rc.cmpd.get.pubchem <- function(
  ramclustObj = NULL,
  all.properties = TRUE,
  get.synonyms = TRUE,
  get.descriptions = TRUE
) {

  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikey slot found, please 'annotate' first", '\n')
  }
  
  params <- c(
    "all.properties" = all.properties,
    "get.synonyms" = get.synonyms,
    "get.descriptions" = get.descriptions
  )
  
  # create list for storing values
  cmpd.props <- as.list(vector(length = length(ramclustObj$cmpd)))
  names(cmpd.props) <- ramclustObj$cmpd
  if(all.properties) {
    fp <- vector(length = length(ramclustObj$cmpd))
  }
  
  if(is.null(ramclustObj$msfinder.structure)) {
    msfinder <- FALSE
  } else {
    msfinder <- TRUE
  }
  
  
  for(i in 1:length(ramclustObj$inchikey)) {  ## change back to 1
    if(is.na(ramclustObj$inchikey[i])) {next}
    
    time.start <- Sys.time()
    
    # Sys.sleep(0.2)
    #  for(i in 1:4) {
    cat("cmpd", i, '\n')
    
    cid <- NA
    ## get CID from inchikey
    if(msfinder) {
      if(is.null(ramclustObj$msfinder.structure[[i]]$resources)) next
      if(!grepl("PubChem=", ramclustObj$msfinder.structure[[i]]$resources)) next
      cid <- unlist(strsplit(ramclustObj$msfinder.structure[[i]]$resources, ","))
      cid <- cid[grep("PubChem=", cid)]
      cid <- gsub("PubChem=", "", cid)
      cid <- unlist(strsplit(cid, ";"))[1]
    }
    
    if(is.na(cid)) {
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", 
        ramclustObj$inchikey[i], 
        "/cids/JSON")
      out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return(NULL)})
      if(is.null(out)) {
        bond.block <- unlist(strsplit(ramclustObj$inchikey[i], "-"))[1]
        link <- paste0(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", 
          bond.block, 
          "/cids/JSON")
        out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return(NULL)})
        
      }
      
      
      if(is.null(out)) next
      if(is.null(out$IdentifierList)) next
      if(is.null(out$IdentifierList$CID)) next
      cid <- out$IdentifierList$CID
      cid <- as.numeric(cid)
      cid <- cid[which.min(cid)]
    }
    
    if(all.properties) {
      ## get properties for compound - this can take several seconds per compound.  
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        cid,
        "/property/",
        "CanonicalSMILES,",
        "IsomericSMILES,",
        "InChI,",
        "InChIKey,",
        "XLogP,",
        "TPSA,",
        "Complexity,",
        "Charge,",
        "HBondDonorCount,",
        "HBondAcceptorCount,",
        "RotatableBondCount,",
        "HeavyAtomCount,",
        "AtomStereoCount,",
        "DefinedAtomStereoCount,",
        "UndefinedAtomStereoCount,",
        "BondStereoCount,",
        "DefinedBondStereoCount,",
        "UndefinedBondStereoCount,",
        "CovalentUnitCount,",
        "Volume3D,",
        "XStericQuadrupole3D,",
        "YStericQuadrupole3D,",
        "ZStericQuadrupole3D,",
        "FeatureCount3D,",
        "FeatureAcceptorCount3D,",
        "FeatureDonorCount3D,",
        "FeatureAnionCount3D,",
        "FeatureCationCount3D,",
        "FeatureRingCount3D,",
        "FeatureHydrophobeCount3D,",
        "ConformerModelRMSD3D,",
        "EffectiveRotorCount3D,",
        "ConformerCount3D,",
        "Fingerprint2D",
        "/JSON")
    } else {
      ## get properties for compound - this can take several seconds per compound.  
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        cid,
        "/property/",
        "CanonicalSMILES,",
        "IsomericSMILES,",
        "InChI,",
        "InChIKey,",
        "XLogP,",
        "TPSA,",
        "Complexity,",
        "Charge,",
        "HBondDonorCount,",
        "HBondAcceptorCount,",
        "RotatableBondCount,",
        "HeavyAtomCount",
        "/JSON")
    }
    
    
    
    out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return("NULL")})
    if(is.null(out)) next
    if(is.null(out$PropertyTable)) next
    
    prop.names <- names(out$PropertyTable$Properties)
    prop.vals <- as.vector(out$PropertyTable$Properties)
    names(prop.vals) <- prop.names
    cmpd.props[[i]] <- prop.vals
    
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
  
  
  
  # 
  # if(is.null(ramclustObj$pubchem.url)) {
  #   ramclustObj$pubchem.url <- rep(NA, nrow(pubchem))
  #   do <- which(!is.na(pubchem[,"CID"]))
  #   ramclustObj$pubchem.url[do] <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem[do,"CID"])
  # } else  {
  #   fix <- which(is.na(ramclustObj$pubchem.url) & !is.na(pubchem$pubchem.url))
  #   if(length(fix) > 0) {
  #     ramclustObj$pubchem.url[fix] <- pubchem$pubchem.url[fix]
  #   }
  # }
  # 
  # 
  # if(is.null(ramclustObj$smiles)) {ramclustObj$smiles <- pubchem$smiles} else {
  #   fix <- which(is.na(ramclustObj$smiles) & !is.na(pubchem$smiles))
  #   if(length(fix) > 0) {
  #     ramclustObj$smiles[fix] <- pubchem$smiles[fix]
  #   }
  # }
  # if(is.null(ramclustObj$inchi)) {ramclustObj$inchi <- pubchem$inchi} else {
  #   fix <- which(is.na(ramclustObj$inchi) & !is.na(pubchem$inchi))
  #   if(length(fix) > 0) {
  #     ramclustObj$inchi[fix] <- pubchem$inchi[fix]
  #   }
  # }
  # if(is.null(ramclustObj$cid)) {ramclustObj$cid <- pubchem$CID} else {
  #   fix <- which(is.na(ramclustObj$cid) & !is.na(pubchem$CID))
  #   if(length(fix) > 0) {
  #     ramclustObj$cid[fix] <- pubchem$CID[fix]
  #   }
  # }
  # 
  # if(!is.null(ramclustObj$inchi) & !is.null(ramclustObj$smiles)) {
  #   get.inchi <- which(!is.na(ramclustObj$smiles) & is.na(ramclustObj$inchi))
  #   if(length(get.inchi) > 0) {
  #     for(i in get.inchi) {
  #       tmp <- unlist(webchem::cs_convert(ramclustObj$smiles[i], from = "smiles", to = "inchi"))
  #       if(length(tmp) == 1) {
  #         if(grepl("InChI=", tmp)) {
  #           ramclustObj$inchi[i] <- tmp
  #         }
  #       } else {
  #         if(grepl("InChI=", tmp)) {
  #           ramclustObj$inchi[i] <- tmp[1]
  #         }
  #       }
  #       rm(tmp)
  #     }
  #   }
  # }
  # 
  # ramclustObj$history$pubchem<- paste(
  #   "Smiles structures were retrieved for each inchikey without a structure using the Pubchem API (Djoumbou 2016) called from RAMClustR using the getSmilesInchi function."
  # )
  
  tmp <- data.frame(cmpd.props[[1]], stringsAsFactors = FALSE)
  for(i in 1:length(ramclustObj$inchikey)) {
    n <- names(cmpd.props[[i]])
    v <- cmpd.props[[i]]
    tmp[i,n] <- v
  }
  
  if(all.properties) {
    fp <- tmp$Fingerprint2D
    tmp <- tmp[, (1:ncol(tmp))[-which(names(tmp) == "Fingerprint2D")]]
    ramclustObj$pubchem.fp    <- fp
  }
  
  st <-tmp[,1:5]
  tmp <- tmp[,-(1:5)]
  
  ramclustObj$pubchem.cmpds <- st
  ramclustObj$pubchem.props <- tmp
  
  
  if(get.synonyms) {
    cat('retreiving synonyms', '\n')
    syns <- as.list(rep(NA, length(ramclustObj$ann)))
    for(i in 1:length(syns)) {
      if(is.na(st$InChIKey[i])) next
      
      time.start <- Sys.time()
      
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", 
        paste(st$InChIKey[i], collapse = ","), 
        "/synonyms/JSON")
      out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return(NULL)})
      if(is.null(out$InformationList$Information)) next
      syns[[i]] <- unlist(out$InformationList$Information$Synonym)
      
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
    ramclustObj$pubchem.synonyms <- syns
  }
  
  ## this may be improved by direcly referencing other websites such as CHEBI or HMDB.
  if(get.descriptions) {
    cat('retreiving descriptions', '\n')
    desc <- as.list(rep(NA, length(ramclustObj$ann)))
    for(i in 1:length(desc)) {
      
      time.start <- Sys.time()
      
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", 
        paste(st$InChIKey[i], collapse = ","), 
        "/description/JSON")
      out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return(NULL)})
      if(is.null(out$InformationList$Information$Description)) next
      tmp.desc <- out$InformationList$Information$Description
      use <- which(!is.na(out$InformationList$Information$Description))[1]
      if(length(use) == 0) next
      desc[[i]] <- out$InformationList$Information$Description[use]
      if(is.null(out$InformationList$Information$DescriptionURL)) next
      desc[[i]] <- paste(desc[[i]], out$InformationList$Information$DescriptionURL[use])
      
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
    ramclustObj$pubchem.description <- desc
  }
  
  if(is.null) {ramclustObj$params <- list()}
  ramclustObj$params$rc.cmpd.get.pubchem <- params
  
  ramclustObj$history$pubchem <- paste0(
    "Pubchem data was retrieved using the PUG rest interface (Kim 2019)."
  )
  
  return(ramclustObj)
  
}


