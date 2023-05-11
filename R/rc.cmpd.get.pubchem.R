#' rc.cmpd.get.pubchem
#'
#' use pubchem rest and view APIs to retrieve structures, CIDs (if a name or inchikey is given), synonyms, and optionally vendor data, when available. 
#' @details useful for moving from chemical name to digital structure representation. greek letters are assumed to be 'UTF-8' encoded, and are converted to latin text before searching.   if you are reading in your compound name list, do so with 'encoding' set to 'UTF-8'. 
#' @param ramclustObj RAMClust Object input.  if used, ramclustObj$CID, ramclustObj$inchikey, and ramclustObj$ann are used as input, in that order, and ramclustObj is returned with $pubchem slot appended.
#' @param cmpd.names character vector.  i.e. c("caffeine", "theobromine", "glucose")
#' @param cmpd.cid numeric integer vector.  i.e. c(2519, 5429, 107526)
#' @param cmpd.inchikey character vector.  i.e. c("RYYVLZVUVIJVGH-UHFFFAOYSA-N", "YAPQBXQYLJRXSA-UHFFFAOYSA-N", "GZCGUPFRVQAUEE-SLPGGIOYSA-N")
#' @param cmpd.smiles character vector.  i.e. c("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CN1C=NC2=C1C(=O)NC(=O)N2C")
#' @param use.parent.cid logical.  If TRUE, the CID for each supplied name/inchikey is used to retrieve its parent CID (i.e. the parent of sodium palmitate is palmitic acid).  The parent CID is used to retrieve all other names, properties.
#' @param manual.entry logical.  if TRUE, user input is enabled for compounds not matched by name. A browser window will open with the pubchem search results in your default browser. 
#' @param get.vendors logical.  if TRUE, vendor data is returned for each compound with a matched CID.  Includes vendor count and vendor product URL, if available
#' @param priority.vendors charachter vector.  i.e. c("MyFavoriteCompany", "MySecondFavoriteCompany").  If these vendors are found, the URL returned is from priority vendors. Priority is given by order input by user. 
#' @param get.properties logical.  if TRUE, physicochemical property data are returned for each compound with a matched CID.
#' @param get.synonyms = TRUE. logical.  if TRUE, retrieve pubchem synonyms.  returned to $synonyms slot
#' @param find.short.lipid.name = TRUE. logical.  If TRUE, and get.synonyms = TRUE, looks for lipid short hand names in synonyms list (i.e. PC(36:6)). returned to $short.name slot.  Short names are assigned only if assign.short.names = TRUE.
#' @param find.short.synonym = TRUE. logical.  If TRUE, and get.synonyms = TRUE, looks for lipid short synonyms, with prioritization for names with fewer numeric characters (i.e. database accession numbers or CAS numbers). returned to $short.name slot.  Short names are assigned only if assign.short.names = TRUE.
#' @param max.name.length = 20.  integer.  If names are longer than this value, short names will be searched for, else, retain original name.
#' @param assign.short.name = TRUE.  If TRUE, short names from find.short.lipid.name and/or find.short.synonym = TRUE, short names are assigned the be the default annotation name ($ann slot), and original annotations are moved to $long.name slot.
#' @param all.props logical.  If TRUE, all pubchem properties (https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865567) are returned.  If false, only a subset (faster).
#' @param get.bioassays logical. If TRUE, return a table summarizing existing bioassay data for that CID. 
#' @param get.pathways logical.  If TRUE, return a table of metabolic pathways for that CID.
#' @param write.csv logical.  If TRUE, write csv files of all returned pubchem data. 
#' @param search.name character.  optional name to assign to pubchem search to name output .csv files.   
#' @return returns a list with one or more of $pubchem (compound name and identifiers) - one row in dataframe per CID; $properties contains physicochemical properties - one row in dataframe per CID; $vendors contains the number of vendors for a given compound and selects a vendor based on 'priority.vendors' supplied, or randomly choses a vendor with a HTML link - one row in dataframe per CID;  $bioassays contains a summary of bioassay activity data from pubchem - zero to many rows in dataframe per CID
#' @author Corey Broeckling
#' 
#' @export 
#' 

rc.cmpd.get.pubchem <- function(
    ramclustObj = NULL,
    search.name = NULL,
    cmpd.names = NULL,
    cmpd.cid = NULL,
    cmpd.inchikey = NULL,
    cmpd.smiles = NULL,
    use.parent.cid = FALSE,
    manual.entry = FALSE,
    get.vendors = FALSE,
    priority.vendors = c("Sigma Aldrich", "Alfa Chemistry", "Acros Organics", "VWR", 
                         "Alfa Aesar", "molport", "Key Organics", "BLD Pharm"),
    get.properties = TRUE,
    all.props = FALSE,
    get.synonyms = TRUE,
    find.short.lipid.name = TRUE,
    find.short.synonym = TRUE,
    max.name.length = 30,
    assign.short.name = TRUE,
    get.bioassays = TRUE,
    get.pathways = TRUE,
    write.csv = TRUE
    
) {
  
  ## function to close failed pubchem queries to prevent 
  ## all connections are in use error
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
  
  
  ## test connection to pubchem servers
  html <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/58-08-2/property/inchikey/JSON"
  out <- tryCatch(
    {
      jsonlite::read_json(html)
    },
    error=function(cond) {
      closePubchemConnections()
      stop("pubchem rest connection could not be established. This may be due to:", '\n',
           "  -  lack of internet access", '\n',
           "  -  pubchem server is down", '\n',
           "  -  pubchem server has blocked access for this IP address (try restarting your R and Rstudio session)",
           '\n')
    },
    warning=function(cond) {
      warning("pubchem rest triggered a warning.", '\n')
    },
    finally={
      cons <- suppressWarnings(showConnections(all = TRUE)); rm(cons)
    }
  )
  rm(out)
  
  if(!is.null(ramclustObj)) {
    cmpd.names <- ramclustObj$ann
    if(!is.null(ramclustObj$inchikey)) {
      cmpd.inchikey <- ramclustObj$inchikey
    }
  }
  
  
  ## check that lengths of cmpd.* make sense, and standardize
  l <- c("cmpd.names" = length(cmpd.names), 
         "cmpd.cid" = length(cmpd.cid), 
         "cmpd.inchikey" = length(cmpd.inchikey),
         "cmpd.smiles" = length(cmpd.smiles)
  )
  if(l["cmpd.names"] != max(l) & l["cmpd.names"] > 0) {
    stop(
      paste(" Length cmpd.names = ", l['cmpd.names'], "while length", names(l)[which.max(l)], "=", max(l), '\n',
            " - these must be either exactly the same length or one must be NULL", sep = " ")
    )
  }
  
  if(l["cmpd.cid"] != max(l) & l["cmpd.cid"] > 0) {
    stop(
      paste(" Length cmpd.cid = ", l['cmpd.cid'], "while length", names(l)[which.max(l)], "=", max(l), '\n',
            " - these must be either exactly the same length or one must be NULL", sep = " ")
    )
  }
  
  if(l["cmpd.inchikey"] != max(l) & l["cmpd.inchikey"] > 0) {
    stop(
      paste(" Length cmpd.inchikey = ", l['cmpd.inchikey'], "while length", names(l)[which.max(l)], "=", max(l), '\n',
            " - these must be either exactly the same length or one must be NULL", sep = " ")
    )
  }
  
  if(l["cmpd.smiles"] != max(l) & l["cmpd.smiles"] > 0) {
    stop(
      paste(" Length cmpd.smiles = ", l['cmpd.smiles'], "while length", names(l)[which.max(l)], "=", max(l), '\n',
            " - these must be either exactly the same length or one must be NULL", sep = " ")
    )
  }
  
  if(l["cmpd.names"] < max(l)) cmpd.names <- rep(NA, max(l))
  if(l["cmpd.cid"] < max(l)) cmpd.cid <- rep(NA, max(l))
  if(l["cmpd.inchikey"] < max(l)) cmpd.inchikey <- rep(NA, max(l))
  if(l["cmpd.smiles"] < max(l)) cmpd.smiles <- rep(NA, max(l))
  
  ## store original data in data.frame
  d <- data.frame("user.cmpd" = cmpd.names, 
                  "user.cid" = cmpd.cid,
                  "user.inchikey" = cmpd.inchikey,
                  "use.smiles" = cmpd.smiles)
  
  pubchem <- list()
  
  ## clean up text
  cmpd.names <- trimws(cmpd.names)
  cmpd.names[which(nchar(cmpd.names) < 1)] <- NA
  cmpd.names <- gsub(" ", "%20", cmpd.names)
  cmpd.inchikey <- trimws(cmpd.inchikey)
  cmpd.inchikey[which(nchar(cmpd.inchikey) < 1)] <- NA
  cmpd.cid <- trimws(cmpd.cid)
  cmpd.cid[which(nchar(cmpd.cid) < 1)] <- NA
  cmpd.smiles <- trimws(cmpd.smiles)
  cmpd.smiles[which(nchar(cmpd.smiles) < 1)] <- NA
  
  
  greek <- read.csv(paste(find.package("RAMClustR"), "/params/greek.csv", sep=""), header=TRUE, encoding = "UTF-8", stringsAsFactors = FALSE)
  
  
  
  for(i in 1:nrow(greek)) {
    cmpd.names <- gsub(greek[i,2], greek[i,1], cmpd.names)
  }  
  
  if(!all(is.na(cmpd.smiles))) {
    html.smiles <- cmpd.smiles
    html.translation <- rbind(
      c("[", "%5B"),
      c("]", "%5D"),
      c("@", "%40"),
      c("/", "%2F"),
      c("\\", "%5C"),
      c("=", "%3D")
    )
    
    for(i in 1:length(html.smiles)) {
      if(is.na(html.smiles[i])) {next}
      html.smiles[i] <- as.character(html.smiles[i])
      html.smiles[i] <- trimws(unlist(strsplit(html.smiles[i], "|", fixed = TRUE))[1])
      for(j in 1:nrow(html.translation)) {
        html.smiles[i] <- gsub(html.translation[j,1], html.translation[j,2], html.smiles[i], fixed = TRUE)
      }
    }
    
  } else {
    html.smiles <- cmpd.smiles
  }
  
  
  
  ## get missing CIDs from inchikeys first
  ## if more than one inchikey per compound, lowest value CID is used
  do <- which(is.na(cmpd.cid) & !is.na(cmpd.inchikey))
  if(length(do) > 0) {
    cat("getting cid from inchikey", '\n')
    do <- cmpd.inchikey[do]
    do.l <- split(do, ceiling(seq_along(do)/1))
    for(i in 1:length(do.l)) {
      keep <- which(!do.l[[i]]=="NA")
      # cat(do.l[[i]][keep], '\n')
      Sys.sleep(0.2)
      if(length(keep)==0) next
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",
                     paste0(do.l[[i]][keep], collapse = ","),
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
      tmp <- tmp[order(tmp[,"CID"], decreasing = TRUE),]
      
      for(x in 1:nrow(tmp)) {
        rp <- which(cmpd.inchikey == tmp[x,"InChIKey"])[1]
        if(!is.na(rp)) {cmpd.cid[rp] <- tmp[x, "CID"]}
      }
    }
  }
  
  
  ## get missing CIDs from smiles next
  do.ind <- which(is.na(cmpd.cid) & !is.na(cmpd.smiles))
  do <- html.smiles[do.ind]
  if(length(do) > 0) {
    for(i in 1:length(do)) {
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", do[i], "/cids/TXT")
      tmp <- tryCatch(
        {
          as.integer(as.character(readLines(html)))
        },
        error=function(cond) {
          closePubchemConnections()
          stop("pubchem rest connection could not be established. This may be due to:", '\n',
               "  -  lack of internet access", '\n',
               "  -  pubchem server is down", '\n',
               "  -  pubchem server has blocked access for this IP address (try restarting your R and Rstudio session)",
               '\n')
        },
        warning=function(x) {
          closePubchemConnections()
          NA
        },
        finally={
          closePubchemConnections()
          cons <- suppressWarnings(showConnections(all = TRUE)); rm(cons)
        }
        
      )
      if(!is.na(tmp)) cmpd.cid[do.ind[i]] <- tmp 
    }
  }
  
  ## get missing CIDs from names next
  do.ind <- which(is.na(cmpd.cid) & !is.na(cmpd.names))
  do <- cmpd.names[do.ind]
  # cat(do)
  if(length(do) > 0) {
    cat("getting cid from names", '\n')
    for(i in 1:length(do)) {
      Sys.sleep(0.2)
      # cat(do[i])
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                     do[i],
                     "/property/", "inchikey", "/JSON")
      # html <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/8007-46-3/property/inchikey/JSON"
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
      tmp <- tmp[order(tmp[,"CID"], decreasing = TRUE),]
      cmpd.cid[do.ind[i]] <- tmp[1, "CID"]
      # cat(i, tmp[1, "CID"], do[i], '\n')
    }
  }
  
  
  ##
  cid <- cmpd.cid
  
  ### offer opportunity to revise CIDs with manual input, if need be. use name only (inchikey should be unambiguous)
  if(manual.entry){
    missing <- which(is.na(cid) & !is.na(cmpd.names)) 
    
    cat(length(missing), "missing CIDs", '\n',
        "  enter '1' to skip these compounds", '\n', 
        "  enter '2' to enter CIDs manually", '\n')
    readback <- readline()
    
    if(readback == '2') {
      for(i in missing) {
        cat("please enter CID number or hit 'enter' to skip to next (no CID found)", '\n',
            "If you wish to quit manual entry, enter 'q' to return current Output only.", '\n',
            cmpd.names[i], '\n')
        Sys.sleep(0.2)
        utils::browseURL(paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", cmpd.names[i]))
        readback <- readline()
        if(readback == "q") {break}
        cmpd.cid[i] <- readback
        
      }
    }
  }
  
  cmpd.cid[which((cmpd.cid == "NA"))] <- NA
  cid <- cmpd.cid
  d <- data.frame(d, "cmpd.cid" = cmpd.cid, stringsAsFactors = FALSE)
  
  ## find.parent.cid, if TRUE
  if(use.parent.cid) {
    cat("getting parent cid from cid", '\n')
    parent.cid <- cmpd.cid
    # do.ind <- which(!is.na(cmpd.cid) & !is.na(cmpd.names))
    do.ind <- which(!is.na(cmpd.cid))
    for(i in do.ind) {
      Sys.sleep(0.2)
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     cmpd.cid[i],
                     "/cids/TXT?cids_type=parent")
      out <- tryCatch(
        {
          readLines(html)
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
        }
      )
      if(is.na(out[1])) next
      out <- sort(as.numeric(out))
      parent.cid[i] <- out[1]
      rm(out)
    }
    cid <- parent.cid
    d <- data.frame(d, "parent.cid" = parent.cid, stringsAsFactors = FALSE)
  }
  
  d <- data.frame(d, 'cid' = cid, stringsAsFactors = FALSE)
  #  pubchem$compounds <- d
  
  cid.l <- split(cid, ceiling(seq_along(cid)/1))
  
  ## pubchem URL
  do <- which(!is.na(d$cid))
  urls <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", d$cid[do])
  d$pubchem.url <- rep(NA, nrow(d))
  d$pubchem.url[do] <- urls
  
  ## get pubchem name
  pubchem.name <- rep(NA, length(cid))
  cat("getting pubchem compound name from cid", '\n')
  for(i in 1:length(cid.l)) {
    keep <- which(!cid.l[[i]]=="NA")
    if(length(keep) == 0) next
    # urls <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", d$CID[do])
    # d$pubchem.url <- rep("", nrow(d))
    # d$pubchem.url[do] <- urls
    html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                   paste0(cid.l[[i]][keep], collapse = ","),
                   "/description/", "JSON")
    # out <- jsonlite::read_json(html)
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
      }
    )
    if(is.na(out[1])) next
    
    for(j in 1:length(out$InformationList$Information)) {
      tmp <- out$InformationList$Information[[j]]
      if(is.null(tmp$Title)) next
      tmp.title <- tmp$Title
      pubchem.name[which(cid == tmp$CID)] <- tmp.title
    }
  }
  
  d <- data.frame(d, "pubchem.name" = pubchem.name, stringsAsFactors = FALSE)
  pubchem$pubchem <- d
  
  ## Get properties
  if(get.properties) {
    cat("getting physicochemical properties and structure representations from cid", '\n')
    if(all.props) {
      props <- c(
        'MolecularFormula',
        'MolecularWeight',
        'CanonicalSMILES',
        'IsomericSMILES',
        'InChI',
        'InChIKey',
        'IUPACName',
        'XLogP',
        'ExactMass',
        'MonoisotopicMass',
        'TPSA',
        'Complexity',
        'Charge',
        'HBondDonorCount',
        'HBondAcceptorCount',
        'RotatableBondCount',
        'HeavyAtomCount',
        'IsotopeAtomCount',
        'AtomStereoCount',
        'DefinedAtomStereoCount',
        'UndefinedAtomStereoCount',
        'BondStereoCount',
        'DefinedBondStereoCount',
        'UndefinedBondStereoCount',
        'CovalentUnitCount',
        'Volume3D',
        'XStericQuadrupole3D',
        'YStericQuadrupole3D',
        'ZStericQuadrupole3D',
        'FeatureCount3D',
        'FeatureAcceptorCount3D',
        'FeatureDonorCount3D',
        'FeatureAnionCount3D',
        'FeatureCationCount3D',
        'FeatureRingCount3D',
        'FeatureHydrophobeCount3D',
        'ConformerModelRMSD3D',
        'EffectiveRotorCount3D',
        'ConformerCount3D',
        'Fingerprint2D'
      )
    } else {
      props <- c(
        'MolecularFormula',
        'MolecularWeight',
        'CanonicalSMILES',
        'IsomericSMILES',
        'InChI',
        'InChIKey',
        'XLogP',
        'ExactMass',
        'MonoisotopicMass',
        'TPSA',
        'HBondDonorCount',
        'HBondAcceptorCount'
      )
    }
    
    properties <- d[,0]
    for(i in 1:length(cid.l)) {
      Sys.sleep(0.2)
      keep <- which(!cid.l[[i]]=="NA")
      if(length(keep) == 0) next
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(cid.l[[i]][keep], collapse = ","),
                     "/property/", paste0(props, collapse =","), "/JSON")
      # out <- jsonlite::read_json(html)
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
        }
      )
      if(is.na(out)) next
      for(j in 1:length(out$PropertyTable$Properties)) {
        tmp <- data.frame(out$PropertyTable$Properties[[j]], stringsAsFactors = FALSE)
        properties[which(cid == tmp$CID),names(tmp)] <- tmp
      }
    }
    dimnames(properties)[[2]][1] <- "cid"
    pubchem$properties <- properties
  }
  
  if(get.vendors) {
    cat("getting vendor data from cid", '\n')
    vendors <- d[,"cid", drop = FALSE]
    do <- which(!is.na(d$cid))
    urls <- paste0(pubchem$pubchem$pubchem.url[do], "#section=Chemical-Vendors")
    n.vendors <- rep(NA, nrow(d))
    vendor.urls <- rep(NA, nrow(d))
    pubchem.vendors.url <-rep(NA, nrow(d))
    for(i in do) {
      Sys.sleep(0.2)
      out <- tryCatch(
        {
          jsonlite::read_json(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/",
                                     cid[i],
                                     "/JSON"))
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
        }
      )
      if(is.na(out)) next
      cats <- sapply(1:length(out$SourceCategories$Categories), FUN = function(x) {
        out$SourceCategories$Categories[[x]]$Category
      })
      cat.select <- which(cats == "Chemical Vendors")
      if(length(cat.select) == 0) next
      n.vendors[i] <- length(out$SourceCategories$Categories[[cat.select]]$Sources)
      
      vendor.names <- sapply(1:length(out$SourceCategories$Categories[[cat.select]]$Sources),
                             FUN = function(x) {
                               out$SourceCategories$Categories[[cat.select]]$Sources[[x]]$"SourceName"
                             }
      )
      for(j in priority.vendors) {
        use <- agrep(j, vendor.names,  max.distance = 0.2)
        if(length(use) > 0) {
          vendor.url <- out$SourceCategories$Categories[[cat.select]]$Sources[[use[1]]]$SourceRecordURL[1]
          if(length(vendor.url) == 0) vendor.url <- NA
          vendor.urls[i] <- vendor.url
          break
        }
      }
      if(is.na(vendor.urls[i]))  {
        vendor.url <- out$SourceCategories$Categories[[cat.select]]$Sources[[sample((1:length(vendor.names)), 1)]]$"SourceRecordURL"[1]
        if(length(vendor.url) == 0) vendor.url <- ""
        vendor.urls[i] <- vendor.url
      }
      pubchem.vendors.url[do] <- urls
    }
    
    vendors <- data.frame(vendors, 
                          "pubchem.vendors.url" = pubchem.vendors.url,
                          "n.vendors" = n.vendors, 
                          "vendor.url" = vendor.urls, 
                          stringsAsFactors = FALSE)
    pubchem$vendors <- vendors
    
  }
  
  if(get.synonyms) {
    cat("getting synonym data from cid", '\n')
    
    pubchem$synonyms <- as.list(rep(NA, length(ramclustObj$cmpd)))
    names(pubchem$synonyms) <- ramclustObj$cmpd
    
    if(find.short.lipid.name | find.short.synonym) {
      pubchem$short.name <- rep(NA, nrow(d))
    }
    
    for(i in 1:length(cid.l)) {
      keep <- which(!cid.l[[i]]=="NA")
      if(length(keep) == 0) next
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(cid.l[[i]][keep], collapse = ","),
                     "/synonyms/JSON")
      #out <- jsonlite::read_json(html)
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
        }
      )
      if(is.na(out)) next
      
      for(j in 1:length(out$InformationList$Information)) {
        on <- which(cid ==  cid.l[[i]][keep[j]])
        if(length(on)==0) next
        if(length(on) > 1) {
          if(!(any(ls() == "multi.cid"))) {
            multi.cid <- rep(1, 1)
            names(multi.cid)[1] <- cid.l[[i]][keep[j]]
            on <- on[1]
          } else {
            if(any(names(multi.cid) == cid.l[[i]][keep[j]])) {
              which.multi <- which(names(multi.cid) == cid.l[[i]][keep[j]])
              multi.cid[which.multi] <- multi.cid[which.multi] + 1
              on <- on[multi.cid[which.multi]]
            } else {
              multi.cid <- c(multi.cid, 1)
              names(multi.cid)[length(multi.cid)] <- cid.l[[i]][keep[j]]
              on <- on[1]
            }
          }
        }
        tmp <- out$InformationList$Information[[j]]
        if(is.null(tmp$Synonym)) next
        syns <- unlist(tmp$Synonym)
        pubchem$synonyms[[on]] <- syns
        
        if(find.short.lipid.name) {
          if(length(pubchem$synonyms[[on]]) > 1) {
            tars <- pubchem$synonyms[[on]][which(
              stringr::str_detect(pubchem$synonyms[[on]], "\\([0-9]{1,2}\\:[0-9]{1,2}\\)") && (nchar(pubchem$synonyms[[on]]) <= max.name.length)
            )]
            if(length(tars) > 0) {
              nc <- nchar(tars)
              pubchem$short.name[on] <- tars[which.min(nc)]
            }
          }
        }
        
        if(find.short.synonym & is.na(pubchem$short.name[on])) {
          if(length(pubchem$synonyms[[on]]) > 1) {
            letter.rich <- nchar(pubchem$synonyms[[on]])
            nchars <- nchar(pubchem$synonyms[[on]])
            pattern <- "[[:digit:]]+"
            numbers <- sapply(1:length(pubchem$synonyms[[on]]), FUN = function(k) {
              nchar(paste(unlist(regmatches(pubchem$synonyms[[on]][k], 
                                            gregexpr(pattern, pubchem$synonyms[[on]][k]))), 
                          collapse = "", sep = ""))
            })
            letter.rich <- 1- (numbers/nchars)
            use <- which(nchars < max.name.length & letter.rich > 0.6)
            if(length(use) == 0) {
              use <- which(nchars < max.name.length & letter.rich > 0.3)
            }
            if(length(use) == 0) {
              use <- which(nchars < max.name.length)
            }
            if(length(use) == 0) {
              use <- which.min(nchars)
            }
            if(length(use) == 0) next
            use <- use[1]
            pubchem$short.name[on] <- pubchem$synonyms[[on]][use]
          }
        }
      }
      
      
      
      if(assign.short.name & !is.null(ramclustObj)) {
        ramclustObj$original.ann <- ramclustObj$ann
        ramclustObj$ann[which(!(is.na(pubchem$short.name)))] <- pubchem$short.name[which(!(is.na(pubchem$short.name)))]
      }
    }
  }
  
  
  if(get.bioassays) {
    cat("getting bioasssay from cid", '\n')
    for(i in 1:length(cid.l)) {
      keep <- which(!cid.l[[i]]=="NA")
      if(length(keep) == 0) next
      url <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        paste(cid.l[[i]][keep], collapse = ","),
        "/assaysummary/CSV"
      )
      
      bd <- tryCatch(
        {
          read.csv(url)
        },
        error=function(cond) {
          closePubchemConnections()
          return( data.frame("cid" = rep(NA, 0)))
        },
        warning=function(cond) {
          closePubchemConnections()
          return( data.frame("cid" = rep(NA, 0)))
        },
        finally={
          closePubchemConnections()
        }
      )
      if(nrow(bd)==0) next
      if(any(ls()=="bioassays")) {
        # cat('registered TRUE')
        bioassays <- rbind(bioassays, bd)
        
      } else {
        # cat('registered FALSE')
        bioassays <- bd       
      }
    }
    if(any(ls()=="bioassays")) {
      dimnames(bioassays)[[2]][which(dimnames(bioassays)[[2]] == "CID")] <- "cid"
      pubchem$bioassays <- bioassays
    } else {
      bioassays <- data.frame("cid" = rep(NA, 0))
      pubchem$bioassays <- bioassays
    }
  }
  
  if(get.pathways) {
    cat("getting pathway from cid", '\n')

    for(i in 1:length(cid)) {
      pathway <- get.pathways(cid = cid[i])
      pathway <- cbind("cid" = rep(cid[i]), pathway)
      if(!any(ls()=="pathways")) pathways <- pathway[0,]
      if(is.na(pathway[1,'name'])) next
      pathways <- rbind(pathways, pathway)
      
    }
    pubchem$pathways <- pathways
  }
  
  for(i in 1:length(pubchem)) {
    if(!is.data.frame(pubchem[[i]])) next
    if(nrow(pubchem[[i]]) == length(ramclustObj$cmpd)) {
      pubchem[[i]] <- data.frame("cmpd" = ramclustObj$cmpd, pubchem[[i]])
    }
  }
  
  if(write.csv) {
    if(is.null(search.name)) {search.name = "pubchem.data"}
    dir.create(search.name)
    write.csv(pubchem[[1]], file = paste0(search.name,"/",names(pubchem)[1], ".csv"))
    if(length(pubchem)>1) {
      for(i in 2:length(pubchem)) {
        if(is.data.frame(pubchem[[i]])) {
          write.csv(pubchem[[i]], 
                    file = paste0(search.name,"/",names(pubchem)[i], ".csv"),
                    row.names = FALSE)
        }
      }
    }
  }
  
  if(!is.null(ramclustObj)) {
    ramclustObj$pubchem <- pubchem
    return(ramclustObj)
  } else {
    return(pubchem)
  }
  
}
