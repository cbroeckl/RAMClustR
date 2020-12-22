#' get.pubchem.data
#'
#' use pubchem rest and view APIs to retreive structures, CIDs (if a name or inchikey is given), and optionally vendor data, when available. 
#' @details useful for moving from chemical name to digital structure represtation. greek letters are assumed to be 'UTF-8' encoded, and are converted to latin text before searching.   if you are reading in your compound name list, do so with 'encoding' set to 'UTF-8'. 
#' @param cmpd.names character vector.  i.e. c("caffeine", "theobromine", "glucose")
#' @param cmpd.cid numeric integer vector.  i.e. c(2519, 5429, 107526)
#' @param cmpd.inchikey character vector.  i.e. c("RYYVLZVUVIJVGH-UHFFFAOYSA-N", "YAPQBXQYLJRXSA-UHFFFAOYSA-N", "GZCGUPFRVQAUEE-SLPGGIOYSA-N")
#' @param use.parent.cid logical.  If TRUE, the CID for each supplied name/inchikey is used to retreive its parent CID (i.e. the parent of sodium palmitate is palmitic acid).  The parent CID is used to retrieve all other names, properties.
#' @param manual.entry logical.  if TRUE, user input is enabled for compounds not matched by name. A browser window will open with the pubchem search results in your default browser. 
#' @param get.vendors logical.  if TRUE, vendor data is returned for each compound with a matched CID.  Includes vendor count and vendor product URL, if available
#' @param priority.vendors charachter vector.  i.e. c("MyFavoriteCompany", "MySecondFavoriteCompany").  If these vendors are found, the URL returned is from priority vendors. Priority is given by order input by user. 
#' @param get.properties logical.  if TRUE, physicochemical property data are returned for each compound with a matched CID.
#' @param all.props logical.  If TRUE, all pubchem properties (https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865567) are returned.  If false, only a subset (faster).
#' @param get.bioassays logical. If TRUE, return a table summarizing existing bioassay data for that CID. 

#' @return returns a list with one or more of $puchem (compound name and identifiers) - one row in dataframe per CID; $properties contains pysicochemical properties - one row in dataframe per CID; $vendors contains the number of vendors for a given compound and selects a vendor based on 'priortity.vendors' supplied, or randomly choses a vendor with a HTML link - one row in dataframe per CID;  $bioassays contains a summary of bioassay activity data from pubchem - zero to many rows in dataframe per CID
#' @author Corey Broeckling
#' 
#' @export 
#' 

get.pubchem.data <- function(
  cmpd.names = NULL,
  cmpd.cid = NULL,
  cmpd.inchikey = NULL,
  use.parent.cid = TRUE,
  manual.entry = FALSE,
  get.vendors = TRUE,
  priority.vendors = c("Sigma Aldrich", "Alfa Chemistry", "Acros Organics", "VWR", 
                       "Alfa Aesar", "molport", "Key Organics", "BLD Pharm"),
  get.properties = TRUE,
  all.props = TRUE,
  get.bioassays = TRUE
  
) {
  
  
  ## test connection to pubchem servers
  html <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/58-08-2/property/inchikey/JSON"
  out <- tryCatch(
    {
      jsonlite::read_json(html)
    },
    error=function(cond) {
      stop("pubchem rest connection could not be established. This may be due to:", '\n',
           "  -  lack of internet acces", '\n',
           "  -  puchem server is down", '\n',
           "  -  pubchem server has blocked access for this IP address (try restarting your R and Rstudio session)",
           '\n')
    },
    warning=function(cond) {
      warning("pubchem rest triggered a warning.", '\n')
    },
    finally={
      
    }
  )
  rm(out)
  
  ## check that lengths of cmpd.* make sense, and standardize
  l <- c("cmpd.names" = length(cmpd.names), 
         "cmpd.cid" = length(cmpd.cid), 
         "cmpd.inchikey" = length(cmpd.inchikey)
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
  
  if(l["cmpd.names"] < max(l)) cmpd.names <- rep(NA, max(l))
  if(l["cmpd.cid"] < max(l)) cmpd.cid <- rep(NA, max(l))
  if(l["cmpd.inchikey"] < max(l)) cmpd.inchikey <- rep(NA, max(l))
  
  ## store orignal data in data.frame
  d <- data.frame("user.cmpd" = cmpd.names, 
                  "user.cid" = cmpd.cid,
                  "user.inchikey" = cmpd.inchikey)
  pubchem <- list()
  
  ## clean up text
  cmpd.names <- trimws(cmpd.names)
  cmpd.names[which(nchar(cmpd.names) < 1)] <- NA
  cmpd.names <- gsub(" ", "-", cmpd.names)
  cmpd.inchikey <- trimws(cmpd.inchikey)
  cmpd.inchikey[which(nchar(cmpd.inchikey) < 1)] <- NA
  cmpd.cid <- trimws(cmpd.cid)
  cmpd.cid[which(nchar(cmpd.cid) < 1)] <- NA
  
  
  greek <- read.csv(paste(find.package("RAMClustR"), "/params/greek.csv", sep=""), header=TRUE, encoding = "UTF-8", stringsAsFactors = FALSE)
  
  for(i in 1:nrow(greek)) {
    cmpd.names <- gsub(greek[i,2], greek[i,1], cmpd.names)
  }  
  
  
  
  ## get missing CIDs from inchikeys first
  ## if more than one inchikey per compound, lowest value CID is used
  do <- which(is.na(cmpd.cid) & !is.na(cmpd.inchikey))
  if(length(do) > 0) {
    cat("getting cid from inchikey", '\n')
    do <- cmpd.inchikey[do]
    do.l <- split(do, ceiling(seq_along(do)/50))
    for(i in 1:length(do.l)) {
      Sys.sleep(0.5)
      keep <- which(!do.l[[i]]=="NA")
      if(length(keep)==0) next
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",
                     paste0(do.l[[i]][keep], collapse = ","),
                     "/property/", "inchikey", "/JSON")
      out <- tryCatch(
        {
          jsonlite::read_json(html)
        },
        error=function(cond) {
          return(NA)
        },
        warning=function(cond) {
          return(NA)
        },
        finally={
          
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
  
  ## get missing CIDs from names next
  do.ind <- which(is.na(cmpd.cid) & !is.na(cmpd.names))
  do <- cmpd.names[do.ind]
  if(length(do) > 0) {
    cat("getting cid from names", '\n')
    for(i in 1:length(do)) {
      Sys.sleep(0.5)
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                     do[i],
                     "/property/", "inchikey", "/JSON")
      # html <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/8007-46-3/property/inchikey/JSON"
      out <- tryCatch(
        {
          jsonlite::read_json(html)
        },
        error=function(cond) {
          return(NA)
        },
        warning=function(cond) {
          return(NA)
        },
        finally={
        
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
        cat(cmpd.names[i], '\n')
        cat("please enter CID number or hit 'enter' to skip to next (no CID found)", '\n',
            "If you wish to quit manual entry, enter 'q' to return current ouput only.", '\n')
        Sys.sleep(0.2)
        browseURL(paste0("https://www.ncbi.nlm.nih.gov/pccompound/?term=", cmpd.names[i]))
        readback <- readline()
        if(readback == "q") {break}
        cid[i] <- readback
        
      }
    }
  }
  cid[which((cid == "NA"))] <= NA
  d <- data.frame(d, "cmpd.cid" = cmpd.cid, stringsAsFactors = FALSE)
  
  ## find.parent.cid, if TRUE
  if(use.parent.cid) {
    cat("getting parent cid from cid", '\n')
    parent.cid <- cmpd.cid
    do.ind <- which(is.na(cmpd.cid) & !is.na(cmpd.names))
    for(i in 1:length(do.ind)) {
      Sys.sleep(0.25)
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     cmpd.cid[i],
                     "/cids/TXT?cids_type=parent")
      out <- tryCatch(
        {
          readLines(html)
        },
        error=function(cond) {
          return(NA)
        },
        warning=function(cond) {
          return(NA)
        },
        finally={
          
        }
      )
      if(is.na(out[1])) next
      out <- sort(as.numeric(out))
      parent.cid[do.ind[i]] <- out[1]
    }
    cid <- parent.cid
    d <- data.frame(d, "parent.cid" = parent.cid, stringsAsFactors = FALSE)
  }
  
  d <- data.frame(d, 'cid' = cid, stringsAsFactors = FALSE)
  #  pubchem$compounds <- d
  
  cid.l <- split(cid, ceiling(seq_along(cid)/50))
  
  ## pubchem URL
  do <- which(!is.na(d$cid))
  urls <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", d$cid[do])
  d$pubchem.url <- rep("", nrow(d))
  d$pubchem.url[do] <- urls
  
  ## get pubchem name
  pubchem.name <- rep(NA, length(cid))
  cat("getting parent pubchem compound name from cid", '\n')
  for(i in 1:length(cid.l)) {
    
    keep <- which(!cid.l[[i]]=="NA")
    urls <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", d$CID[do])
    d$pubchem.url <- rep("", nrow(d))
    d$pubchem.url[do] <- urls
    html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                   paste0(cid.l[[i]][keep], collapse = ","),
                   "/description/", "JSON")
    out <- jsonlite::read_json(html)
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
      keep <- which(!cid.l[[i]]=="NA")
      
      html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(cid.l[[i]][keep], collapse = ","),
                     "/property/", paste0(props, collapse =","), "/JSON")
      out <- jsonlite::read_json(html)
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
    urls <- paste0(urls[do], "#section=Chemical-Vendors")
    n.vendors <- rep(0, nrow(d))
    vendor.urls <- rep(NA, nrow(d))
    for(i in do) {
      Sys.sleep(0.25)
      out <- tryCatch(
        {
          jsonlite::read_json(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/",
                                     cid[i],
                                     "/JSON"))
        },
        error=function(cond) {
          # message(paste("URL invalid:", cmpd.names[i]))
          return(NA)
        },
        warning=function(cond) {
          # message(cmpd.names[i], " failed; " )
          return(NA)
        },
        finally={
         
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
          if(length(vendor.url) == 0) vendor.url <- ""
          vendor.urls[i] <- vendor.url
          break
        }
      }
      if(is.na(vendor.urls[i]))  {
        vendor.url <- out$SourceCategories$Categories[[cat.select]]$Sources[[sample((1:length(vendor.names)), 1)]]$"SourceRecordURL"[1]
        if(length(vendor.url) == 0) vendor.url <- ""
        vendor.urls[i] <- vendor.url
      }
      vendors$pubchem.vendors.url[do] <- urls
    }
    
    vendors <- data.frame(vendors, "n.vendors" = n.vendors, 
                          "vendor.url" = vendor.urls, 
                          stringsAsFactors = FALSE)
    pubchem$vendors <- vendors
    
  }
  
  if(get.bioassays) {
    cat("getting bioasssay from cid", '\n')
    for(i in 1:length(cid.l)) {
      keep <- which(!cid.l[[i]]=="NA")
      
      url <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        paste(cid.l[[i]], collapse = ","),
        "/assaysummary/CSV"
      )
      
      bd <- tryCatch(
        {
          read.csv(url)
        },
        error=function(cond) {
          return( data.frame("cid" = rep(NA, 0)))
        },
        warning=function(cond) {
          return( data.frame("cid" = rep(NA, 0)))
        },
        finally={
          
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
  
  return(pubchem)
}
