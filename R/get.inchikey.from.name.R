#' get.inchikey.from.name
#'
#' use rpubchem and Chemmimer tools to retreive inchikey from chemical names. 
#' @details useful for moving from chemical name to digital structure represtation. greek letters are assumed to be 'UTF-8' encoded, and are converted to latin text before searching.   if you are reading in your compound name list, do so with 'encoding' set to 'UTF-8'. 
#' @param cmpd.names character vector.  i.e. c("caffeine", "theobromine", "glucose")
#' @return returns a data.frame with three columns: name (supplied in cmpd.names), CID (pubchem cid), and inchikey. one row for each compound, compounds are processed in order supplied. 
#' @author Corey Broeckling
#' @export 
get.inchikey.from.name <- function(
  cmpd.names = NULL,
  citation.weights = FALSE
) {
  
  require(rpubchem)
  require(ChemmineR)
  require(RCurl)
  require(XML)
  require(data.table)

  ## from rpubchem
  .check.cas <- function(cas)
  {
    if (is.null(cas) || is.na(cas))
      return(NA)

    ## Input: character vector of CAS RNs
    ## Output: logical vector indicating valid CAS RNs

    # Check each element of CAS vector against CAS format with regex.
    cas.format <- regexpr("\\d{2,7}-\\d\\d-\\d", cas, perl=TRUE) > 0 & !is.na(cas)

    # If format matches, do checksum validation.
    cas[cas.format] <- sapply(cas[cas.format], function(x) {
      # remove non-numeric
      x <- gsub("[^0-9]", "", x)

      # list of integers
      names(x) <- x
      xl <- lapply(strsplit(x, ""), as.integer)

      # checksum validation
      sapply(xl, function(y) {
        cas.length <- length(y)
        actual.check.digit <- y[cas.length]
        y <- y[-cas.length]
        expected.check.digit <- sum(rev(y) * seq_along(y)) %% 10L
        expected.check.digit == actual.check.digit
      })
    })

    # return TRUE if format matches and checksum validated
    ifelse(cas.format, cas, FALSE)
  }

  ## from rpubchem
  get.synonyms <- function(name, idtype = NULL, quiet=TRUE)
  {
    ## Input: character vector of compound names
    ## Output: data.frame with matched names, PubChem CIDs, synonyms and CAS flag
    ##
    ## API Documentation: https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    ##
    ## USAGE POLICY: Please note that PUG REST is not designed for very large volumes
    ## (millions) of requests. We ask that any script or application not make more
    ## than 5 requests per second, in order to avoid overloading the PubChem servers.
    ## If you have a large data set that you need to compute with, please contact us
    ## for help on optimizing your task, as there are likely more efficient ways to
    ## approach such bulk queries.

    curlHandle <- getCurlHandle()
    out <- data.frame(stringsAsFactors=FALSE)

    for (compound in name) {
      tryCatch(
        {
          field = NULL
          if (is.null(idtype)) {
            field <- "name="
            endpoint <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/synonyms/XML"
          } else if (idtype == 'inchikey') {
            field <- "inchikey="
            endpoint <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/synonyms/XML"
          } else if (idtype == 'cid') {
            field <- "cid="
            endpoint <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/XML"
          } else stop("Invalid idtype specified")

          res <- dynCurlReader()
          curlPerform(##postfields=paste0(field, compound),
            url=sprintf(endpoint, URLencode(compound)),
            ##post=1L,
            curl=curlHandle, writefunction = res$update)
          doc <- xmlInternalTreeParse(res$value())
          rootNode <- xmlName(xmlRoot(doc))
          if (rootNode == "InformationList") {
            xpathApply(doc, "//x:Information", namespaces="x", function(x) {
              cid <- xpathSApply(x, "./x:CID", namespaces="x", xmlValue)
              synonym <- xpathSApply(x, "./x:Synonym", namespaces="x", xmlValue)
              df <- data.frame(Name=compound, CID=cid, Synonym=synonym, stringsAsFactors=FALSE)
              out <<- rbindlist(list(out, df))
            })
          } else if (rootNode == "Fault") {
            fault <- xpathApply(doc, "//x:Details", namespaces="x", xmlValue)
            if (!quiet) {
              print(paste(compound, fault[[1]], sep=": "))
            }
          }
        },
        error=function(e) {
          print(e)
        },
        finally=Sys.sleep(0.2) # See usage policy.
      )
    }

    # CAS validation
    if (nrow(out) > 0)
      out$CAS <- .check.cas(out$Synonym)

    # Cleanup
    rm(curlHandle)
    gc()
    return(out)
  }

  greek <- read.csv(paste(find.package("RAMClustR"), "/params/greek.csv", sep=""), header=TRUE, encoding = "UTF-8", stringsAsFactors = FALSE)
  
  for(i in 1:nrow(greek)) {
    cmpd.names <- gsub(greek[i,2], greek[i,1], cmpd.names)
  }  
  # short.names <- grep(";", cmpd.names)
  # if(length(short.names)>0) {
  #   for(i in 1:length(short.names)) {
  #     names <- c(cmpd.names, unlist(strsplit(short.names[i], ";"))[1])
  #     if(length(unlist(strsplit(short.names[i], ";"))) == 2) {
  #       cmpd.names <- c(cmpd.names, paste(unlist(strsplit(short.names[i], ";"))[2], unlist(strsplit(short.names[i], ";"))[1]))
  #     }
  #   }
  # }

  
  
  cmpd.names <- trimws(cmpd.names)
  out <- data.frame('cmpd.name' = cmpd.names,
                    'CID' = rep(NA, length(cmpd.names)),
                    'inchikey' = rep(NA, length(cmpd.names)))

  cat("looking up pubchem CIDs", '\n')

  syns <- rpubchem::get.synonyms(name = cmpd.names, quiet = FALSE)
  if(nrow(syns) == 0) {
    return(as.data.frame(matrix(nrow = 0, ncol = 4)))
    stop()
  }
  for(i in 1:length(cmpd.names)) {
    keep <- which(syns[,"Name"] == cmpd.names[i])
    if(length(keep) == 0) next
    CIDs <- table(syns[keep,'CID'])
    out[i,"CID"] <- as.numeric(names(CIDs[which.max(CIDs)]))
  }
  library(ChemmineR)
  cat("looking up inchikey from CID", '\n')
  do <- out[!(is.na(out[,"CID"])),"CID"]
  compounds <- ChemmineR::getIds(unique(do))
  
  CIDs <- sapply(1:length(compounds), FUN = function(x) {
    as.numeric(compounds[[x]]@header["Molecule_Name"])
  })
  inchikeys <- sapply(1:length(compounds), FUN = function(x) {
    as.character(compounds[[x]]@datablock["PUBCHEM_IUPAC_INCHIKEY"], "-")
  })
  # pubchem_names <- sapply(1:length(compounds), FUN = function(x) {
  #   as.character(compounds[[x]]@datablock["PUBCHEM_IUPAC_TRADITIONAL_NAME"])
  # })
  for(i in 1:nrow(out)) {
    if(is.na(out[i,"CID"])) next
    inchikey <- grep(out[i,"CID"], CIDs)
    if(length(inchikey) == 0) next
    out[i, "inchikey"] <- inchikeys[inchikey]
  }
  
  if(citation.weights) {
    cat("counting pubmed citations...", '\n')
    cites <- vector(mode = "numeric")
    tmp.cmpd <- out[,1]
    for(x in c("[", "]", "(", ")", "--", "#")) {
      tmp.cmpd <- gsub(x, "-", tmp.cmpd, fixed = TRUE)
    }
    for(i in 1:nrow(out)) {
      cites <- c(cites, as.numeric(easyPubMed::get_pubmed_ids(
        paste('"', tmp.cmpd[i], '"', ' [All Fields]', sep =  ''))$Count))
    }
    
    out <- data.frame(out, 'citations' = cites)
    weights <- 0.5 + ((log10(cites+1)/max(log10(cites+1)))/2)
    
    out <- data.frame(out, 'citations' = cites, 'weights' = weights)
    
  }
  
  return(out)
}
