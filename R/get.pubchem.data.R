#' get.pubchem.data
#'
#' use pubchem rest and view APIs to retreive structures, CIDs (if a name or inchikey is given), and optionally vendor data, when available. 
#' @details useful for moving from chemical name to digital structure represtation. greek letters are assumed to be 'UTF-8' encoded, and are converted to latin text before searching.   if you are reading in your compound name list, do so with 'encoding' set to 'UTF-8'. 
#' @param cmpd.names character vector.  i.e. c("caffeine", "theobromine", "glucose")
#' @param cmpd.cid numeric integer vector.  i.e. c(2519, 5429, 107526)
#' @param cmpd.inchikey character vector.  i.e. c("RYYVLZVUVIJVGH-UHFFFAOYSA-N", "YAPQBXQYLJRXSA-UHFFFAOYSA-N", "GZCGUPFRVQAUEE-SLPGGIOYSA-N")
#' @param manual.entry logical.  if TRUE, user input is enabled for compounds not matched by name. A browser window will open with the pubchem search results in your default browser. 
#' @param vendor.data logical.  if TRUE, vendor data is returned for each compound with a matched CID.  Includes vendor count and vendor product URL, if available
#' @param priority.vendors charachter vector.  i.e. c("MyFavoriteCompany", "MySecondFavoriteCompany").  If these vendors are found, the URL returned is from priority vendors. Priority is given by order input by user. 
#' @return returns a data.frame with several columns: "compound", "CID", "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", "pubchem.url", "pubchem.vendors.url", "n.vendors", "vendor.url"
#' @author Corey Broeckling
#' @export 
get.pubchem.data <- function(
  cmpd.names = NULL,
  cmpd.cid = NULL,
  cmpd.inchikey = NULL,
  manual.entry = FALSE,
  vendor.data = TRUE,
  priority.vendors = c("Sigma Aldrich", "Alfa Chemistry", "Acros Organics", "VWR", "Alfa Aesar")
) {
  
  cmpd.names <- trimws(cmpd.names)
  cmpd.names <- cmpd.names[which(nchar(cmpd.names)>1)]
  cmpd.names <- gsub(" ", "-", cmpd.names)
  
  greek <- read.csv(paste(find.package("RAMClustR"), "/params/greek.csv", sep=""), header=TRUE, encoding = "UTF-8", stringsAsFactors = FALSE)
  
  for(i in 1:nrow(greek)) {
    cmpd.names <- gsub(greek[i,2], greek[i,1], cmpd.names)
  }  
  
  syns <- as.list(cmpd.names)
  for(i in 1:length(cmpd.names)) {
    #  for(i in (i+1):length(cas)) {
    Sys.sleep(0.6)
    
    out <- tryCatch(
      {
        jsonlite::read_json(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                   cmpd.names[i],
                                   "/synonyms/JSON"))
      },
      error=function(cond) {
        message(paste("URL invalid:", cmpd.names[i]))
        return(NA)
      },
      warning=function(cond) {
        message(cmpd.names[i], " failed; " )
        return(NA)
      },
      finally={
      }
    )    
    syns[[i]]  <- out
  }
  
  cid <- rep(NA, length(syns))
  for(x in 1:length(syns)) {
    if(is.na(syns[[x]])) next
    if(!is.na(syns[[x]]$InformationList$Information[[1]]$CID)) {
      cid[x] <- syns[[x]]$InformationList$Information[[1]]$CID
    }
  }
  
  
  if(manual.entry){
    missing <- which(is.na(cid)) 
    
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
  
  d <- data.frame('compound' = cmpd.names, 'cid' = cid, stringsAsFactors = FALSE)
  cid.u <- unique(d$cid)
  cid.l <- split(cid.u, ceiling(seq_along(cid.u)/50))
  
  for(i in 1:length(cid.l)) {
    keep <- which(!cid.l[[i]]=="NA")
    
    html <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                   paste0(cid.l[[i]][keep], collapse = ","),
                   "/property/CanonicalSMILES,IsomericSMILES,InChI,InChIKey,/JSON")
    out <- jsonlite::read_json(html)
    for(j in 1:length(out$PropertyTable$Properties)) {
      tmp <- data.frame(out$PropertyTable$Properties[[j]], stringsAsFactors = FALSE)
      d[which(d$cid == tmp$CID),names(tmp)] <- tmp
    }
  }
  
  d <- d[,-which(names(d) == 'cid')]
  
  if(vendor.data) {
    do <- which(!is.na(d$CID))
    urls <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", d$CID[do])
    d$pubchem.url <- rep("", nrow(d))
    d$pubchem.url[do] <- urls
    d$pubchem.vendors.url <- rep("", nrow(d))
    urls <- paste0(urls, "#section=Chemical-Vendors")
    d$pubchem.vendors.url[do] <- urls
    n.vendors <- rep(0, nrow(d))
    vendor.urls <- rep(NA, nrow(d))
    for(i in 1:nrow(d)) {
      if(is.na(d$CID[i])) next
      out <- tryCatch(
        {
          jsonlite::read_json(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/",
                                     d$CID[i],
                                     "/JSON"))
        },
        error=function(cond) {
          message(paste("URL invalid:", cmpd.names[i]))
          return(NA)
        },
        warning=function(cond) {
          message(cmpd.names[i], " failed; " )
          return(NA)
        },
        finally={
        }
      )
      cats <- sapply(1:length(out$SourceCategories$Categories), FUN = function(x) {
        out$SourceCategories$Categories[[x]]$Category
      })
      cat.select <- which(cats == "Chemical Vendors")
      if(length(cat.select) == 0) next
      n.vendors[i] <- length(out$SourceCategories$Categories[[cat.select]]$Sources)
      
      vendors <- sapply(1:length(out$SourceCategories$Categories[[cat.select]]$Sources),
                        FUN = function(x) {
                          out$SourceCategories$Categories[[cat.select]]$Sources[[x]]$"SourceName"
                        }
      )
      for(j in priority.vendors) {
        use <- agrep(j, vendors,  max.distance = 0.2)
        if(length(use) > 0) {
          vendor.url <- out$SourceCategories$Categories[[cat.select]]$Sources[[use[1]]]$SourceRecordURL[1]
          if(length(vendor.url) == 0) vendor.url <- ""
          vendor.urls[i] <- vendor.url
          break
        }
      }
      if(is.na(vendor.urls[i]))  {
        vendor.url <- out$SourceCategories$Categories[[cat.select]]$Sources[[sample((1:length(vendors)), 1)]]$"SourceRecordURL"[1]
        if(length(vendor.url) == 0) vendor.url <- ""
        vendor.urls[i] <- vendor.url
      }
      
    }
    
    d <- data.frame(d, "n.vendors" = n.vendors, 
                    "vendor.url" = vendor.urls, 
                    stringsAsFactors = FALSE)
    
  }

  return(d)
}
