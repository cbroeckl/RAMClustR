#' get.synonyms()
#'
#' Use chemical translation service to retreive synonyms for ramclustR inchikeys
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param get.db logical - should the HMDB, LipidMaps, CheBI, and Pubchem names be retreived when available?
#' @param lipid.short.hand - logical - should the stringr pacakge be used to look for lipid short hand nomenclature in the synonyms? Only used when update.names = TRUE
#' @param update.names logical - should the ramclustObj$ann slot (annotation) be updated based on the above? selection of which synonym to chose is difficult to automate well - new name may not be the most commonly used. 
#' @details this function uses the chemical translation service (http://cts.fiehnlab.ucdavis.edu/), HMDB, LipidMaps, and PubChem databases to retreive synonymns and compound names ror a given inchikey).  Lipid shorthand (i.e. PC(36:6)) can be identified and used when available.  Precendence for naming is lipid.short.hand > HMDB > LipidMaps > Pubchem > original assignment. 
#' @return an updated ramclustR object, with the RC$synonyms slot containing a list with length equal to the number of compounds in the dataset.  Each list element represents a character vector of all the synonyms returns (or NA, for compounds with no inchikey)  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references http://cts.fiehnlab.ucdavis.edu/static/download/CTS2-MS2015.pdf 
#' @importFrom jsonlite fromJSON
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms', 'MSFinder', 'chemical translation service', 'cts'
#' @author Corey Broeckling
#' @export 


get.synonyms <- function(ramclustObj = NULL,
                         get.db = TRUE,
                         update.names = TRUE,
                         lipid.short.hand = TRUE
) {
  
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikeys found in this ramclustR object")
  }
  
  synonyms<-as.list(rep(NA, length(ramclustObj$ann)))
  names(synonyms)<-ramclustObj$cmpd
  
  if(any(names(ramclustObj) == "msfinder.structure")) {
    for(x in 1:length(ramclustObj$cmpd)) {
      if(is.data.frame(ramclustObj$msfinder.structure[[x]])) {
        res <- ramclustObj$msfinder.structure[[x]][1,"resources"]
        res <- unlist(strsplit(res, ","))
        if(length(res) > 0) {
          if(length(synonyms[[x]]) > 1) {
            synonyms[[x]] <- c(synonyms[[x]], res)
          }
          if(length(synonyms[[x]]) == 1 ) {
            if(is.na(synonyms[[x]])) {
              synonyms[[x]] <- res
            } else {
              synonyms[[x]] <- c(synonyms[[x]], res)
            }
          }
        }
      }
    }
  }
  
  
  if(get.db) {
    cat("referencing web content: please be patient", '\n')
    
    hmdb.name <- rep(NA, length(ramclustObj$cmpd))
    hmdb.url <- rep(NA, length(ramclustObj$cmpd))
    for(x in 1:length(ramclustObj$cmpd)) {
      # for(x in x:length(ramclustObj$cmpd)) {
      tmp <- synonyms[[x]]
      if(any(grep("HMDB=", tmp))) {
        tmp <- tmp[grep("HMDB=", tmp)[1]]
        tmp <- unlist(strsplit(tmp, ";"))[1]
        tmp <- gsub("HMDB=", "", tmp)
      } else {next}
      hmdb.url[x] <- paste0("http://www.hmdb.ca/metabolites/", tmp)
      link <- paste0("http://www.hmdb.ca/metabolites/", tmp, ".xml")
      # tryCatch(suppressWarnings(out<-readLines(link)), error = function(x) {NA}, finally = NA)
      out<- tryCatch(suppressWarnings(xml2::read_xml(link)), error = function(x) {return(NA)})
      if(!is.na(out)) {
        baz <- xml2::xml_find_all(out, ".//name")
        hmdb.name[x] <- xml2::xml_text(baz, "name")[1]
      } else {hmdb.url[x] <- NA}
    }
    ramclustObj$hmdb.name <- hmdb.name
    ramclustObj$hmdb.url <- hmdb.url
    
    lm.name <- rep(NA, length(ramclustObj$cmpd))
    lm.url <- rep(NA, length(ramclustObj$cmpd))
    for(x in 1:length(ramclustObj$cmpd)) {
      tmp <- synonyms[[x]]
      if(any(grep("LipidMAPS=", tmp))) {
        tmp <- tmp[grep("LipidMAPS=", tmp)[1]]
        tmp <- unlist(strsplit(tmp, ";"))[1]
        tmp <- gsub("LipidMAPS=", "", tmp)
      } else {next}
      link <- paste0("https://www.lipidmaps.org/rest/compound/lm_id/", tmp, "/name/")
      lm.url[x] <- paste0("https://www.lipidmaps.org/data/LMSDRecord.php?LMID=", tmp)
      # chemnames <- unlist(fromJSON(out))
      #tryCatch(read_xml(link), error = function(x) {return(NA)})
      tryCatch(out<- suppressWarnings(readLines(link)), error = function(x) {return(NA)})
      if(!is.na(out)){
        chemnames <- unlist(jsonlite::fromJSON(out))
        if(any(names(chemnames) == "name")) lm.name[x] <- chemnames["name"]
      }
    }
    ramclustObj$lm.name <- lm.name
    ramclustObj$lm.url <- lm.url
    
    pubchem.name <- rep(NA, length(ramclustObj$cmpd))
    pubchem.url <- rep(NA, length(ramclustObj$cmpd))
    for(x in 1:length(ramclustObj$cmpd)) {
      tmp <- synonyms[[x]]
      if(any(grep("PubChem=", tmp))) {
        tmp <- tmp[grep("PubChem=", tmp)[1]]
        tmp <- unlist(strsplit(tmp, ";"))[1]
        tmp <- gsub("PubChem=", "", tmp)
      } else {next}
      pubchem.url[x] <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", tmp, "#section=Top")
      link <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", tmp, "/description/JSON")
      out<- tryCatch(suppressWarnings(readLines(link)), error = function(x) {return("stuff.and.nonsense")})
      if(out[1] != "stuff.and.nonsense") {
        out <- jsonlite::fromJSON(out)$InformationList$Information[,"Title"]
        out <- out[-which(is.na(out))]
        if(length(out) > 0) {
          pubchem.name[x] <- out[1]
        }else {pubchem.url[x] <- NA}
      }
    }
    ramclustObj$pubchem.name <- pubchem.name
    if(any(!is.na(pubchem.url))) {ramclustObj$pubchem.url <- pubchem.url}
    
    
    chebi.url <- rep(NA, length(ramclustObj$cmpd))
    for(x in 1:length(ramclustObj$cmpd)) {
      tmp <- synonyms[[x]]
      if(any(grep("ChEBI=", tmp))) {
        tmp <- tmp[grep("ChEBI=", tmp)[1]]
        tmp <- unlist(strsplit(tmp, ";"))[1]
        tmp <- gsub("ChEBI=", "", tmp)
      } else {next}
      chebi.url[x] <- paste0("https://www.ebi.ac.uk/chebi/searchId.do?chebiId=", tmp)
    }
    ramclustObj$chebi.url <- chebi.url
    
  }
  
  for(i in 1:length(ramclustObj$ann)) {
    if(!is.na(ramclustObj$inchikey[i])) {
      link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/synonyms/", ramclustObj$inchikey[i])
      # link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/synonyms/", 
      #                unlist(strsplit(ramclustObj$inchikey[i], "-"))[1]
      #                )
      #tryCatch(read_xml(link), error = function(x) {return(NA)})
      tryCatch(suppressWarnings(out<-readLines(link)), error = function(x) {return(NA)})
      syns<-unlist(jsonlite::fromJSON(out))
      
      link <- paste0("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/Chemical%20Name/", ramclustObj$inchikey[i])
      # link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/synonyms/", 
      #                unlist(strsplit(ramclustObj$inchikey[i], "-"))[1]
      #                )
      tryCatch(suppressWarnings(out<-readLines(link)), 
               error = function(x) {return(NA)})
      chemnames <- unlist(jsonlite::fromJSON(out))
      if(any(grepl("result", names(chemnames)))) {
        chemnames <- as.vector(chemnames[grepl("result", names(chemnames))])
        syns <- c(chemnames, syns)
      }
      
      
      
      if(!is.null(syns)) {
        syns <- unique(c(ramclustObj$ann[i], syns))
        syns <- syns[order(nchar(syns))]
        if(any(names(ramclustObj) == "msfinder.structure")) {   
          if(is.data.frame(ramclustObj$msfinder.structure[[i]])){
            if(nrow(ramclustObj$msfinder.structure[[i]]) == 1) {
              res <- ramclustObj$msfinder.structure[[i]][1,"resources"]
              res <- unlist(strsplit(res, ","))
              if(length(synonyms[[i]]) > 1) {
                syns <- c(synonyms[[i]], res)
              }
            }
          }  
          synonyms[[i]] <- syns
        }
      }
    }
    
  }
  
  if(lipid.short.hand){
    lmname <- rep(NA, length(ramclustObj$cmpd))
    for(x in 1:length(lmname)) {
      if(length(synonyms[[x]]) <= 1) next
      tars <- synonyms[[x]][stringr::str_detect(synonyms[[x]], "\\([0-9]{1,2}\\:[0-9]{1,2}\\)")]
      if(length(tars) == 0) next
      nc <- nchar(tars)
      lmname[x] <- tars[which.min(nc)]
    }
    
    ramclustObj$short.name <- lmname
  }
  
  if(update.names) {
    if(any(names(ramclustObj) == "short.name")) {
      new.name <- ramclustObj$short.name
    } else {new.name <- rep(NA, length(ramclustObj$ann))}
    if(any(names(ramclustObj) == "hmdb.name")) {
      do <- which(is.na(new.name))
      new.name[do] <- hmdb.name[do]
    }
    if(any(names(ramclustObj) == "lm.name")) {
      do <- which(is.na(new.name))
      new.name[do] <- lm.name[do]
    }
    if(any(names(ramclustObj) == "pubchem.name")) {
      do <- which(is.na(new.name))
      new.name[do] <- lm.name[do]
    }
    
    if(any(is.na(new.name))) {
      do <- which(is.na(new.name))
      new.name[do] <- ramclustObj$ann[do]
    }
    
    ramclustObj$ann <- new.name
  }
  
  ramclustObj$history <- paste(ramclustObj$history, 
                               "Synonyms were retreived using the RAMClustR get.synonyms function within RAMClustR, calling the Chemical Translation Service (Wohlgemuth 2010) as well as database API for Pubchem (Kim 2019), LipidMaps (Fahy 2007), and HMDB (Wishart 2018), when appropriate.")
  ramclustObj$synonyms <- synonyms
  return(ramclustObj)
}



