#' foodb2msfinder
#'
#' convenience function for converting FoodDB database export format to MSFinder custom database import format. Before running this, please have downloaded .csv files from FoodDB with the appropriate Display Field Headers (see details)
#'
#' @param foodb.files default = NULL, if path is set, will read automatically.  If NULL, direcory selection by user. 
#' @param out.dir default = NULL.  Can set to exiseting directory with full path name.  If NULL, direcory selection by user.
#' @param out.name default = "FoodDB_for_MSFinder.txt".
#' @details Input file(s) should be csv formatted, with required headers of 'Name',	'Smiles',	'Inchikey',	'Chemical formula', and 'Mono mass' - case sensitive.  Output will be in tab delimited text format in directory of choice.
#' @return  Nothing is returned - output file written to directory set by 'out.dir' and name set by 'out.name'
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export 

fooddb2msfinder <- function(
  foodb.files = NULL,
  out.dir = NULL,
  out.name = "FoodDB_for_MSFinder.txt"
) {
  
  if(is.null(foodb.files)) {
    foodb.files <- utils::choose.files (caption = "Select input files from FooDB")
  }
  
  if(is.null(out.dir)) {
    out.dir <- utils::choose.dir(caption = "Select output directory for MSFinder DB format")
  }
  
  if(length(foodb.files) == 0) {
    stop ('no files chosen!')
  }
  
  required.names <- c('Name',	'Smiles',	'Inchikey',	'Chemical formula', 'Mono mass')
  
  out <- data.frame(matrix(nrow = 0, ncol= length(required.names)))
  names(out) <- required.names
  
  
  for(i in 1:length(foodb.files)) {
    d <- read.csv(foodb.files[i], header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    missing <- !(required.names %in% names(d))
    if(any(missing)) {
      stop(paste("The file", 
                 foodb.files[i], 
                 "is missing the header(s)",
                 required.names[missing], 
                 '\n', "note that headers are case sensitive"))
    }
    d[,"Inchikey"] <- gsub("InChIKey=", "", d[,"Inchikey"])
    d <- d[,required.names]
    tmp <- d
    tmp[,] <- NA
    out <- rbind(out, tmp)
    old.inchikeys <- unique(out[,"Inchikey"])
    new.inchikeys <- unique(d[,"Inchikey"])
    keep.inchikeys <- new.inchikeys[!(new.inchikeys %in% old.inchikeys)]
    for(j in 1:length(keep.inchikeys)) {
      out[which(is.na(out[,1]))[1],] <- d[grep(keep.inchikeys[j], d[,"Inchikey"])[1],]
    }
    out <- out[-which(is.na(out[,1])),]
    cat("after", 
        foodb.files[i], 
        nrow(out), 
        "unique compounds", '\n')
    rm(tmp); rm(d); gc()
  }
  
  nas <- which(is.na(out), arr.ind = TRUE)
  nas <- nas[,1]
  if(length(nas)>0) {
    cat("removing", 
        length(nas), 
        "compound(s) due to missing values in data", 
        '\n')
    out <- out[-nas,]
  }

  short.inchi <- sapply(1:nrow(out), FUN = function(x) {
    unlist(strsplit(out[x,"Inchikey"], "-"))[1]
  })
  
  msfinder <- data.frame(
    "Title" = out[, "Name"],
    "InChIKey" = out[,"Inchikey"],
    "Short InChIKey" = short.inchi,
    "PubChem CID" = rep("-", length(short.inchi)),
    "Exact mass" = out[,"Mono mass"], 
    "Formula" = out[,"Chemical formula"], 
    "smiles" = out[,"Smiles"], 
    "Database ID" = out[,"Name"], 
    check.names = FALSE, stringsAsFactors = FALSE
  )
  utils::write.table(msfinder, file = paste(out.dir, "/", out.name, sep = ""), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("output file can be found here:",
      '\n', 
      "  ", 
      paste(out.dir, "/", out.name, sep = ""))
  
}

