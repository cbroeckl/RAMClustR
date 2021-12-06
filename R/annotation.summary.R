#' annotation.summary()
#'
#' Write a .csv file containing a summary of the annotations in the ramclustR object. 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param outfile file path/name of output csv summary file.  if NULL (default) will be exported to spectra/annotaionSummary.csv
#' @details this function exports a csv file summarizing annotation evidence for each compound
#' @return nothing
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @concept ramlclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept xcms
#' @concept MSFinder
#' @author Corey Broeckling
#' @export 


annotation.summary<-function(ramclustObj = NULL,
                             outfile = NULL
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(!is.null(outfile)) {
    f<-basename(outfile)
    p<-dirname(outfile)
    if(!dir.exists(p)) {
      dir.create(p)
    }
  } else {
    outfile <- paste0(getwd(), "/spectra/annotationSummary.csv")
    if(!dir.exists('spectra')) {dir.create('spectra')}
  }
  
  out<- data.frame("cmpd" = ramclustObj$cmpd,
                   "rt" = ramclustObj$clrt,
                   "annotation" = ramclustObj$ann,
                   "ann.confidence" = ramclustObj$annconf,
                   "median signal" = as.vector(apply(ramclustObj$SpecAbund, 2, "median", na.rm = TRUE))) 
  
  if(any(names(ramclustObj) == "cmpd.use")) {
    out<- data.frame(out, "qc.cv.acceptable" = ramclustObj$cmpd.use)
  }
  if(any(names(ramclustObj) == "qc.cv.cmpd.full")) {
    out<- data.frame(out, "qc.cv" = ramclustObj$qc.cv.cmpd.full)
  }
  
  if(any(names(ramclustObj) == "M")) {
    out<- data.frame(out, "inferred M" = ramclustObj$M)
  }
  if(any(names(ramclustObj) == "zmax")) {
    out<- data.frame(out, "zmax" = ramclustObj$zmax)
  }
  if(any(names(ramclustObj) == "formula")) {
    out<- data.frame(out, "molecular formula" = ramclustObj$msfinder.formula)
  } else {
    if(any(names(ramclustObj) == "msfinder.formula")) {
      out<- data.frame(out, "MSFinder inferred formula" = ramclustObj$msfinder.formula)
    }
  }
  if(any(names(ramclustObj) == "pubchem.cid")) {
    out<- data.frame(out, "pubchem.cid" = ramclustObj$pubchem.cid)
  }
  if(any(names(ramclustObj) == "inchikey")) {
    out<- data.frame(out, "inchikey" = ramclustObj$inchikey)
  }
  if(any(names(ramclustObj) == "inchi")) {
    if(any(names(ramclustObj) == "pubchem")) {
      if((any(names(ramclustObj$pubchem) == "properties"))) {
        if((any(names(ramclustObj$pubchem$properties) == "InChI"))) {
          r.inch <- which(is.na(ramclustObj$inchi))
          if(length(r.inch) > 0) {
            ramclustObj$inchi[r.inch] <- ramclustObj$pubchem$properties[r.inch,"InChI"]
          }
        }
      }
    }
    out<- data.frame(out, "inchi" = ramclustObj$inchi)
  }
  if(any(names(ramclustObj) == "hmdb.url")) {
    # out<- data.frame(out, "hmdb.url" = paste0("<a href='", ramclustObj$hmdb.url, "'>", ramclustObj$hmdb.name, "</a>")) 
    out<- data.frame(out, "hmdb.url" = ramclustObj$hmdb.url)
  }
  if(any(names(ramclustObj) == "lm.url")) {
    out<- data.frame(out, "lm.url" = ramclustObj$lm.url)
  }
  if(any(names(ramclustObj) == "pubchem.url")) {
    if(any(names(ramclustObj) == "pubchem")) {
      if((any(names(ramclustObj$pubchem) == "pubchem"))) {
        if((any(names(ramclustObj$pubchem$properties) == "pubchem.url"))) {
          r.pcurl <- which(is.na(ramclustObj$pubchem.url))
          if(length(r.pcurl) > 0) {
            ramclustObj$inchi[r.pcurl] <- ramclustObj$pubchem$properties[r.pcurl,"pubchem.url"]
          }
        }
      }
    }
    out<- data.frame(out, "pubchem.url" = ramclustObj$pubchem.url)
  }
  if(any(names(ramclustObj) == "chebi.url")) {
    out<- data.frame(out, "chebi.url" = ramclustObj$chebi.url)
  }
  if(any(names(ramclustObj) == "synonyms")) {
    if(any(names(ramclustObj) == "pubchem")) {
      if((any(names(ramclustObj$pubchem) == "synonyms"))) {
        r.syn <- which(is.na(ramclustObj$synonyms))
        if(length(r.syn) > 0) {
          for(z in r.syn) {
            ramclustObj$synonyms[[z]] <- ramclustObj$pubchem$synonyms[[z]]
          }
        }
      }
    }
    out<- data.frame(out, "synonyms" = sapply(1:length(ramclustObj$synonyms), 
                                              FUN = function(x) {
                                                paste(ramclustObj$synonyms[[x]], collapse = " __ ")
                                              }))
  }
  if(any(names(ramclustObj) == "classyfire")) {
    out<- data.frame(out, ramclustObj$classyfire)
  }
  
  write.csv(out, file = outfile, row.names = FALSE)
  
  # ## library(xlsx)
  # out.wb <- createWorkbook()
  # sheet <- createSheet(out.wb, sheetName = "annotation summary")
  # addDataFrame(out, sheet)
  # rows   <- createRow(sheet, 1)
  # cells <- createCell(rows, 1:ncol(out))
  # 
  # url.cols <- grep("url", names(out))
  # for(i in 1:nrow(out)) {
  #   for(j in url.cols) {
  #     if(is.na(out[i, j])) {next}
  #   }
  #   addHyperlink(cell = cells[[i,j]], address = as.character(out[i,j]))
  # }
  # saveWorkbook(out.wb, file = "ann.summ.xlsx")
}

