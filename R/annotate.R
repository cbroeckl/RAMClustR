#' evalute ramSearch, MSFinder mssearch, MSFinder Structure, MSFinder Formula, and findmain output to annotate spectra of ramclustR object
#'
#' After running RAMSearch (msp) and MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param msfinder.dir full path to MSFinder directory - used for naming refinement
#' @param standardize.names logical: if TRUE, use inchikey for standardized chemical name lookup (http://cts.fiehnlab.ucdavis.edu/)
#' @details this function imports the output from the MSFinder program to annotate the ramclustR object
#' @return an updated ramclustR object, with the RC$ann and RC$ann.conf slots updated to annotated based on output from 1. ramsearch output, 2. msfinder mssearch, 3. msfinder predicted structure, 4. msfinder predicted formula, and 5. interpretMSSpectrum inferred molecular weight, with listed order as priority.  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @references http://cts.fiehnlab.ucdavis.edu/static/download/CTS2-MS2015.pdf 
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms', 'MSFinder'
#' @author Corey Broeckling
#' @export 


annotate<-function(ramclustObj = RC,
                   msfinder.dir = "K:/software/MSFinder/MS-FINDER program ver. 2.20",
                   standardize.names = TRUE,
                   delay.time = 0
) {
  
  if(!dir.exists(msfinder.dir)) {
    stop("msfinder directory does not exist: please set 'msfinder.dir' option as your full msfinder directory path")
  }
  
  use.short.inchikey = TRUE
  
  sfile<-list.files(paste0(msfinder.dir, "/Resources"), pattern = "ExistStructureDB_vs")
  if(length(sfile)==0) {
    stop("no structure DB file found in msfinder.dir / Resources")
  }
  
  if(length(sfile) > 1) {
    sfile <- sfile[which.max(as.numeric(gsub(".esd", "", gsub("ExistStructureDB_vs", "", sfile))))]
  }
  
  d<-read.delim2(paste0(msfinder.dir, "/Resources/", sfile), header = TRUE, na.strings = "N/A", quote = "", stringsAsFactors = FALSE)
  
  if(any(names(ramclustObj)=="M")) {
    findmain = TRUE
  } else {findmain = FALSE}
  
  if(any(names(ramclustObj) == "msfinder.formula")) {
    formula = TRUE
  } else {formula = FALSE}
  
  if(any(names(ramclustObj) == "msfinder.structure")) {
    structure = TRUE
  } else {structure = FALSE}
  
  if(any(names(ramclustObj) == "msfinder.mssearch.details")) {
    mssearch = TRUE
  } else {mssearch = FALSE}
  
  if(!any(names(ramclustObj) == "inchikey")) {ramclustObj$inchikey <- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "inchi")) {ramclustObj$inchi<- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "smiles"))  {ramclustObj$smiles <- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "dbid")) {ramclustObj$dbid<- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "synonyms")) {
    ramclustObj$synonyms <- as.list(rep(NA, length(ramclustObj$cmpd)))
  }
  
  
  if(mssearch) {
    for(i in 1:length(ramclustObj$ann)) {
      if((nrow(ramclustObj$msfinder.mssearch.details[[i]]$summary)>0) & (ramclustObj$cmpd[i] == ramclustObj$ann[i]))  {
        if(ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"totalscore"] >=4 ) {
          ramclustObj$inchikey[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"inchikey"]
          ramclustObj$smiles[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"smiles"]
          ramclustObj$ann[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"name"]
          ramclustObj$annconf[i]<-2
          ramclustObj$dbid[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"resources"]
        }
      }
    }
  }
  
  if(structure) {
    for(i in 1:length(ramclustObj$ann)) {
      if( is.data.frame(ramclustObj$msfinder.structure[[i]]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
        
        tmpinch<-ramclustObj$msfinder.structure[[i]][1, "inchikey"]
        tmpinch.short<-unlist(strsplit(tmpinch, "-"))[1]
        if(use.short.inchikey) {
          drow<-grep(tmpinch.short, d[,"InChIKey"])
        } else {
          drow<-grep(tmpinch, d[,"InChIKey"])
        }
        # d[drow,"InChIKey"]
        # tmp<- ramclustObj$msfinder.structure[[i]][1, "inchikey"]
        
        ramclustObj$inchikey[i]<-ramclustObj$msfinder.structure[[i]][1,"inchikey"]
        ramclustObj$smiles[i]<-ramclustObj$msfinder.structure[[i]][1,"smiles"]
        ramclustObj$annconf[i]<-2
        ramclustObj$dbid[i]<-ramclustObj$msfinder.structure[[i]][1,"resources"]
        
        if(length(drow) == 0) {
          ramclustObj$ann[i]<-ramclustObj$msfinder.structure[[i]][1, "name"]
        }
        
        if(length(drow) == 1) {
          ramclustObj$ann[i]<-d[drow, "Title"]
        }
        
        if(length(drow) > 1) {
          n<-d[drow, "Title"]
          nl<-nchar(n)
          ramclustObj$ann[i]<-n[which.min(nl)]
          ramclustObj$annconf[i]<-2
        }
        
        
      }
    }
  }
  
  if(formula) {
    for(i in 1:length(ramclustObj$ann)) {
      if( !is.na(ramclustObj$msfinder.formula[[i]]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
        ramclustObj$ann[i]<-ramclustObj$msfinder.formula[i]
        ramclustObj$annconf[i]<-3
        ramclustObj$dbid[i]<-ramclustObj$msfinder.formula.details[[i]][1,"resourcenames"]
      }
    }
  }
  
  if(findmain) {
    for(i in 1:length(ramclustObj$ann)) {
      if( !is.na(ramclustObj$M[i]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
        ramclustObj$ann[i]<-paste("M =", ramclustObj$M[i])
      }
    }
  }
  
  ramclustObj$inchikey[which(ramclustObj$inchikey == "undefined")]<-NA
  
  if(standardize.names) {
    cat("using chemical translation service - requires interet access and may take a few minutes to complete", '\n')
    
    inchikey2inchi<-which(!is.na(ramclustObj$inchikey) & is.na(ramclustObj$inchi))
    for(i in inchikey2inchi) {
      Sys.sleep(delay.time)
      if(!is.na(ramclustObj$inchikey[i])) {
        
        link <- paste0("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/InChI Code/", ramclustObj$inchikey[i])
        out<-NA
        start<-Sys.time()
        while(is.na(out[1])) {
          tryCatch(suppressWarnings(out<-readLines(link)), error = function(x) {NA}, finally = NA)
          if(as.numeric(Sys.time() - start) > 5) {
            ramclustObj$inchi[[i]] <- NA
            break
          }
        }
        inchis<-unlist(fromJSON(out)$result)
        if(length(inchis) == 0) {
          ramclustObj$inchi[i] <- NA
        }
        if(length(inchis)>=1) {
          ramclustObj$inchi[i] <- inchis[1]
        }
      }
    }
    
    inchi2smiles<-which(!is.na(ramclustObj$inchi) & is.na(ramclustObj$smiles))
    if(length(inchi2smiles) > 0) {
      for(i in inchi2smiles) {
        inchi<-ramclustObj$inchi[i]
        m<-parse.inchi(inchi)[[1]]
        s<-get.smiles(m)
        rm(m)
        ramclustObj$smiles[i]<-s
      }
    }
    
    for(i in 1:length(ramclustObj$ann)) {
      Sys.sleep(delay.time)
      if(!is.na(ramclustObj$inchikey[i])) {
        
        link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/synonyms/", ramclustObj$inchikey[i])
        out<-NA
        start<-Sys.time()
        while(is.na(out[1])) {
          tryCatch(suppressWarnings(out<-readLines(link)), error = function(x) {NA}, finally = NA)
          stop<-Sys.time()
          if(as.numeric(stop - start) > 5) {
            ramclustObj$synonyms[[i]] <- NA
            break
          }
        }
        syns<-unlist(fromJSON(out))
        if(length(syns) == 0) {
          ramclustObj$synonyms[[i]] <- NA
        }
        if(length(syns)>=1) {
          nc<-nchar(syns)
          syns<-syns[order(nc, decreasing = FALSE)]
          ramclustObj$synonyms[[i]] <- syns
        }
      #}}
        
        if(is.na(ramclustObj$inchi[i])) {          
          link <- paste0("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/InChI Code/", ramclustObj$inchikey[i])
          out<-NA
          while(is.na(out[1])) {
            tryCatch(suppressWarnings(out<-readLines(link)), error = function(x) {NA}, finally = NA)
          }
          inchi<-as.character(unlist(fromJSON(out))["result"])
          if(length(inchi) == 0) {
            ramclustObj$inchi[[i]] <- NA
          }
          if(length(inchi)>=1) {
            nc<-nchar(inchi)
            inchi<-inchi[order(nc, decreasing = FALSE)]
            ramclustObj$inchi[i] <- inchi[1]
          }
        }
      }
    }
  }
  
  
  
  ## modify compound names to make them unique
  nt <- table(ramclustObj$ann)
  ramclustObj$inchikey[which(is.na(ramclustObj$inchikey))] <- "NA"
  while(any(nt > 1)) {
    do<-which(nt>1)[1]
    mtch<-which(ramclustObj$ann == names(nt)[do])
    if(length(unique(ramclustObj$inchikey[mtch]))==1){
      ramclustObj$ann[mtch] <- paste(ramclustObj$ann[mtch], c(1:length(mtch)), sep = "__")
    } else {
      for(j in 1:length(mtch)) {
        if(any(names(ramclustObj)=="synonyms")) {
          cur<-which(ramclustObj$synonyms[[mtch[j]]] == names(do))
          if(length(cur)>0) {
            if(length(ramclustObj$synonyms[[mtch[j]]]) > cur){
              ramclustObj$ann[mtch[j]] <- ramclustObj$synonyms[[mtch[j]]][cur+1]
            } else {ramclustObj$ann[mtch[j]] <- paste(ramclustObj$ann[mtch[j]], j, sep = "__")}
          } else {ramclustObj$ann[mtch[j]] <- paste(ramclustObj$ann[mtch[j]], j, sep = "__")}
        }
      }
    }
    nt <- table(ramclustObj$ann)
  }
  
  ramclustObj$inchikey[which(ramclustObj$inchikey == "NA")] <- NA
  
  
  
  return(ramclustObj)
}

