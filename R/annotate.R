#' evalute ramSearch, MSFinder mssearch, MSFinder Structure, MSFinder Formula, and findmain output to annotate spectra of ramclustR object
#'
#' After running MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param msfinder.dir full path to MSFinder directory - used for naming refinement
#' @param standardize.names logical: if TRUE, use inchikey to name lookup (http://cts.fiehnlab.ucdavis.edu/)
#' @details this function imports the output from the MSFinder program to annotate the ramclustR object
#' @return an updated ramclustR object, with new slots at $msfinder.mssearch.details and $msfinder.mssearch.scores
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms', 'MSFinder'
#' @author Corey Broeckling
#' @export 


annotate<-function(ramclustObj = RC,
                            msfinder.dir = "K:/software/MSFinder/MS-FINDER program ver. 2.20",
                            standardize.names = TRUE
                            ) {
  
  if(!dir.exists(msfinder.dir)) {
    stop("msfinder directory does not exist: please set 'msfinder.dir' option as your full msfinder directory path")
  }
  
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
  }
  
  if(any(names(ramclustObj) == "msfinder.formula")) {
    formula = TRUE
  }
  
  if(any(names(ramclustObj) == "msfinder.structure")) {
    structure = TRUE
  }
  
  if(any(names(ramclustObj) == "msfinder.mssearch.details")) {
    mssearch = TRUE
  }
  
  ramclustObj$inchikey <- rep(NA, length(ramclustObj$cmpd))
  ramclustObj$smiles <- rep(NA, length(ramclustObj$cmpd))
  ramclustObj$dbid<- rep(NA, length(ramclustObj$cmpd))
  
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
        
        drow<-grep(ramclustObj$msfinder.structure[[i]][1, "inchikey"], d[,"InChIKey"])
        
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
  
  
  if(standardize.names) {
    cat("using chemical translation service - requires interet access and may take a few minutes to complete", '\n')
    require(jsonlite)
    for(i in 1:length(ramclustObj$ann)) {
      if(!is.na(ramclustObj$inchikey[i])) {
        link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/convert/InChIKey/Chemical%20Name/", ramclustObj$inchikey[i])
        suppressWarnings(out<-readLines(link))
        names<-unlist(fromJSON(out)$result)
        if(length(names) == 0) {
          names<-ramclustObj$ann[i]
        }
        if(length(names)>1) {
          nc<-nchar(names)
          names<-names[which.min(nc)]
        }
        
        if(nchar(names) > 20) {
          
          link <- paste0("http://cts.fiehnlab.ucdavis.edu/service/synonyms/", ramclustObj$inchikey[i])
          suppressWarnings(out<-readLines(link))
          syns<-unlist(fromJSON(out))
          if(length(syns) == 0) {
            syns<-names
          }
          if(length(syns)>1) {
            nc<-nchar(syns)
            syns<-syns[which.min(nc)]
          }
          if(nchar(syns) < nchar(names)) {
            names<-syns
          }
        }
        ramclustObj$ann[i] <- names
      }
    }
    
  }
  
  
  return(ramclustObj)
}
  