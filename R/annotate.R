#' evalute ramSearch, MSFinder mssearch, MSFinder Structure, MSFinder Formula, and findmain output to annotate spectra of ramclustR object
#'
#' After running RAMSearch (msp) and MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param msfinder.dir full path to MSFinder directory - used for naming refinement.
#' @param standardize.names logical: if TRUE, use inchikey for standardized chemical name lookup (http://cts.fiehnlab.ucdavis.edu/)
#' @param min.msms.score numerical: what is the minimum MSFinder similarity score acceptable.  default = 6.5
#' @param database.priority character.  Formula assignment prioritization based on presence in one or more (structure) databases.  Can be set to a single or multiple database names.  must match database names as they are listed in MSFinder precisily. Can also be set to 'all' (note that MSFinder reports all databases matched, not just databases in MSFinder parameters).  If any database is set, the best formula match to any of those databases is selected, rather than the best formula match overall.  If NULL, this will be set to include all selected databases (from ramclustObj$msfinder.dbs, retreieved from search output during import.msfinder.formulas(), when available) or 'all'.  
#' @param rescore.structure logical.  If TRUE, uses an internal scoring method, rather than default MSFinder method for structure scores.
#' @param taxonomy.inchi vector or data frame.  Only when rescore.structure = TRUE.  user can supply a vector of inchikeys.  If used, structures which match first block of inchikey receive an addition score value of 1.  
#' @param database.score logical. Only when rescore.structure = TRUE. default = TRUE.  Should we use the MSFinder database score in assigning structures, as MSFinder does? If false the database score is subtracted out.  The database score increases the chance of returning biologically relevent compound matches, but biases the annotation toward well described/databased compounds. 
#' @param citation.score logical. get.inchikey.from.name() returns a data frame which (optionally) includes a citation count and derived weights.  we can use these in taxonomy scoring to bias annotation toward compounds that are commonly referred to in pubmed literature. 
#' @param form.structure.scoring logical. default = TRUE.  Should a combined score using the product of the formula and structure scores be used?  If FALSE, structure score alone is used in annotations. 
#' @param find.inchikey logical.  default = TRUE. use chemical translation service to try to look up inchikey for chemical name.
#' @param reset logical.  If TRUE, removes any previously assigned annotations.  

#' @details this function imports the output from the MSFinder program to annotate the ramclustR object

#' @return an updated ramclustR object, with the at $msfinder.formula, $msfinder.formula.score,  $ann, and $ann.conf slots updated to annotated based on output from 1. ramsearch output, 2. msfinder mssearch, 3. msfinder predicted structure, 4. msfinder predicted formula, and 5. interpretMSSpectrum inferred molecular weight, with listed order as priority.  

#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @references http://cts.fiehnlab.ucdavis.edu/static/download/CTS2-MS2015.pdf 
#'
#' @concept ramlclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept xcms
#' @concept MSFinder
#' 
#' @author Corey Broeckling
#' @export 

annotate<-function(ramclustObj = NULL,
                   msfinder.dir = "C:/MSFinder/MSFINDER ver 3.24",
                   standardize.names = FALSE,
                   min.msms.score = 6.5,
                   database.priority = NULL,
                   # database.score = TRUE,
                   citation.score = TRUE, 
                   # rescore.structure = TRUE,
                   # form.structure.scoring = TRUE,
                   find.inchikey = TRUE,
                   taxonomy.inchi = NULL,
                   reset = TRUE
) {
  
  
  if(!dir.exists(msfinder.dir)) {
    stop("msfinder directory does not exist: please set 'msfinder.dir' option as your full msfinder directory path")
  }
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(reset) {
    ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$annconf <- rep(4, length(ramclustObj$cmpd))
    ramclustObj$ann <- ramclustObj$cmpd
    ramclustObj$inchikey <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$inchi <- NULL
    ramclustObj$smiles <- NULL
    ramclustObj$dbid <- NULL
    ramclustObj$synonyms <- NULL
    ramclustObj$classyfire <- NULL
    ramclustObj$rs.lib <- NULL
    ramclustObj$rs.specn <- NULL
    ramclustObj$rs.libn <- NULL
    ramclustObj$rs.mf <- NULL
    ramclustObj$rs.rmf <- NULL
    ramclustObj$rs.prob <- NULL
  }
  
  if(!is.null(taxonomy.inchi)) {
    if(is.data.frame(taxonomy.inchi)) {
      tax.df <- taxonomy.inchi
      taxonomy.inchi <- tax.df[,"inchikey"]
      if(citation.score) {
        if(!any(names(tax.df) == "weights")) {
          stop("no citation weights provided, consider using 
               'get.inchikey.from.name(, citation.weights = TRUE)'
               or setting 'citation.score = FALSE' in the 'annotate' function",  '\n')
          cat('it didnt stop??')
        }
      }
    }
  }
  
  
  
  
  if(!is.null(taxonomy.inchi)) {
    if(is.factor(taxonomy.inchi)) {
      taxonomy.inchi <- as.character(taxonomy.inchi)
    }
    taxonomy.inchi <- sapply(1:length(taxonomy.inchi), FUN = function(x) {
      as.character(unlist(strsplit(taxonomy.inchi[x], "-"))[1])
    })
  }
  
  use.short.inchikey = TRUE
  
  sfile<-list.files(paste0(msfinder.dir, "/Resources"), pattern = "ExistStructureDB_vs")
  # if(length(sfile)==0) {
  #   stop("no structure DB file found in msfinder.dir / Resources")
  # }
  # 
  if(length(sfile) > 1) {
    sfile <- sfile[which.max(as.numeric(gsub(".esd", "", gsub("ExistStructureDB_vs", "", sfile))))]
  }
  
  if(length(sfile) == 0) {
    
    sfile<-list.files(paste0(msfinder.dir, "/Resources"), pattern = "MsfinderStructureDB-VS")
    if(length(sfile)==0) {
      stop("no structure DB file found in msfinder.dir / Resources")
    }
    
    if(length(sfile) > 1) {
      r <- grep("_bin", sfile, fixed = TRUE)
      if(length(r) > 0) {
        sfile <- sfile[-r]
      }
      if(length(sfile) > 1) {
        sfile <- sfile[which.max(as.numeric(gsub(".esd", "", gsub("MsfinderStructureDB-VS", "", sfile))))] 
      }
    }
  }
  
  ## this is not stable.  vers 3.30 of MSFinder is now using a binary format. 
  reference.data<-read.delim2(paste0(msfinder.dir, "/Resources/", sfile), header = TRUE, na.strings = "N/A", quote = "", stringsAsFactors = FALSE)
  
  if(any(names(ramclustObj)=="rs.out")) {
    ramsearch = TRUE
  } else {ramsearch = FALSE}
  
  if(any(names(ramclustObj)=="M")) {
    findmain = TRUE
  } else {findmain = FALSE}
  
  if(any(names(ramclustObj) == "msfinder.formula.details")) {
    formula = TRUE
    #use.formula <- rep(FALSE, length(ramclustObj$msfinder.formula.details))
  } else {formula = FALSE}
  
  if(any(names(ramclustObj) == "msfinder.structure.details")) {
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
  
  if(!is.null(ramclustObj$msfinder.dbs) & is.null(database.priority)) {
    database.priority <- ramclustObj$msfinder.dbs
  }
  
  if(ramsearch) {
    ##these items will be filled and added to the RC object
    out <- ramclustObj$rs.out 
    if(is.null(ramclustObj$inchikey)) {
      ramclustObj$inchikey <- rep(NA, length(ramclustObj$cmpd))
    }
    
    ramclustObj$rs.spec	<-as.list(rep("", max(ramclustObj$featclus)))
    ramclustObj$rs.lib	<-rep("", max(ramclustObj$featclus))
    ramclustObj$rs.specn	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.libn	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.mf	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.rmf	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.prob	<-as.numeric(rep(-1, max(ramclustObj$featclus)))
    
    ##pull relevent line numbers for all spectra
    ##name line is first, so call that directly, range between name 1 and name 2 is the
    ##range of spectrum 1
    name<-which(regexpr('Matched Spectrum:', out)==1)
    ann<-which(regexpr('Annotation:', out)==1)
    origname<-which(regexpr('Original Name', out)==1)
    
    if(any((origname-name)!=2)) stop("please don't edit the output from ramsearch manually")
    if(any(ramclustObj$cmpd!=sub("Original Name: ", "", out[origname]))) {stop("compound names/order differ between ramclust object and ramsearch output")}  
    
    ##if length of name is not equal to length of ramclust Object 'cmpd' slot, something is wrong:
    if(length(name)/as.integer(as.character(ramclustObj$ExpDes[[2]]["MSlevs",1]))!= length(ramclustObj$cmpd)) stop("number of spectra in ramsearch output different than number of compounds in ramclust object")
    
    ##now pull relevent info our for each spectrum in output
    for(i in 1:length(ann)) {
      md<-out[name[i]:min((name[i+1]-1), length(out), na.rm=TRUE)]
      cname<-sub("Original Name: ", "", md[grep("Original Name: ", md)])
      ind<-as.numeric(sub("C", "", cname))
      mf<-ramclustObj$rs.mf[ind]
      newmf<-as.integer(sub("Match Factor / Dot Product: ", "", md[grep("Match Factor / Dot Product: ", md)]))
      if(ramclustObj$cmpd[ind] != cname) {
        stop(paste("something is amiss with compound ", i, ": the names do not match", sep=""))
      }
      if((nchar(sub("Annotation: ", "", md[2]))>0) & max(newmf, -1, na.rm=TRUE) >= mf ) {
        ramclustObj$ann[ind] 	<- as.character(sub("Annotation: ", "", md[grep("Annotation: ", md)]))
        ramclustObj$rs.specn[ind] 	<- as.character(sub("Matched Spectrum: ", "", md[grep("Matched Spectrum: ", md)]))
        ramclustObj$annconf[ind] 	<- as.integer(sub("Confidence: ", "", md[grep("Confidence: ", md)]))
        ramclustObj$annnotes[ind] <- as.character(sub("Comments: ", "", md[grep("Comments: ", md)]))
        ramclustObj$rs.lib[ind]	<- as.character(sub("Library: ", "", md[grep("Library: ", md)]))
        ramclustObj$rs.libn[ind]	<- as.integer(sub("Library Id: ", "", md[grep("Library Id: ", md)]))
        ramclustObj$rs.mf	[ind]	<- as.integer(sub("Match Factor / Dot Product: ", "", md[grep("Match Factor / Dot Product: ", md)]))
        ramclustObj$rs.rmf[ind]	<- as.integer(sub("Rev Match Factor / Rev Dot: ", "", md[grep("Rev Match Factor / Rev Dot: ", md)]))
        #ramclustObj$rs.prob[i]	<- sub("", "", md[grep("", md)])
        tmp.inchi <- as.character(gsub("InChIKey: ", "", md[grep("InChIKey: ", md)]))
        tmp.inchi <- gsub("InChIKey=", "", tmp.inchi)
        if(nchar(tmp.inchi) > 0) {ramclustObj$inchikey[i] <- tmp.inchi; rm(tmp.inchi)} else {ramclustObj$inchikey[i] <- NA}
        if(length(grep("Library Match Num Peaks:", md))==1) {
          ramclustObj$rs.spec[[ind]]	<- matrix(as.numeric(unlist(strsplit(md[(grep("Library Match Num Peaks:", md)+1)], " "))), ncol=2, byrow=TRUE)
        }
      }
    }
    ramclustObj$ann[which(nchar(ramclustObj$ann)<1)]<-ramclustObj$cmpd[which(nchar(ramclustObj$ann)<1)]
    
  }
  
  if(mssearch) {
    spec.formula.warnings <- vector(length = 0, mode = 'numeric')
    if(is.null(ramclustObj$msfinder.formula)) {
      ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$msfinder.structure.details))
    }
    
    msfinder.mssearch.score<-as.numeric(sapply(1:length(ramclustObj$msfinder.mssearch.details), FUN = function(x) {
      if(is.null(nrow(ramclustObj$msfinder.mssearch.details[[x]]$summary))) {
        NA 
      } else {
        ramclustObj$msfinder.mssearch.details[[x]]$summary[1,"totalscore"]
      }
    }
    )
    )
    
    ramclustObj$msfinder.mssearch.score<-msfinder.mssearch.score
    
    
    for(i in 1:length(ramclustObj$ann)) {
      if((nrow(ramclustObj$msfinder.mssearch.details[[i]]$summary)>0) & (ramclustObj$cmpd[i] == ramclustObj$ann[i]))  {
        if(ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"totalscore"] >= min.msms.score ) {
          ramclustObj$inchikey[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"inchikey"]
          ramclustObj$smiles[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"smiles"]
          ramclustObj$ann[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"name"]
          ramclustObj$annconf[i]<-2
          ramclustObj$dbid[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"resources"]
          tar.inchikey <- unlist(strsplit(ramclustObj$inchikey[i], "-"))[1]
          ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
          form <- reference.data[which(reference.data[,"Short.InChIKey"] == tar.inchikey)[1], "Formula"]
          ramclustObj$msfinder.formula[i] <- form
        }
      }
    }
  }
  
  if(structure) {
    
    
    tmpdb <- as.list(rep("", length(ramclustObj$msfinder.structure.details)))
    summary <- as.list(rep("", length(ramclustObj$msfinder.structure.details)))
    dbs <- vector(length = 0, mode = 'character')
    for(x in 1:length(ramclustObj$msfinder.structure.details)) {
      # cat(x, " ")
      if(!is.list(ramclustObj$msfinder.structure.details[[x]])) next
      if(length(ramclustObj$msfinder.structure.details[[x]]) > 0) {
        all <- ramclustObj$msfinder.structure.details[[x]][[1]][["structures"]][0,]
        for(i in 1:length(ramclustObj$msfinder.structure.details[[x]])) {
          all <- rbind(all, ramclustObj$msfinder.structure.details[[x]][[i]][["structures"]])
        }
      }
      
      if(nrow(all) > 0) {
        tmpdbs <- paste(all[,"resources"], collapse = ",")
        if(tmpdbs == "") next
        tmpdbs <- unlist(strsplit(tmpdbs, ","))
        tmpdbs <- unlist(strsplit(tmpdbs, ";"))
        tmpdbs <- tmpdbs[grepl("=", tmpdbs)]
        tmpdbs <- strsplit(tmpdbs, "=")
        if(length(tmpdbs)>0) {
          tmpdbs <- sapply(1:length(tmpdbs), FUN = function(x) {tmpdbs[[x]][1]})
          dbs <- unique(c(dbs, tmpdbs))
        }
        tmpdb[[x]] <- tmpdbs
      }
    }
    # rm(tmpdb)
    
    if(!(is.null(ramclustObj$msfinder.formula.dbs))) {
      r <- which(grepl("IsUserDefinedDB", ramclustObj$msfinder.formula.dbs))
      if(length(r)>0) {
        ramclustObj$msfinder.formula.dbs[r] <- "Database ID"
      }
    }
    
    suppressWarnings(
      if (is.null(database.priority)) {
        if(!is.null(ramclustObj$msfinder.formula.dbs)) {
          if(length(ramclustObj$msfinder.formula.dbs)>0) {
            database.priority <- ramclustObj$msfinder.formula.dbs
          } else {
            database.priority <- dbs
          }
        } else {
          database.priority <- dbs }
        
        if (database.priority == "all") {
          database.priority <- dbs
        }
      }
    )
    
    if(any(grepl('custom', database.priority))) {
      database.priority <- gsub("custom", "Database ID", database.priority)
    }
    
    ## use case insensitive matching to try to reconcile mismatches
    dbmatch <- database.priority %in% dbs
    fix <- which(!dbmatch)
    if(length(fix) > 0) {
      for(i in fix) {
        case.insensitive.match <- grep(database.priority[i], dbs, ignore.case = TRUE)
        if(length(case.insensitive.match) == 1) {
          dbs[case.insensitive.match] <- database.priority[i]
        }
      }
    }
    
    ## should case sensitive fail, remove mismatches
    dbmatch <- database.priority %in% dbs
    fix <- which(!dbmatch)
    if(length(fix) > 0) {
      warning("The following databases do not match and will be removed from consideration: ", '\n',
              "  ", database.priority[fix])
      database.priority <- database.priority[-fix]
      
    }
    
    
    search.dbs <- dbs
    priority.dbs <- database.priority
    
    
    ## create template
    template <- matrix(nrow = 0, ncol = 0)
    for(i in 1:length(ramclustObj$ann)) {
      db.m <- data.frame(lapply(1:length(database.priority), FUN = function(x) grepl(database.priority[x], ramclustObj$msfinder.formula.details[[i]][,"resourcenames"])))
      db.m <- apply(db.m, 1, FUN = 'any')
      if(!any(db.m))  next
      form.tab <- ramclustObj$msfinder.formula.details[[i]][db.m,]
      if(is.null(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[1]]])) next
      for(j in 1:nrow(form.tab)) {
        if(is.null(names(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]]))) next
        template.col.names <- c(paste0("f.", names(form.tab)), names(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]]))
        template <- (matrix(nrow = 0, ncol = (ncol(form.tab) + ncol(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]]))))
        if(ncol(template > 0)) break
      }
      if(ncol(template > 0)) break
    }
    
    
    
    for(i in 1:length(ramclustObj$ann)) {
      # cat(i, " ")
      if(is.na(ramclustObj$msfinder.formula[i])) {
        # keep <- vector(length = 0)
        # for(j in 1:length(ramclustObj$msfinder.structure.details[[i]])) {
        #   # if(!is.na(keep)) break
        #   for(k in 1:length(database.priority)) {
        #     # if(!is.na(keep)) break
        #     if(is.list(ramclustObj$msfinder.structure.details[[i]])) {
        #       if(any(grepl(database.priority[k], ramclustObj$msfinder.structure.details[[i]][[j]][["structures"]][,"resources"]))) {
        #         keep <- c(keep, j)
        #         # if(!is.na(keep)) break
        #       }
        #     }
        #   }
        # }
        
        db.m <- data.frame(lapply(1:length(database.priority), FUN = function(x) grepl(database.priority[x], ramclustObj$msfinder.formula.details[[i]][,"resourcenames"])))
        db.m <- apply(db.m, 1, FUN = 'any')
        if(!any(db.m)) next
        form.tab <- ramclustObj$msfinder.formula.details[[i]][db.m,]
        tmp <- template
        
        # tmp <- data.frame(tmp)
        for(j in 1:nrow(form.tab)) {
          if(!any(names(ramclustObj$msfinder.structure.details[[i]]) == form.tab$name[j])) next
          if(nrow(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]]) == 0) next
          for(k in 1:nrow(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]])) {
            if(is.null(dimnames(tmp))) {
              dimnames(tmp)[[2]] <- c(paste0("f.", names(form.tab)), names(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]]))
            }
            nr <- as.vector(c(unlist(form.tab[j,]), unlist(ramclustObj$msfinder.structure.details[[i]][[form.tab$name[j]]][["structures"]][k,]))) 
            names(nr) <- NULL
            tmp <- rbind(tmp, nr)
          }
        }
        
        dimnames(tmp)[[2]] <- template.col.names
        tmp <- data.frame(tmp, stringsAsFactors = FALSE)
        for(x in c("f.totalscore", "totalhrrulesscore", "totalbondcleavagescore", "totalmassaccuracyscore", 
                   "totalfragmentlinkagescore", "totalbonddissociationenergyscore", 
                   "databasescore", "substructureassignmentscore", "rtsimilarityscore", 
                   "risimilarityscore", "totalscore")) {
          tmp[,x] <- as.numeric(tmp[,x])
        }
        
        if(nrow(tmp)>0) {
          for(j in c("rtsimilarityscore", "risimilarityscore")) {
            if(max(tmp[,j]) <=0) {
              tmp[,j] <- 0
            }
          }
        }
        
        use <- c("totalhrrulesscore", "totalbondcleavagescore", "totalmassaccuracyscore", 
                 "totalfragmentlinkagescore", "totalbonddissociationenergyscore", 
                 "databasescore", "substructureassignmentscore", "rtsimilarityscore", 
                 "risimilarityscore", "f.totalscore") 
        
        # if(database.score) {
        #   use <- use[use!="databasescore"]
        # } 
        
        
        # if(rescore.structure) {
        #   tmp$totalscore <- rowSums(tmp[,use])
        # }
        
        ## if compound inchikey (id) is in taxonomy.inchi vector, add 1 to total score.
        if(!is.null(taxonomy.inchi)) {
          taxonomy.score <- as.numeric(tmp[,"id"] %in% taxonomy.inchi)
          if(citation.score & any(taxonomy.score > 0)) {
            if(any(ls()=="tax.df")) {
              for(z in 1:nrow(tmp)) {
                if(taxonomy.score[z]!=1) next
                mtch <- grep(tmp[z,"id"], taxonomy.inchi)[1]
                if(length(mtch) == 0) {
                  # cat('BAD: i', i, "; z", z, tmp[z,"id"], "  ", paste(taxonomy.inchi, collapse = " "))
                  next
                }
                # cat('GOOD: i', i, "; z", z, cat(mtch, collapse = " "))
                taxonomy.score[z] <-tax.df[mtch,"weights"]
                
              }
            }
          }
          tmp$taxonomy.score <- taxonomy.score
          tmp$totalscore <- tmp$totalscore + (tmp$totalscore * 0.2 *  taxonomy.score)
        }
        
        # if(!database.score) {tmp[,"totalscore"] <- tmp[,"totalscore"] - (0.5*tmp[,"databasescore"])}
        # tmp[,"combined.score"] <- tmp[,"f.totalscore"] * tmp[,"totalscore"]
        # if(form.structure.scoring) {
        #   tmp <- tmp[order(tmp$combined.score, decreasing = TRUE),]
        # } else {
        #   tmp <- tmp[order(tmp$totalscore, decreasing = TRUE),]
        # }
        #
        # }
        
        summary[[i]] <- tmp
        tar.inchikey <- tmp$id[1]
        ramclustObj$msfinder.formula[i] <- tmp$f.name[1]
        ramclustObj$msfinder.formula.score[i] <- as.numeric(tmp$f.totalscore[1])
        best <- which(ramclustObj$msfinder.structure.details[[i]][[tmp$f.name[1]]][["structures"]][,"inchikey"] == tmp$inchikey[1])
        ramclustObj$msfinder.structure[[i]] <- ramclustObj$msfinder.structure.details[[i]][[tmp$f.name[1]]][["structures"]][best,]
        ramclustObj$msfinder.structure.fragments[[i]] <- {
          ramclustObj$msfinder.structure.details[[i]][[tmp$f.name[1]]][["fragments"]][[tmp$id[1]]]
        }
        rm(tmp)
      }
    } 
    names(summary) <- ramclustObj$cmpd
    ramclustObj$msfinder.summary <- summary
    # ramclustObj$msfinder.structure<- as.list(rep(NA, length(ramclustObj$msfinder.structure.details)))
    # names(ramclustObj$msfinder.structure) <- ramclustObj$cmpd
    
    ## look up smiles and inchi from MSFinder reference data using inchikey
    
    ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
    inchikey <- grep("inchikey", names(reference.data),  ignore.case = TRUE)
    
    if(length(inchikey) == 2) {
      ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
      inchikey.short<- inchikey[grep("short", names(reference.data)[inchikey], ignore.case = TRUE)]
      inchikey <- inchikey[which(inchikey != inchikey.short)]
    }
    
    if(length(inchikey) > 2) {
      stop("too many inchikey columns in MSFinder table - please report error to ", utils::maintainer('RAMClustR'),  '\n')
    }
    
    for(i in 1:length(ramclustObj$ann)) {
      if(i > length(ramclustObj$msfinder.structure)) next
      if( is.data.frame(ramclustObj$msfinder.structure[[i]]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
        
        tmpinch<-ramclustObj$msfinder.structure[[i]][1, "inchikey"]
        tmpinch.short<-unlist(strsplit(tmpinch, "-"))[1]
        if(use.short.inchikey & exists("inchikey.short")) {
          ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
          drow<-grep(tmpinch.short, reference.data[,inchikey.short])
        } else {
          ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
          drow<-grep(tmpinch, reference.data[,inchikey])
        }
        # d[drow,]
        # tmp<- ramclustObj$msfinder.structure[[i]][1, "inchikey"]
        
        ramclustObj$inchikey[i]<-ramclustObj$msfinder.structure[[i]][1,"inchikey"]
        ramclustObj$ann[i]<-ramclustObj$msfinder.structure[[i]][1,"name"]
        ramclustObj$smiles[i]<-ramclustObj$msfinder.structure[[i]][1,"smiles"]
        ramclustObj$annconf[i]<-2
        ramclustObj$dbid[i]<-ramclustObj$msfinder.structure[[i]][1,"resources"]
        
        if(length(drow) == 0) {
          ramclustObj$ann[i]<-ramclustObj$msfinder.structure[[i]][1, "name"]
        }
        
        if(length(drow) == 1) {
          ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
          ramclustObj$ann[i]<-reference.data[drow, "Title"]
        }
        
        if(length(drow) > 1) {
          ### ### ###  MOVE AWAY FROM REFERENCE.DATA  ### ### ### 
          n<-reference.data[drow, "Title"]
          nl<-nchar(n)
          ramclustObj$ann[i]<-n[which.min(nl)]
          ramclustObj$annconf[i]<-2
        }
        
        
      }
    }
  }
  
  if(find.inchikey) {
    
    do <- which( (ramclustObj$ann != ramclustObj$cmpd) &  is.na(ramclustObj$inchikey))
    if(length(do)>0) {
      fill.inchis <- get.inchikey.from.name(
        cmpd.names = ramclustObj$ann[do], 
        citation.weights = FALSE)
      if(nrow(fill.inchis) == 0) break
      fill.inchis <- fill.inchis[!is.na(fill.inchis$inchikey),]
      for(i in 1:length(do)) {
        mtch <- which(ramclustObj$ann[do[i]] == fill.inchis[,"cmpd.name"])
        # if(length(mtch) > 0) {break}
        if(length(mtch)==0) {next}
        ramclustObj$inchikey[do[i]] <- fill.inchis[mtch[1], "inchikey"]
      }
    }
  }
  
  
  if(formula) {
    
    if(is.null(ramclustObj$msfinder.formula)) {
      ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$msfinder.structure.details))
    }
    if(is.null(ramclustObj$msfinder.formula.score)) {
      ramclustObj$msfinder.formula.score <- rep(NA, length(ramclustObj$msfinder.structure.details))
    }
    
    dbs <- sapply(1:length(ramclustObj$msfinder.formula.details), FUN = function(x) {
      if(nrow(ramclustObj$msfinder.formula.details[[x]]) > 0) {
        paste(ramclustObj$msfinder.formula.details[[x]][, "resourcenames"], 
              collapse = ",")
      }
      else {
        NA
      }
    }
    )
    dbs <- paste(dbs, collapse = ",")
    dbs <- unlist(strsplit(dbs, ","))
    dbs <- unique(dbs[which(nchar(dbs) > 0)])
    
    if(any(database.priority == "Database ID")) {database.priority <- database.priority[-which(database.priority == "Database ID")]}
    if(any(database.priority == "custom")) {database.priority <- database.priority[-which(database.priority == "custom")]}
    if(length(database.priority) == 0) {database.priority <- NULL}
    
    suppressWarnings(if(!is.null(database.priority)) {
      if (database.priority == "all") {
        database.priority <- dbs
      }
    })
    
    ## use case insensitive matching to try to reconcile mismatches
    dbmatch <- database.priority %in% dbs
    fix <- which(!dbmatch)
    if(length(fix) > 0) {
      for(i in fix) {
        case.insensitive.match <- grep(database.priority[i], dbs, ignore.case = TRUE)
        if(length(case.insensitive.match) == 1) {
          dbs[case.insensitive.match] <- database.priority[i]
        }
      }
    }
    
    ## should case sensitive fail, remove mismatches
    dbmatch <- database.priority %in% dbs
    fix <- which(!dbmatch)
    if(length(fix) > 0) {
      warning("The following databases do not match and will be removed from consideration: ", '\n',
              "  ", database.priority[fix])
      database.priority <- database.priority[-fix]
      
    }
    
    if(is.null(ramclustObj$msfinder.formula)) {
      ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$msfinder.formula.details))
    }
    
    for(x in 1:length(ramclustObj$msfinder.formula)) { 
      if(is.na(ramclustObj$msfinder.formula[x])) {
        if (is.null(database.priority)) {
          ramclustObj$msfinder.formula[x] <- ramclustObj$msfinder.formula.details[[x]][1, "name"]
        } else {
          f <- NA
          while (is.na(f)) {
            for (i in 1:nrow(ramclustObj$msfinder.formula.details[[x]])) {
              if(nrow(ramclustObj$msfinder.formula.details[[x]])==0) next 
              for (j in 1:length(database.priority)) {
                if (grepl(database.priority[j], ramclustObj$msfinder.formula.details[[x]][i, 
                                                                                          "resourcenames"])) {
                  f <- i
                  if (!is.na(f)) {
                    break
                  }
                }
                if (!is.na(f)) {
                  break
                }
              }
              if (!is.na(f)) {
                break
              }
            }
            if (i == nrow(ramclustObj$msfinder.formula.details[[x]])) {
              f <- 1
            }
          }
          if (!is.na(f)) {
            ramclustObj$msfinder.formula[x] <- ramclustObj$msfinder.formula.details[[x]][f, "name"]
          }
          else {
            ramclustObj$msfinder.formula[x] <- NA
          }
        }
        
      }
    }
    
    for(x in 1:length(ramclustObj$msfinder.formula)) {
      df <- ramclustObj$msfinder.formula.details[[x]]
      f <- which(ramclustObj$msfinder.formula.details[[x]][, "name"] == ramclustObj$msfinder.formula[x])
      if(length(f) > 0) {
        ramclustObj$msfinder.formula.score[x] <- as.numeric(df[f, "totalscore"])
      } else  {
        ramclustObj$msfinder.formula.score[x] <- NA
      }
    }
    
  }
  
  for(i in 1:length(ramclustObj$ann)) {
    # cat("length: ", length(ramclustObj$msfinder.formula), '\n')
    if(is.null(ramclustObj$msfinder.formula)) stop("ramclustObj$msfinder.formula is null", '\n')
    if(any(!is.na(ramclustObj$msfinder.formula.details[[i]])) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
      ramclustObj$ann[i]<-ramclustObj$msfinder.formula[i]
      ramclustObj$annconf[i] <- 3
      ramclustObj$dbid[i] <- ramclustObj$msfinder.formula.details[[i]][1,"resourcenames"]
    }
  }
  
  
  if(findmain) {
    for(i in 1:length(ramclustObj$ann)) {
      if( !is.na(ramclustObj$M[i]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
        ramclustObj$ann[i]<-paste("M =", ramclustObj$M[i])
      }
    }
  }
  
  if(length(which(ramclustObj$inchikey == "undefined")) > 0) {
    ramclustObj$inchikey[which(ramclustObj$inchikey == "undefined")]<-NA
  }
  
  
  if(standardize.names) {
    cat("using pubchem PUGrest to retrieve compound names from inchikeys", '\n')
    for(i in 1:length(ramclustObj$inchikey)) {
      Sys.sleep(0.2)
      if(is.na(ramclustObj$inchikey[i])) {next}
      
      # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InChIKey/BSYNRYMUTXBXSQ-UHFFFAOYSA-N/synonyms/TXT
      syns <- tryCatch(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InChIKey/",
                                        ramclustObj$inchikey[i], 
                                        "/synonyms/TXT")),
                       error = function(x) {vector(length = 0)}, 
                       finally = function(x) {vector(length = 0)})
      if(length(syns) == 0) next
      
      if(length(syns) >=2) {
        ramclustObj$synonyms[[i]] <- syns
      }
      ramclustObj$ann[i] <- syns[1]
      
    }
    
  }
  
  ## modify compound names to make them unique
  nt <- table(ramclustObj$ann)
  ramclustObj$inchikey[which(is.na(ramclustObj$inchikey))] <- "NA"
  while(any(nt > 1)) {
    do<-which(nt>1)[1]
    mtch<-which(ramclustObj$ann == names(nt)[do])
    if(length(unique(ramclustObj$inchikey[mtch]))==1){
      ramclustObj$ann[mtch] <- paste(ramclustObj$ann[mtch], c(1:length(mtch)), sep = "_")
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
  
  ramclustObj$history$annotate <- paste( 
                               "Annotations were assigned using the RAMClustR annotate function.", 
                               " Annotation priority was assigned from higest priority to lowest:", 
                               if(any(names(ramclustObj) == "rs.lib")) {" RAMsearch, "},
                               if(any(names(ramclustObj) == "msfinder.mssearch.details")) {" MSFinder spectrum search, "},
                               if(any(names(ramclustObj) == "msfinder.structure.details")) {" MSFinder structure, "},
                               if(any(names(ramclustObj) == "msfinder.formula.details")) {" MSFinder formula, "},
                               if(any(names(ramclustObj) == "M.ann")) {" interpretMSSpectrum M."},
                               sep = ""
  )
  
  
  if(structure) {
    search.dbs <- search.dbs[-which(search.dbs == "NA")]
    ramclustObj$history$annotate2 <- paste(" MSFinder strucutures were considered from databases including", 
                                 paste(search.dbs, collapse = " "), ".", 
                                 " Database priority was set to ", 
                                 paste(priority.dbs, collapse = " "), ".",
                                 sep = "")
  }
  
  if(mssearch) {
    if(length(spec.formula.warnings) > 0) {
      warning(" The following compounds have spectral match molecular formulas", '\n',
              "   which are inconsistent with the de novo MSFinder formula results:", '\n',
              paste(ramclustObj$cmpd[spec.formula.warnings], collapse = '\n'))
    }
  }
  return(ramclustObj)
}


