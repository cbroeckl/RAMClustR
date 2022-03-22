#' evaluate ramSearch, MSFinder mssearch, MSFinder Structure, MSFinder Formula, and findmain output to annotate spectra of ramclustR object
#'
#' After running RAMSearch (msp) and MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param msfinder.dir full path to MSFinder directory - used for naming refinement
#' @param standardize.names logical: if TRUE, use inchikey for standardized chemical name lookup (http://cts.fiehnlab.ucdavis.edu/)
#' @param min.msms.score numerical: what is the minimum MSFinder similarity score acceptable.  default = 3.5
#' @param database.priority character.  Formula assignment prioritization based on presence in one or more databases.  Can be set to a single or multiple database names.  must match database names as they are listed in MSFinder precisely. Can also be set to 'all' (note that MSFinder reports all databases matched, not just selected databases).  If any database is set, the best formula match to that (those) database(s) is selected, rather than the best formula match overall.  
#' @param any.database.priority logical.  First priority in formula assignment is based on any of the 'database.priority' values.  Secondary priority from all other databases (determined in original MSFinder search) if TRUE.  If false, formula assignment score from MSFinder used independent of structure search results.
#' @param reset logical.  If TRUE, removes any previously assigned annotations.  
#' @details this function imports the output from the MSFinder program to annotate the ramclustR object
#' @return an updated ramclustR object, with the at $msfinder.formula, $msfinder.formula.score,  $ann, and $ann.conf slots updated to annotated based on output from 1. ramsearch output, 2. msfinder mssearch, 3. msfinder predicted structure, 4. msfinder predicted formula, and 5. interpretMSSpectrum inferred molecular weight, with listed order as priority.  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @references http://cts.fiehnlab.ucdavis.edu/static/download/CTS2-MS2015.pdf 
#' @importFrom utils edit read.csv read.delim read.delim2
#' @concept ramlclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept xcms
#' @concept MSFinder
#' @author Corey Broeckling
#' @export 


change.annotation<-function(ramclustObj = NULL,
                   msfinder.dir = "C:/MSFinder/MSFINDER ver 3.22",
                   standardize.names = FALSE,
                   min.msms.score = 3.5,
                   database.priority = "all",
                   any.database.priority = TRUE,
                   reset = TRUE
) {
  
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(!dir.exists(msfinder.dir)) {
    stop("msfinder directory does not exist: please set 'msfinder.dir' option as your full msfinder directory path")
  }
  
  if(reset) {
    ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$msfinder.structure.details))
    ramclustObj$ann <- ramclustObj$cmpd
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
    # if(length(sfile)==0) {
    #   stop("no structure DB file found in msfinder.dir / Resources")
    # }
    # 
    if(length(sfile) > 1) {
      sfile <- sfile[which.max(as.numeric(gsub(".esd", "", gsub("MsfinderStructureDB-VS", "", sfile))))]
    }
  }
  
  reference.data<-read.delim2(paste0(msfinder.dir, "/Resources/", sfile), header = TRUE, na.strings = "N/A", quote = "", stringsAsFactors = FALSE)
  
  if(any(names(ramclustObj)=="M")) {
    findmain = TRUE
  } else {findmain = FALSE}
  
  if(any(names(ramclustObj) == "msfinder.formula.details")) {
    formula = TRUE
    use.formula <- rep(FALSE, length(ramclustObj$msfinder.formula.details))
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
          form <- reference.data[which(reference.data[,"Short.InChIKey"] == tar.inchikey)[1], "Formula"]
          ramclustObj$msfinder.formula[i] <- form
        }
      }
    }
  }

  if(structure) {
    
    ### need to revisit this - not sure it is behaving appropriately.  
    tmpdb <- as.list(rep("", length(ramclustObj$msfinder.structure.details)))
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
    
    suppressWarnings(if (!is.null(database.priority)) {
      if (database.priority == "all") {
        database.priority <- dbs
      }
    })
    
    if(any(grepl('custom', database.priority))) {
      database.priority <- gsub("custom", "Database ID", database.priority)
    }
    
    dbmatch <- dbs %in% database.priority
    while (any(!dbmatch)) {
      dbmatch <- database.priority %in% dbs
      for (i in which(!dbmatch)) {
        close <- agrep(database.priority[i], dbs, max.distance = 0.2)
        fix <- readline(prompt = cat(database.priority[i], 
                                     "does not match any database names", "please type one of the following names or 'q' to quit:", 
                                     "\n", "\n", dbs, "\n"))
        if (fix == "q") {
          stop("function ended")
        }
        database.priority[i] <- fix
      }
    }
    
    search.dbs <- dbs
    priority.dbs <- database.priority
    
    use.formula <- sapply(1:length(ramclustObj$msfinder.structure.details), FUN = function(x) {
      if(any(database.priority %in% tmpdb[[x]])) {
        FALSE
      } else {TRUE}
    }
    )
    for(i in 1:length(ramclustObj$msfinder.formula)) {
      # cat(i, " ")
      if(is.na(ramclustObj$msfinder.formula[i])) {
        keep <- NA
        for(j in 1:length(ramclustObj$msfinder.structure.details[[i]])) {
          if(!is.na(keep)) break
          for(k in 1:length(database.priority)) {
            if(!is.na(keep)) break
            if(is.list(ramclustObj$msfinder.structure.details[[i]])) {
              if(any(grepl(database.priority[k], ramclustObj$msfinder.structure.details[[i]][[j]][["structures"]][,"resources"]))) {
                keep <- j
                if(!is.na(keep)) break
              }
            }
          }
        }
        if(is.na(keep)) next
        tab <- ramclustObj$msfinder.structure.details[[i]][[keep]][["structures"]]
        tar.inchikey <- tab[1,"id"]
        form <- reference.data[which(reference.data[,"Short.InChIKey"] == tar.inchikey)[1], "Formula"]
        ramclustObj$msfinder.formula[i] <- form
        # ramclustObj$msfinder.formula.score[i] <-  ## MAY NEED TO FIX THIS - CAN ALSO DO THIS ALL LATER
        rm(tab); rm(tar.inchikey); rm(form); rm(keep)
      }
    } 
    
    ramclustObj$msfinder.structure<- as.list(rep(NA, length(ramclustObj$msfinder.structure.details)))
    names(ramclustObj$msfinder.structure) <- names(ramclustObj$msfinder.structure.details)

    for(x in 1:length(ramclustObj$msfinder.structure)) {
      break.out <- FALSE
      ###################################################################################################################
      if(!is.na(ramclustObj$msfinder.formula[x])) {
        f <- ramclustObj$msfinder.formula[x]
        if(is.list(f)) stop(cat("on", x, "formula should not be a list"))
        if(!(any(names(ramclustObj$msfinder.structure.details[[x]]) == f))) {
          break.out <- TRUE
          spec.formula.warnings <- c(spec.formula.warnings, x)
        }
        if(break.out) {next}
        if(nrow(ramclustObj$msfinder.structure.details[[x]][[f]]$structures)>0) {
          dbm <- rep(FALSE, nrow(ramclustObj$msfinder.structure.details[[x]][[f]]$structures))
          for(i in 1:length(database.priority)) {
            tmpm <- grepl(database.priority[i], ramclustObj$msfinder.structure.details[[x]][[f]]$structures[,"resources"])
            dbm[which(tmpm)] <- TRUE
          }
          best <- which(dbm)[1]
          ramclustObj$msfinder.structure[[x]] <- ramclustObj$msfinder.structure.details[[x]][[f]][["structures"]][best,]
          ramclustObj$msfinder.structure.fragments[[x]] <- {
            ramclustObj$msfinder.structure.details[[x]][[f]][["fragments"]][[ramclustObj$msfinder.structure[[x]][best,"id"]]]
          }
        }
      }
    }
    
    if(any.database.priority) {
      for(x in 1:length(ramclustObj$msfinder.structure)) {
        # cat(x, " ")
        
        if(is.na(ramclustObj$msfinder.formula[x])) {
          
          if(is.list(ramclustObj$msfinder.structure.details[[x]])) {
            nms <- names(ramclustObj$msfinder.structure.details[[x]]) 
            do.form <- nms[-which(nms == "Spectral DB search")]
            keep <- sapply(do.form, FUN = function(y) {
              (nrow(ramclustObj$msfinder.structure.details[[x]][[y]][["structures"]]) > 0)
            }
            )
            if(any(keep)) {
              # stop('found one!')
              best <- which(keep)[1]
              ramclustObj$msfinder.formula[x] <- names(best)
              ramclustObj$msfinder.structure[[x]] <- ramclustObj$msfinder.structure.details[[x]][[names(best)]]$structures[1,]
              ramclustObj$msfinder.structure.fragments[[x]] <- ramclustObj$msfinder.structure.details[[x]][[names(best)]]$fragments[[1]]
            }
          }
        }
      }
    }
    
    inchikey <- grep("inchikey", names(reference.data),  ignore.case = TRUE)
    
    if(length(inchikey) == 2) {
      inchikey.short<- inchikey[grep("short", names(reference.data)[inchikey], ignore.case = TRUE)]
      inchikey <- inchikey[which(inchikey != inchikey.short)]
    }
    
    if(length(inchikey) > 2) {
      stop("too many inchikey columns in MSFinder table - please report error to ", utils::maintainer('RAMClustR'),  '\n')
    }
    
    for(i in 1:length(ramclustObj$ann)) {
      if( is.data.frame(ramclustObj$msfinder.structure[[i]]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
        
        tmpinch<-ramclustObj$msfinder.structure[[i]][1, "inchikey"]
        tmpinch.short<-unlist(strsplit(tmpinch, "-"))[1]
        if(use.short.inchikey & exists("inchikey.short")) {
          drow<-grep(tmpinch.short, reference.data[,inchikey.short])
        } else {
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
          ramclustObj$ann[i]<-reference.data[drow, "Title"]
        }
        
        if(length(drow) > 1) {
          n<-reference.data[drow, "Title"]
          nl<-nchar(n)
          ramclustObj$ann[i]<-n[which.min(nl)]
          ramclustObj$annconf[i]<-2
        }
        
        
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
    
    suppressWarnings(if (!is.null(database.priority)) {
      if (database.priority == "all") {
        database.priority <- dbs
      }
    })
    
    dbmatch <- dbs %in% database.priority
    
    while (any(!dbmatch)) {
      dbmatch <- database.priority %in% dbs
      for (i in which(!dbmatch)) {
        close <- agrep(database.priority[i], dbs, max.distance = 0.2)
        fix <- readline(prompt = cat(database.priority[i], 
                                     "does not match any database names", "please type one of the following names or 'q' to quit:", 
                                     "\n", "\n", dbs, "\n"))
        if (fix == "q") {
          stop("function ended")
        }
        database.priority[i] <- fix
      }
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
    if( !is.na(ramclustObj$msfinder.formula[[i]]) && (ramclustObj$cmpd[i] == ramclustObj$ann[i]) )  {
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
    cat("using chemical translation service - requires interet access and may take a few minutes to complete", '\n')
    
    inchikey2inchi<-which(!is.na(ramclustObj$inchikey) & is.na(ramclustObj$inchi))
    for(i in inchikey2inchi) {
      if(!is.na(ramclustObj$inchikey[i])) {
        
        link <- paste0("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/InChI Code/", 
                       # unlist(strsplit(ramclustObj$inchikey[i], "-"))[1])
                       ramclustObj$inchikey[i])
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
        if(length(inchis) == 0 ) {
          ramclustObj$inchi[i] <- NA
        }
        if(length(inchis)>=1) {
          ramclustObj$inchi[i] <- inchis[1]
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
  
  ramclustObj$history <- paste(ramclustObj$history, 
                               " Annotations were assigned using the RAMClustR annotate function.", 
                               " Annotation priority was assigned from higest priority to lowest:", 
                               if(any(names(ramclustObj) == "rs.lib")) {" RAMsearch, "},
                               if(any(names(ramclustObj) == "msfinder.mssearch.details")) {" MSFinder spectrum search, "},
                               if(any(names(ramclustObj) == "msfinder.structure.details")) {" MSFinder structure, "},
                               if(any(names(ramclustObj) == "msfinder.formula.details")) {" MSFinder formula, "},
                               if(any(names(ramclustObj) == "M.ann")) {" interpretMSSpectrum M."},
                               sep = ""
  )
  search.dbs <- search.dbs[-which(search.dbs == "NA")]
  ramclustObj$history <- paste(ramclustObj$history, " MSFinder strucutures were considered from databases including", 
                               paste(search.dbs, collapse = " "), ".", 
                               " Database priority was set to ", 
                               paste(priority.dbs, collapse = " "), ".",
                               " Any.database.priority was set to ", any.database.priority, ".",
                               sep = "")
  if(length(spec.formula.warnings) > 0) {
    warning(" The following compounds have spectral match molecular formulas", '\n',
            "   which are inconsistent with the de novo MSFinder formula results:", '\n',
            paste(ramclustObj$cmpd[spec.formula.warnings], collapse = '\n'))
  }
  return(ramclustObj)
}

