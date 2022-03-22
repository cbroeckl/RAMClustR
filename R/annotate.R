#' evaluate ramSearch, MSFinder mssearch, MSFinder Structure, MSFinder Formula, and findmain output to annotate spectra of ramclustR object
#'
#' After running RAMSearch (msp) and MSFinder on .mat or .msp files, import the spectral search results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param standardize.names logical: if TRUE, use inchikey for standardized chemical name lookup (http://cts.fiehnlab.ucdavis.edu/)
#' @param min.msms.score numerical: what is the minimum MSFinder similarity score acceptable.  default = 6.5
#' @param database.priority character.  Formula assignment prioritization based on presence in one or more (structure) databases.  Can be set to a single or multiple database names.  must match database names as they are listed in MSFinder precisely. Can also be set to 'all' (note that MSFinder reports all databases matched, not just databases in MSFinder parameters).  If any database is set, the best formula match to any of those databases is selected, rather than the best formula match overall.  If NULL, this will be set to include all selected databases (from ramclustObj$msfinder.dbs, retrieved from search output during import.msfinder.formulas(), when available) or 'all'.  
#' @param database.priority.factor numeric, between 0 and 1.  0.1 by default.  The proportion by which scores for structures not in priority database are assessed
#' @param taxonomy.inchi vector or data frame.  Only when rescore.structure = TRUE.  user can supply a vector of inchikeys.  If used, structures which match first block of inchikey retain full score, while all other structures are penalized.  
#' @param taxonomy.inchi.factor numeric, between 0 and 1.  0.1 by default.  The proportion by which scores for structures not in taxonomy.inchi vector are assessed
#' @param find.inchikey logical.  default = TRUE. use chemical translation service to try to look up inchikey for chemical name.
#' @param use.ri logical.  default = TRUE.  If retention index is available in ramclustObj (set by 'rc.calibrate.ri') and in library spectra from MSFinder, use RI similiarity to rescore.
#' @param sri numeric.  sigma value for retention index. controls decay rate of retention index curve. decay rate between 0 and 1 exported, and multiplied by spectrum score, totalscore.
#' @param ri.na.factor numeric. between 0 and 1.  0.5 by default.  how should spectrum scores be treated when no retention index is available?  NA values are replaced by retention index similarities of ri.na.factor when use.ri = TRUE.
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
                   standardize.names = FALSE,
                   min.msms.score = 0.8,
                   database.priority = NULL,
                   database.priority.factor = 0.1,
                   find.inchikey = TRUE,
                   taxonomy.inchi = NULL,
                   taxonomy.inchi.factor = 0.1,
                   use.ri = TRUE,
                   sri = 300,
                   ri.na.factor = 0.6,
                   reset = TRUE
) {
  
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  ## reset annotation slots to ensure that all annotation data reflects new processing
  if(reset) {
    ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$formula <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$annconf <- rep(4, length(ramclustObj$cmpd))
    ramclustObj$ann <- ramclustObj$cmpd
    ramclustObj$inchikey <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$inchi <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$smiles <- rep(NA, length(ramclustObj$cmpd))
    ramclustObj$pubchem.cid <- rep(NA, length(ramclustObj$cmpd))
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
  
  
  suppressWarnings(
    have.internet <- !as.logical(system(paste("ping -n 1", "www.google.com"), show.output.on.console = FALSE))
  )
  
  
  ## make sure taxonomy inchikeys are properly formatted
  if(!is.null(taxonomy.inchi)) {
    if(is.data.frame(taxonomy.inchi)) {
      tax.df <- taxonomy.inchi
      taxonomy.inchi <- tax.df[,"inchikey"]
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
  
  ## determine what annotation data we have available
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
  
  if(any(names(ramclustObj) == "sirius")) {
    sirius = TRUE
  } else {sirius = FALSE}
  
  ## add new (or replace existing) slots for holding annotation data from external tools
  if(!any(names(ramclustObj) == "inchikey")) {ramclustObj$inchikey <- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "inchi")) {ramclustObj$inchi<- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "smiles"))  {ramclustObj$smiles <- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "dbid")) {ramclustObj$dbid<- rep(NA, length(ramclustObj$cmpd))}
  if(!any(names(ramclustObj) == "synonyms")) {
    ramclustObj$synonyms <- as.list(rep(NA, length(ramclustObj$cmpd)))
  }
  
  ## establish databases (MSFinder) that are considered of higher priority. 
  if(!is.null(ramclustObj$msfinder.dbs) & is.null(database.priority)) {
    database.priority <- ramclustObj$msfinder.dbs
  }
  
  ## if ramsearch has been imported, annotate those compounds first
  ## this is prioritized due to the manual nature of the process.  
  if(ramsearch) {
    out <- ramclustObj$rs.out
    
    
    ##these items will be filled and added to the RC object
    ramclustObj$rs.spec	<-as.list(rep("", max(ramclustObj$featclus)))
    ramclustObj$rs.lib	<-rep("", max(ramclustObj$featclus))
    ramclustObj$rs.specn	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.libn	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.mf	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.rmf	<-as.integer(rep(-1, max(ramclustObj$featclus)))
    ramclustObj$rs.prob	<-as.numeric(rep(-1, max(ramclustObj$featclus)))
    
    ##pull relevant line numbers for all spectra
    ##name line is first, so call that directly, range between name 1 and name 2 is the
    ##range of spectrum 1
    name<-which(regexpr('Matched Spectrum:', out)==1)
    ann<-which(regexpr('Annotation:', out)==1)
    origname<-which(regexpr('Original Name', out)==1)
    
    if(any((origname-name)!=2)) stop("please don't edit the output from ramsearch manually")
    if(any(ramclustObj$cmpd!=sub("Original Name: ", "", out[origname]))) {stop("compound names/order differ between ramclust object and ramsearch output")}  
    
    ##if length of name is not equal to length of ramclust Object 'cmpd' slot, something is wrong:
    if(length(name)/as.integer(as.character(ramclustObj$ExpDes[[2]]["MSlevs",1]))!= length(ramclustObj$cmpd)) stop("number of spectra in ramsearch output different than number of compounds in ramclust object")
    
    ##now pull relevant info our for each spectrum in output
    for(i in 1:length(ann)) {
      md<-out[name[i]:min((name[i+1]-1), length(out), na.rm=TRUE)]
      cname<-sub("Original Name: ", "", md[grep("Original Name: ", md)])
      ind<-as.numeric(sub("C", "", cname))
      mf<-ramclustObj$rs.mf[ind]
      newmf<-as.integer(sub("Match Factor / Dot Product: ", "", md[grep("Match Factor / Dot Product: ", md)]))
      if(ramclustObj$cmpd[ind] != cname) {
        stop(paste("something is amiss with compound ", i, ": the names do not match", sep=""))
      }
      
      ## note that inchikey returned from GOLM is the metabolite inchikey, not derivative inchikey
      ## inchikey for NIST? This incongruity may cause issues later
      ## i.e. modeling retention time/index from compound properties would benefit from derivative inchikey
      ## while metabolic networks benefit from metabolite inchikey.  both would be valuable.  
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
  
  ## Annotate based on msfinder spectral matches next
  if(mssearch) {
    spec.formula.warnings <- vector(length = 0, mode = 'numeric')
    if(is.null(ramclustObj$msfinder.formula)) {
      ramclustObj$msfinder.formula <- rep(NA, length(ramclustObj$msfinder.structure.details))
    }
    
    
    ## get best score from all spectral matches
    msfinder.mssearch.score<-as.numeric(sapply(1:length(ramclustObj$msfinder.mssearch.details), FUN = function(x) {
      if(is.null(nrow(ramclustObj$msfinder.mssearch.details[[x]]$summary))) {
        NA 
      } else {
        ramclustObj$msfinder.mssearch.details[[x]]$summary[1,"totalscore"]
      }
    }
    )
    )
    
    ## low res EI-GC-MS scoring maxes out at 0.67 for total score, if 
    ## scores suggest this to be the case, divide all scores by 0.67 to bring
    ## max score t0 1.
    
    if(max(msfinder.mssearch.score, na.rm = TRUE) < 0.68) {
      for(x in 1:length(ramclustObj$msfinder.mssearch.details)) {
        if(nrow(ramclustObj$msfinder.mssearch.details[[x]]$summary)==0) next
        ramclustObj$msfinder.mssearch.details[[x]]$summary[,"totalscore"] <- 
          ramclustObj$msfinder.mssearch.details[[x]]$summary[,"totalscore"]/0.67
      }
      msfinder.mssearch.score<-as.numeric(
        sapply(1:length(ramclustObj$msfinder.mssearch.details), 
               FUN = function(x) {
                 if(is.null(nrow(ramclustObj$msfinder.mssearch.details[[x]]$summary))) {
                   NA 
                 } else {
                   ramclustObj$msfinder.mssearch.details[[x]]$summary[1,"totalscore"]
                 }
               }
        )
      )
    }
    
    if(use.ri) {
      if(is.null(ramclustObj$clri)) warning("no retention index available, not using RI weighting")
    }
    if(use.ri & !is.null(ramclustObj$clri)) {
      for(x in 1:length(ramclustObj$msfinder.mssearch.details[[x]])) {
        spec.sim <- ramclustObj$msfinder.mssearch.details[[x]]$summary$totalscore
        ramclustObj$msfinder.mssearch.details[[x]]$summary$spectrum.score <- spec.sim
        ri <- ramclustObj$msfinder.mssearch.details[[x]]$summary$retentionindex
        ri[which(ri == -1)] <- NA
        ri.sim <- round(exp(-(( 
          abs(
            ri - ramclustObj$clri[x]
          )
        )^2)/(2*(sri^2))), 
        
        digits=3 )
        ri.sim[which(is.na(ri.sim))] <- ri.na.factor
        ramclustObj$msfinder.mssearch.details[[x]]$summary$retentionindexsim <- ri.sim
        ramclustObj$msfinder.mssearch.details[[x]]$summary$totalscore <- {
          ramclustObj$msfinder.mssearch.details[[x]]$summary$totalscore * ri.sim
        }
        ramclustObj$msfinder.mssearch.details[[x]]$summary[order(
          ramclustObj$msfinder.mssearch.details[[x]]$summary$totalscore, decreasing = TRUE),]
      }
    }
    
    ramclustObj$msfinder.mssearch.score <- msfinder.mssearch.score
    
    
    for(i in 1:length(ramclustObj$ann)) {
      if((nrow(ramclustObj$msfinder.mssearch.details[[i]]$summary)>0) & (ramclustObj$cmpd[i] == ramclustObj$ann[i]))  {
        if(ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"totalscore"] >= min.msms.score ) {
          ramclustObj$inchikey[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"inchikey"]
          ramclustObj$smiles[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"smiles"]
          ramclustObj$ann[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"name"]
          ramclustObj$annconf[i]<-"2a"
          ramclustObj$dbid[i]<-ramclustObj$msfinder.mssearch.details[[i]]$summary[1,"resources"]
          if(grepl("PUBCHEM", ramclustObj$dbid[i])) {
            cid <- unlist(strsplit(ramclustObj$dbid[i], ";"))
            cid <- cid[grep("PUBCHEM", cid)[1]]
            cid <- trimws(gsub("PUBCHEM", "", cid))
            if(cid != "CID") {
              ramclustObj$pubchem.cid[i] <- cid
            }
          }
        }
      }
    }
    inchikey <- which(!is.na(ramclustObj$inchikey) & is.na(ramclustObj$msfinder.formula))
    inchikey <- unique(ramclustObj$inchikey[inchikey])
    inchikey <- inchikey[which(nchar(inchikey)==27 & stringr::str_count(inchikey, "-")==2)]
    
    if(length(inchikey) > 0) {
      if(use.short.inchikey) {
        inchikey <- sapply(1:length(inchikey), FUN = function(x) {
          unlist(strsplit(inchikey[x], "-"))[1]
        })
        inchikey <- unique(inchikey)
        inchikey <- inchikey[which(nchar(inchikey)==14)]
      }
      form <- rep(NA, length(inchikey))
      tmp.cid <- rep(NA, length(inchikey))
      for(i in 1:length(inchikey)) {
        tmp.pc <- tryCatch(
          jsonlite::read_json(
            paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",
                   inchikey[i],
                   "/property/MolecularFormula/JSON"
            )
          ),  # $PropertyTable$Properties[[1]]$MolecularFormula[1]
          error=function(cond) {
            return(NA)
          },
          warning=function(cond) {
            return(NA)
          })
        form[i] <- tmp.pc$PropertyTable$Properties[[1]]$MolecularFormula[1]
        tmp.cid[i] <- tmp.pc$PropertyTable$Properties[[1]]$CID[1]
        
        Sys.sleep(0.2)
      }
      for(i in 1:length(form)) {
        if(is.na(form[i])) next
        ramclustObj$msfinder.formula[grepl(inchikey[i], ramclustObj$inchikey)] <- form[i]
        ramclustObj$formula[grepl(inchikey[i], ramclustObj$inchikey)] <- form[i]
        ramclustObj$pubchem.cid[grepl(inchikey[i], ramclustObj$inchikey)] <- tmp.cid[i]
      }
    }
  }
  
  
  ## Sirius section
  if(sirius) {
    do.ann <- 1:length(ramclustObj$ann)
    for(i in 1:length(do.ann)) {
      if(ramclustObj$ann[i] != ramclustObj$cmpd[i]) next
      if(is.data.frame(ramclustObj$sirius$formula[[i]])) {  ## if formula empty, do not annotate
        if(any(ramclustObj$sirius$formula[[i]]$annotated)) {  ## if formula is annotated proceed
          use.form <- which(ramclustObj$sirius$formula[[i]]$annotated)
          ramclustObj$formula[i] <- ramclustObj$sirius$formula[[i]]$molecularFormula[use.form]
          ramclustObj$ann[i]     <- ramclustObj$sirius$formula[[i]]$molecularFormula[use.form]
          ramclustObj$annconf[i] <- "4"
          if(is.data.frame(ramclustObj$sirius$structure[[i]])) {  ## if structure is empty, do not annotate further.
            if(any(ramclustObj$sirius$structure[[i]]$annotated)) {   ## if structure is annotated, proceed
              st.use <- which(ramclustObj$sirius$structure[[i]]$annotated)
              
              # if(have.internet) {
              #   cid <- ramclustObj$sirius$structure[[i]][st.use, "pubchemids"]
              #   cid <- unlist(strsplit(cid, ";"))
              #   tmp.pc <- rc.cmpd.get.pubchem(cmpd.cid = cid[1], write.csv = FALSE, get.vendors = FALSE, get.bioassays = FALSE, get.synonyms = FALSE)
              #   ramclustObj$ann[i] <- tmp.pc$pubchem$pubchem.name
              #   ramclustObj$inchikey[i] <- tmp.pc$properties$InChIKey
              #   ramclustObj$annconf[i] <- "2b"
              #   ramclustObj$smiles[i] <- tmp.pc$properties$CanonicalSMILES
              #   ramclustObj$inchi[i] <- tmp.pc$properties$InChI
              #   ramclustObj$formula[i] <- tmp.pc$properties$MolecularFormula
              #   ramclustObj$pubchem.cid <- tmp.pc$pubchem$cid
              #   
              # } else {
              ramclustObj$ann[i] <- ramclustObj$sirius$structure[[i]][st.use, "name"]
              ramclustObj$inchikey[i] <- ramclustObj$sirius$structure[[i]][st.use, "InChIkey2D"]
              ramclustObj$annconf[i] <- "2b"
              ramclustObj$smiles[i] <- ramclustObj$sirius$structure[[i]][st.use, "smiles"]
              ramclustObj$inchi[i] <- ramclustObj$sirius$structure[[i]][st.use, "InChI"]
              ramclustObj$formula[i] <- ramclustObj$sirius$structure[[i]][st.use, "molecularFormula"]
              tmp.cid <- ramclustObj$sirius$structure[[i]][st.use, "pubchemids"]
              if(grepl(";", tmp.cid)) {
                tmp.cid <- unlist(strsplit(tmp.cid, ";"))[1]
              }
              ramclustObj$pubchem.cid <- tmp.cid
              
              #}
              
              
              
            }
          }
        }
      }
    }
  }
  
  ## MSFinder structure section
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
    
    ## set database priority, if not manually set, to the databases selected when running msfinder.
    ## else use the full db list
    suppressWarnings(
      if(is.null(database.priority)) {
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
    
    ## create template to hold called annotations
    template <- matrix(nrow = 0, ncol = 0)
    
    ## annotation section 
    for(i in 1:length(ramclustObj$ann)) {
      # cat(i, "  ")
      # db.m <- data.frame(lapply(1:length(database.priority), FUN = function(x) {
      #   grepl(database.priority[x], ramclustObj$msfinder.formula.details[[i]][,"resourcenames"])
      # }
      # )
      # )
      
      # get all plausible formulas
      forms <- names(ramclustObj$msfinder.structure.details[[i]])
      if(any(forms == "Spectral DB search")) {
        forms <- forms[-which(forms == "Spectral DB search")]
      }
      
      if(length(forms) == 0) next
      
      for(j in 1:length(forms)) {
        if(is.null(ramclustObj$msfinder.structure.details[[i]][forms[j]][[1]]$structures)) next
        tmp <- ramclustObj$msfinder.structure.details[[i]][forms[j]][[1]]$structures
        if(!exists('str.sum')) {
          str.sum <- cbind(data.frame("formula" = rep("", 0)), tmp[0, ])
        } 
        if(nrow(tmp) == 0) next
        form.score <- ramclustObj$msfinder.formula.details[[i]]
        form.score <- as.numeric(form.score[which(form.score$name == forms[j]),"totalscore"])
        str.sum.tmp <- cbind(data.frame(
          "formula" = rep(forms[j], nrow(tmp)),
          "formula.score" = rep(form.score, nrow(tmp)))
          , tmp)
        str.sum <- rbind(str.sum, str.sum.tmp)
      }
      
      
      num.cols <- c("retentiontime", "retentionindex", "totalbondenergy",
                    "totalscore", "totalhrrulesscore", "totalbondcleavagescore", "totalmassaccuracyscore",
                    "totalfragmentlinkagescore", "totalbonddissociationenergyscore", "databasescore",
                    "substructureassignmentscore", "rtsimilarityscore", "risimilarityscore"
      )
      for(k in num.cols) {
        str.sum[,k] <- as.numeric(str.sum[,k])
      }
      
      ## make new totalscore which is product of formula score and structure score
      str.sum$totalscore.structure <- str.sum$totalscore
      str.sum$totalscore <- str.sum$totalscore.structure * str.sum$formula.score
      
      
      ## adjust scores if priority databases are specified
      if(!is.null(priority.dbs) & (nrow(str.sum)>0)) {
        priority <- lapply(1:length(priority.dbs), FUN = function(k) {grepl(priority.dbs[k], str.sum$resources)})
        priority <- data.frame(priority)
        priority <- apply(priority, 1, "any")
        priority.score <- rep(1-database.priority.factor, length(priority))
        priority.score[priority] <- 1
        new.totalscore <- str.sum[,"totalscore"] * priority.score
        totalscore <- str.sum[,"totalscore"]
        str.sum <- data.frame(str.sum, 
                              "db.priority.score" = priority.score,
                              "msfinder.totalscore" = totalscore)
        str.sum$totalscore <- new.totalscore
        
      }
      
      ## adjust scores if taxonomy inchikeys are provided
      if(!is.null(taxonomy.inchi)& (nrow(str.sum)>0)) {
        priority <- sapply(1:nrow(str.sum), FUN = function(k) {any(grepl(str.sum$inchikey[k], taxonomy.inchi))})
        priority.score <- rep(1-database.priority.factor, length(priority))
        priority.score[priority] <- 1
        new.totalscore <- str.sum[,"totalscore"] * priority.score
        if(any(names(str.sum) == "msfinder.totalscore")) {
          totalscore <- str.sum[,"msfinder.totalscore"]
        } else {
          totalscore <- str.sum[,"totalscore"]
        }
        
        str.sum <- data.frame(str.sum, 
                              "tax.priority.score" = priority.score)
        if(!any(names(str.sum) == "msfinder.totalscore")) {
          str.sum <- data.frame(str.sum,"msfinder.totalscore" = totalscore)
        }
        str.sum$totalscore <- new.totalscore
      }
      
      str.sum <- str.sum[order(str.sum$totalscore, decreasing = TRUE),]
      # head(str.sum, n = 4)
      
      
      summary[[i]] <- str.sum
      
      
      if(nrow(str.sum) == 0) next
      
      ramclustObj$ann[i] <- str.sum[1,"name"]
      ramclustObj$annconf[i] <- 2
      ramclustObj$msfinder.formula[i] <- str.sum[1,"formula"]
      ramclustObj$msfinder.formula.score[i] <- str.sum[1,"formula.score"]
      ramclustObj$inchikey[i] <- str.sum[1,"inchikey"]
      ramclustObj$smiles[i] <- str.sum[1,"smiles"]
      ramclustObj$dbid[i] <- str.sum[1,"resources"]
      
      rm(str.sum); rm(tmp)
    }
    
    names(summary) <- ramclustObj$cmpd
    ramclustObj$msfinder.summary <- summary
  }
  
  if(find.inchikey) {
    
    do <- which( (ramclustObj$ann != ramclustObj$cmpd) &  is.na(ramclustObj$inchikey))
    if(length(do)>0) {
      fill.inchis <- rc.cmpd.get.pubchem(cmpd.names = ramclustObj$ann[do],  ## just set cid.l to 1
                                         use.parent.cid = FALSE,
                                         manual.entry = FALSE,
                                         get.vendors = FALSE,
                                         get.properties = TRUE,
                                         all.props = FALSE,
                                         get.bioassays = FALSE,
                                         get.synonyms = FALSE)
      ramclustObj$inchikey[do] <- fill.inchis$properties$InChIKey
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
  
  ramclustObj$inchikey[which(ramclustObj$inchikey == "NA")] <- NA
  
  ## replace 2d inchikey with 3d based on pubchem
  if(have.internet) {
    do.inch <- which(nchar(ramclustObj$inchikey) > 0 & nchar(ramclustObj$inchikey) < 27)
    for(i in 1:length(do.inch)) {
      
      tmp.inchikey <- tryCatch(
        readLines(
          paste0(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",
            ramclustObj$inchikey[do.inch[i]],
            "/property/inchikey/txt"
          )
        ),  # $PropertyTable$Properties[[1]]$MolecularFormula[1]
        error=function(cond) {
          return(NA)
        },
        warning=function(cond) {
          return(NA)
        })
      suppressWarnings(if(!is.na(tmp.inchikey)) {
        ramclustObj$inchikey[do.inch[i]] <- tmp.inchikey[1]
      })
      Sys.sleep(0.2)
    }
  }
  
  
  if(standardize.names & have.internet) {
    cat("using pubchem PUGrest to retrieve compound names from inchikeys", '\n')
    do.inchi <- which(!is.na(ramclustObj$inchikey))
    tmp <- rc.cmpd.get.pubchem(cmpd.inchikey = ramclustObj$inchikey[do.inchi],
                               use.parent.cid = FALSE,
                               manual.entry = FALSE,
                               get.vendors = FALSE,
                               get.properties = FALSE,
                               all.props = FALSE,
                               get.bioassays = FALSE,
                               get.synonyms = FALSE)
    keep <- which(!is.na(tmp$pubchem$pubchem.name))
    ramclustObj$ann[do.inchi[keep]] <- tmp$pubchem$pubchem.name[keep]
  }
  
  ## modify compound names to make them unique
  ramclustObj$ann <- make.unique(ramclustObj$ann)
  
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


