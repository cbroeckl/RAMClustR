#' import.sirius
#'
#' After running Sirius on .ms files, import the annotation results
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param ms.dir optional path to .mat directory. default = "spectra/ms/out" subdirectory in working directory
#' @param ion.mode specify either "N" for negative ionization mode or "P" for positive ionization mode
#' @details this function imports the output from the Sirius program to annotate the ramclustR object
#' @return an updated ramclustR object, with new slots at $msfinder.sirius
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @importFrom utils edit read.csv read.delim read.delim2
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept annotation
#' @author Corey Broeckling
#' @export 

import.sirius <- function (
  ramclustObj = NULL,
  ms.dir = NULL,
  ion.mode = NULL
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  home.dir <-getwd()
  
  if(is.null(ms.dir)) {
    ms.dir <- paste0(getwd(),"/spectra/ms/out/")
  }

  if(!dir.exists(ms.dir)) {
    stop(cat(ms.dir, "does not exist", '\n')) 
  }

  setwd(ms.dir)
  all.dirs <- list.dirs(recursive  = FALSE)
  
  all.dirs.map <- rep(NA, length(ramclustObj$cmpd)) 
  for(i in 1:length(all.dirs.map)) {
    use <- grep(ramclustObj$cmpd[i], all.dirs)
    ## need to modify this to select newest directory.
    if(length(use) !=1) {stop("directory confusion: expected 1 match and found ", length(use), " for cmpd ", ramclustObj$cmpd[i], '\n')}
    all.dirs.map[i] <- all.dirs[use]
  }
  
  sirius.formula <- as.list(rep(NA, length(all.dirs.map)))
  names(sirius.formula) <- ramclustObj$cmpd
  sirius.structure<- sirius.formula
  
  sirius.formula.ids <- utils::read.delim("formula_identifications.tsv")
  sirius.structure.ids <- utils::read.delim("compound_identifications.tsv")
  
  if(is.null(ion.mode)) {
    ion.mode <- ramclustObj$ExpDes[[2]][which(row.names(ramclustObj$ExpDes[[2]])=="msmode"),1]
    if(startsWith(ion.mode, "N") | startsWith(ion.mode, "n")) {
      ion.mode <- "N"
    } else {
      ion.mode <- "P"
    }
  }
  
  if(ion.mode == "P") {
    suppressWarnings(sirius.fingerprint <- utils::read.delim("csi_fingerid.tsv", quote = "\\"))
    suppressWarnings(sirius.canopus <- utils::read.delim2("canopus.tsv", quote="", fill=FALSE))
  } else {
    suppressWarnings(sirius.fingerprint <- utils::read.delim("csi_fingerid_neg.tsv", quote = "\\"))
    suppressWarnings(sirius.canopus <- utils::read.delim("canopus_neg.tsv", quote="", fill=FALSE))
  }
  
  tmp <- matrix(nrow = nrow(sirius.fingerprint), ncol = length(ramclustObj$cmpd))
  dimnames(tmp)[[1]] <- rownames(sirius.fingerprint)
  dimnames(tmp)[[2]] <- ramclustObj$cmpd
  sirius.fingerprint <- cbind(sirius.fingerprint, tmp)
  
  tmp <- matrix(nrow = nrow(sirius.canopus), ncol = length(ramclustObj$cmpd))
  dimnames(tmp)[[1]] <- rownames(sirius.canopus)
  dimnames(tmp)[[2]] <- ramclustObj$cmpd
  sirius.canopus <- cbind(sirius.canopus, tmp)
  
  for(i in 1:length(all.dirs.map)) {
    
    if(i == 1) {
      config <- readLines(paste0(all.dirs.map[i], "/compound.config"))
      config <- config[!startsWith(config, "#")]
      config <- config[which(nchar(config)>0)]
    }
    
    if(!file.exists(paste0(all.dirs.map[i], "/formula_candidates.tsv"))) next
    form <- utils::read.delim(paste0(all.dirs.map[i], "/formula_candidates.tsv"))
    form <- data.frame("annotated" = rep(FALSE, nrow(form)), form)
    use <- grep(basename(all.dirs.map[i]), sirius.formula.ids$id)
    sirius.formula[[i]] <- form
    if(length(use) < 1) next
    if(length(use) > 1) stop("too many formula identifications 1:", all.dirs.map[i], sirius.formula.ids$id, '\n')
    use <- which(form$molecularFormula == sirius.formula.ids$molecularFormula[use])
    if(length(use) < 1 ) next
    if(length(use) > 1 ) stop("too many formula identifications 2:", all.dirs.map[i], sirius.formula.ids$id, '\n')
    form$annotated[use] <- TRUE
    sirius.formula[[i]] <- form
    
    fpts <- list.files(paste0(all.dirs.map[i], "/fingerprints/"))
    fpt <- fpts[grep(form[which(form$annotated), "molecularFormula"], fpts)]
    if(length(fpt)==1) {
      fpt <- as.vector(utils::read.delim(paste0(all.dirs.map[i], "/fingerprints/", fpt), header = FALSE)[,1])
      sirius.fingerprint[,ramclustObj$cmpd[i]] <- fpt
    }
    
    cnps <- list.files(paste0(all.dirs.map[i], "/canopus/"))
    cnp <- fpts[grep(form[which(form$annotated), "molecularFormula"], cnps)]
    if(length(cnp)==1) {
      cnp <- as.vector(utils::read.delim(paste0(all.dirs.map[i], "/canopus/", cnp), header = FALSE)[,1])
      sirius.canopus[,ramclustObj$cmpd[i]] <- cnp
    }
    
    if(!file.exists(paste0(all.dirs.map[i], "/structure_candidates.tsv"))) next
    struc <- suppressWarnings(utils::read.delim(paste0(all.dirs.map[i], "/structure_candidates.tsv")))
    struc <- data.frame("annotated" = rep(FALSE, nrow(struc)), struc)
    sirius.structure[[i]] <- struc
    use <- grep(basename(all.dirs.map[i]), sirius.structure.ids$id)
    if(length(use) < 1)  {
      if(nrow(struc) == 0) {next}
      stop(i)
      }
    if(length(use) > 1) stop("too many structure identifications", '\n')
    use <- which(struc$InChIkey2D == sirius.structure.ids$InChIkey2D[use])
    if(length(use) < 1 ) next
    if(length(use) > 1 ) stop("too many structure identifications", '\n')
    struc$annotated[use] <- TRUE
    sirius.structure[[i]] <- struc

  }
  
  ## organize into list and attach
  sirius <- list(length = 0)
  sirius$formula <- sirius.formula
  sirius$structure <- sirius.structure
  sirius$fingerprint <- sirius.fingerprint
  sirius$canopus <- sirius.canopus
  sirius$config <- config

  ramclustObj$sirius <- sirius
  
  ramclustObj$history$sirius <- paste(
    "Sirius (citations) was used for",
    "formula inference and computational structure assignment.",
    "Results were imported into the RAMClustR object.")
  setwd(home.dir)
  return(ramclustObj)

}

