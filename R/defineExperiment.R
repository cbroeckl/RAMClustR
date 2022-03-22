#' defineExperiment
#'
#' Create an Experimental Design R object for record-keeping and msp output
#'
#'
#' @param csv logical or filepath.  If   csv = TRUE , csv template called "ExpDes.csv" will be written to your working directory.  you will fill this in manually, ensuring that when you save you retain csv format.  ramclustR will then read this file in and and format appropriately.  If csv = FALSE, a pop up window will appear (in windows, at least) asking for input.  If a character string with full path (and file name) to a csv file is given, this will allow you to read in a previously edited csv file. 
#' @param force.skip logical.  If TRUE, ramclustR creates a pseudo-filled ExpDes object to enable testing of functionality. Not recommended for real data, as your exported spectra will be improperly labelled.
#' @return an Exp Des R object which will be used for record keeping and writing spectra data.  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
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

defineExperiment<-function(csv = FALSE, force.skip=FALSE) {
  LCMS <- data.frame("value" = c(chrominst="",
                                 msinst="",
                                 column="",
                                 solvA="",
                                 solvB="",
                                 CE1="",
                                 CE2="",
                                 mstype="",
                                 msmode="",
                                 ionization="",
                                 colgas="",
                                 msscanrange="",
                                 conevolt="",
                                 MSlevs=1), stringsAsFactors = FALSE)
  
  GCMS <- data.frame("value" = c(chrominst="",
                                 msinst="",
                                 column="",
                                 InletTemp="",
                                 TransferTemp="",
                                 mstype="",
                                 msmode="",
                                 ionization="",
                                 msscanrange="",
                                 scantime="",
                                 deriv="",
                                 MSlevs=1), stringsAsFactors = FALSE)
  
  Experiment<-data.frame("Value" = rep("", 5),
                         "Description" = c("experiment name, no spaces",
                                           "species name",
                                           "sample type",
                                           "individual and/or organizational affiliation",
                                           "GC-MS or LC-MS"), 
                         row.names = c("Experiment",
                                       "Species",
                                       "Sample",
                                       "Contributor",
                                       "platform"))
  
  
  
  
  if(force.skip == TRUE) {
    csv <- TRUE
  }
  
  if (is.logical(csv)) {
    if(csv) {
      # cat(find.package("RAMClustR"))
      #out<-read.csv(paste(find.package("RAMClustR"), "/params/params.csv", sep=""), header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)
      out <- structure(list(
        parameter = c("", "Experimental Design: ", "Experiment", 
                      "Species", "Sample", "Contributer", "platform", "", "GC-MS ", 
                      "chrominst", "msinst", "column", "InletTemp", "TransferTemp", 
                      "mstype", "msmode", "ionization", "msscanrange", "scantime", 
                      "deriv", "MSlevs", "", "LC-MS", "chrominst", "msinst", "column", 
                      "solvA", "solvB", "CE1", "CE2", "mstype", "msmode", "ionization", 
                      "colgas", "msscanrange", "conevol", "MSlevs"), 
        Value = c("", 
                  "", "fill", "fill", "fill", "fill", "fill", "", "", "fill", "fill", 
                  "fill", "fill", "fill", "fill", "fill", "fill", "fill", "fill", 
                  "fill", "1", "", "", "fill", "fill", "fill", "fill", "fill", 
                  "fill", "fill", "fill", "fill", "fill", "fill", "fill", "fill", 
                  "1"), 
        Description = c("", "", "experiment name, no spaces", "Genus species from which samples are derived", 
                        "Type of sample (i.e., serum, leaf)", "you or your PI's name", 
                        "either GC-MS or LC-MS", "", "", "model of LC instrument", "model of MS instrument", 
                        "column description", "temperature of inlet", "temperature of GC to MS transfer line", 
                        "Type of mass spectrometer (one of QQQ, TOF, QTOF, Orbi, Q)", 
                        "positive or negative ion mode", "EI, AP, or CI, typically", 
                        "scan range used for acquisition", "time for each full scan spectrum (i.e 0.2 seconds)", 
                        "derivitization used ( \"TMS\", \"None\", \"TBDMS\".)", "number of levels of energy acquired - 1 typically", 
                        "", "", "model of LC instrument", "model of MS instrument", "column description", 
                        "Solvent A composition", "Solvent B composition", "Collision energy of acquisition of MS data", 
                        "Collision energy of acquisition for MSe/idMSMS data (when applicable)", 
                        "Type of mass spectrometer (one of QQQ, TOF, QTOF, Orbi, Q)", 
                        "positive or negative ion mode", "ESI or APCI, typically", "gas used for collisional dissociation", 
                        "scan range used for acquisition", "cone voltage used for acquisition", 
                        "number of levels of energy in XCMS object data - 1 typically"
        )), class = "data.frame", row.names = c(NA, -37L))
      
      if(!force.skip) {
        write.csv(out, file=paste(getwd(), "/ExpDes.csv", sep=""), row.names=FALSE)
        readline(prompt=cat("A file called ExpDes.csv has been written to your working directorty:",
                            '\n', '\n',
                            getwd(), 
                            '\n', '\n',
                            "please replace platform appropriate 'fill' cells with instrument and experiment",
                            '\n', "data and save file.  When complete, press [enter] to continue"
        ))
        csv.in<-read.csv(file=paste(getwd(), "/ExpDes.csv", sep=""), header=TRUE, check.names=FALSE, stringsAsFactors = FALSE) 
      } else {
        csv.in <-out
        csv.in[7,2]<-'LC-MS'
        csv.in[which(csv.in[,1]=="MSlevs"),2]<-1
      }
      design<-data.frame("value" = csv.in[3:7,2], row.names = csv.in[3:7,1], stringsAsFactors = FALSE)
      
      instrument <- NULL
      plat<-as.character(design[5,1])
      if(grepl("LC-MS", plat)) {
        instrument<-"LC-MS"
      }  
      if( grepl("GC-MS", plat )) {
        instrument<-"GC-MS"
      }  
      if(!grepl("LC-MS", plat ) & !grepl("GC-MS", plat ) ) {
        if(grepl('[Gg]',  substring(plat,1,1))) {
          instrument<-"GC-MS"
        }
        if(grepl('[Ll]',  substring(plat,1,1))) {
          instrument<-"LC-MS"
        }
      }
      if(is.null(instrument)) {
        stop("do not regonize instrument platform, please use either 'GC-MS' or 'LC-MS' " )
      }
      
      rowstart<-grep(instrument, csv.in[,1])+1
      rowend<-grep("MSlevs", csv.in[,1])
      rowend<-rowend[which(rowend > rowstart)]
      if(length(rowend)>1) {
        rowend<-rowend[which.min((rowend - rowstart))]
      }
      instrument <- data.frame('value' = csv.in[rowstart:rowend,2], row.names = csv.in[rowstart:rowend,1], stringsAsFactors = FALSE)
      ExpDes <- list("design" = design, "instrument" = instrument)
      
    }  
    
    if(!csv){
      
      suppressWarnings( design<-edit(Experiment))
      
      plat<-as.character(design[5,1])
      if(grepl("LC-MS", plat)) {
        instrument<-"LC-MS"
      }  
      if( grepl("GC-MS", plat )) {
        instrument<-"GC-MS"
      }  
      if(!grepl("LC-MS", plat ) & !grepl("GC-MS", plat ) ) {
        if(grepl('[Gg]',  substring(plat,1,1))) {
          instrument<-"GC-MS"
        }
        if(grepl('[Ll]',  substring(plat,1,1))) {
          instrument<-"LC-MS"
        }
      }
      if(is.null(instrument)) {
        stop("do not regonize instrument platform, please use either 'GC-MS' or 'LC-MS' " )
      }
      
      if(instrument == "LC-MS") platform<-LCMS
      if(instrument == "GC-MS") platform<-GCMS
      
      # rowstart<-grep(instrument, csv.in[,1])+1
      # rowend<-grep("MSlevs", csv.in[,1])
      # rowend<-rowend[which(rowend > rowstart)]
      # if(length(rowend)>1) {
      #   rowend<-rowend[which.min((rowend - rowstart))]
      # }
      # instrument <- data.frame('value' = csv.in[rowstart:rowend,2], row.names = csv.in[rowstart:rowend,1], stringsAsFactors = FALSE)
      # ExpDes <- list("design" = design, "instrument" = instrument)
      # 
      instrument <- platform
      
      suppressWarnings( instrument<-edit(instrument))
      # instrument <- data.frame('value' = csv.in[rowstart:rowend,2], row.names = csv.in[rowstart:rowend,1], stringsAsFactors = FALSE)
      ExpDes<-list("design" = design, "instrument" = instrument)
      
    }
    
  } else  {
    if(file.exists(csv)) { 
      csv.in <- read.csv(csv, header=TRUE, check.names=FALSE)
      design<-data.frame("value" = csv.in[3:7,2], row.names = csv.in[3:7,1], stringsAsFactors = FALSE)
      
      instrument <- NULL
      plat<-as.character(design[5,1])
      if(grepl("LC-MS", plat)) {
        instrument<-"LC-MS"
      }  
      if( grepl("GC-MS", plat )) {
        instrument<-"GC-MS"
      }  
      if(!grepl("LC-MS", plat ) & !grepl("GC-MS", plat ) ) {
        if(grepl('[Gg]',  substring(plat,1,1))) {
          instrument<-"GC-MS"
        }
        if(grepl('[Ll]',  substring(plat,1,1))) {
          instrument<-"LC-MS"
        }
      }
      if(is.null(instrument)) {
        stop("do not regonize instrument platform, please use either 'GC-MS' or 'LC-MS' " )
      }
      
      rowstart<-grep(instrument, csv.in[,1])+1
      rowend<-grep("MSlevs", csv.in[,1])
      rowend<-rowend[which(rowend > rowstart)]
      if(length(rowend)>1) {
        rowend<-rowend[which.min((rowend - rowstart))]
      }
      instrument <- data.frame('value' = csv.in[rowstart:rowend,2], row.names = csv.in[rowstart:rowend,1], stringsAsFactors = FALSE)
      ExpDes <- list("design" = design, "instrument" = instrument)
      
    }
  }
  
  return(ExpDes)
}

