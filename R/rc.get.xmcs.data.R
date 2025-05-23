#' Extract Sample Metadata from an xcms Object
#'
#' This function retrieves sample metadata from a provided xcms object.
#'
#' @param xcmsObj An xcms object containing chromatographic and metadata information.
#'                This object is typically generated during preprocessing of 
#'                metabolomics data using the xcms package.
#'
#' @return A data frame or list containing the sample metadata extracted from the xcms object.
#'         The structure and content of the returned metadata depend on the input xcms object.
#'
#' @details This function is designed to facilitate the extraction of metadata
#'          associated with samples in an xcms object. The metadata may include
#'          information such as sample names, injection order, batch information,
#'          and other experimental details.
#'
#' @concept xcms
#' @export
getSampleMetadata <- function(xcmsObj) {
  if (inherits(xcmsObj@phenoData, "NAnnotatedDataFrame")) {
    return(xcmsObj@phenoData@data)
  } else {
    stop("Unsupported phenoData class")
  }
}

#' rc.get.xcms.data
#'
#' extractor for xcms objects in preparation for normalization and clustering
#'
#' @param xcmsObj xcmsObject: containing grouped feature data for clustering by ramclustR
#' @param MStag character: character string in 'taglocation' to designate files as either MS / DIA(MSe, MSall, AIF, etc) e.g. "01.mzML"
#' @param MSMStag character: character string in 'taglocation' to designate files as either MS / DIA(MSe, MSall, AIF, etc) e.g. "02.mzML"
#' @param taglocation character: "filepaths" by default, "phenoData[,1]" is another option. refers to xcms slot
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param mzdec integer: number of decimal places for storing m/z values
#' @param ensure.no.na logical: if TRUE, any 'NA' values in msint and/or msmsint are replaced with numerical values based on 10 percent of feature min plus noise.  Used to ensure that spectra are not written with NA values.
#' @details This function creates a ramclustObj which will be used as input for clustering.
#' @return  an empty ramclustR object.  this object is formatted as an hclust object with additional slots for holding feature and compound data. details on these found below.
#' @return   $frt: feature retention time, in whatever units were fed in (xcms uses seconds, by default)
#' @return   $fmz: feature retention time, reported in number of decimal points selected in ramclustR function
#' @return   $ExpDes: the experimental design object used when running ramclustR.  List of two dataframes.
#' @return   $MSdata:  the MSdataset provided by either xcms or csv input
#' @return   $MSMSdata: the (optional) DIA(MSe, MSall, AIF etc) dataset provided be either xcms or csv input
#' @return   $xcmsOrd: original xcms order of features, for back-referencing when necessary
#' @return   $msint: weighted.mean intensity of feature in ms level data
#' @return   $msmsint:weighted.mean intensity of feature in msms level data
#'
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export
rc.get.xcms.data <- function(xcmsObj = NULL,
                             taglocation = "filepaths",
                             MStag = NULL,
                             MSMStag = NULL,
                             ExpDes = NULL,
                             mzdec = 3,
                             ensure.no.na = TRUE) {
  MSMSdata <- NULL
  ########
  # If experimental design is NULL:
  if (is.null(ExpDes)) {
    ExpDes <- defineExperiment(force.skip = TRUE)
    warning(
      "\n", "you failed to define your experimental descriptor using 'defineExperiment()'", "\n",
      "RAMClustR must now guess at what you are trying to do ", "\n",
      "and your exported spectra will be labelled incorrectly"
    )
  }

  params <- c(
    "taglocation" = "filepaths",
    "MStag" = MStag,
    "MSMStag" = MSMStag,
    "mzdec" = mzdec,
    "ensure.no.na" = ensure.no.na
  )
  ## add xcms processing history narrative here

  ## check xcms object presence
  if (is.null(xcmsObj)) {
    stop("please supply an xcms object as input", "\n")
  }

  ## check xcms format
  newXCMS <- FALSE
  if (!is.null(xcmsObj)) {
    if (!requireNamespace("xcms", quietly = TRUE)) {
      stop(
        "The use of this function requires package 'xcms'. Please ",
        "install with 'Biobase::install(\"xcms\")'"
      )
    }

    OK <- FALSE

    if (inherits(xcmsObj, "xcmsSet")) {
      OK <- TRUE
    }
    if (inherits(xcmsObj, "XCMSnExp")) {
      OK <- TRUE
      newXCMS <- TRUE
    }
    if (!OK) {
      stop("xcmsObj must reference an object generated by XCMS, of class 'xcmsSet'")
    }
  }

  ## check to see that we can find the MS2 data files

  ##  check to see if we have MS2 data, if not, error or switch to ms1
  if (!is.null(MSMStag)) {
    if (is.null(taglocation)) {
      stop("you must specify the the MStag, MSMStag, and the taglocations")
    }
    if (!any(grepl(taglocation, c("filepaths", "pheno")))) {
      stop("taglocation needs to be one of 'filepaths' or 'pheno'", "\n")
    }
  }

  if (!is.null(MSMStag)) {
    if (taglocation == "filepaths") {
      if (!newXCMS) {
        nfiles <- length(xcmsObj@filepaths)
        if (!is.null(MSMStag)) {
          msmsfiles <- grep(MSMStag, xcmsObj@filepaths, ignore.case = TRUE)
        }
      }
      if (newXCMS) {
        nfiles <- length(xcmsObj@processingData@files)
        if (!is.null(MSMStag)) {
          msmsfiles <- grep(MSMStag, xcmsObj@processingData@files, ignore.case = TRUE)
        }
      }
    }
    if (taglocation == "pheno") {
      if (!newXCMS) {
        nfiles <- length(xcmsObj@phenoData)
        if (!is.null(MSMStag)) {
          msmsfiles <- grep(MSMStag, row.names(xcmsObj@phenoData), ignore.case = TRUE)
        }
      }
      if (newXCMS) {
        if (inherits(xcmsObj, "XcmsExperiment")) {  
          nfiles <- length(xcmsObj)  
          if (!is.null(MSMStag)) {  
              msmsfiles <- grep(MSMStag, MsExperiment::sampleData(xcmsObj)[, 1L], ignore.case = TRUE)  
          }  
        }  
        if (inherits(xcmsObj, "XCMSnExp")) {  
          nfiles <- nrow(getSampleMetadata(xcmsObj))
          if (!is.null(MSMStag)) {  
            msmsfiles <- grep(MSMStag, getSampleMetadata(xcmsObj)[, 1L], ignore.case = TRUE)  
          }  
        }  
      }
      
    }
  } else {
    if (!newXCMS) {
      nfiles <- nrow(xcmsObj@phenoData)
    }
    if (newXCMS) {
      if (inherits(xcmsObj, "XCMSnExp")) {  
          nfiles <- nrow(getSampleMetadata(xcmsObj))  
      }  
      if (inherits(xcmsObj, "XcmsExperiment")) {  
          nfiles <- length(xcmsObj)  
      }
    }
  }

  if (!newXCMS) st <- round(median(xcmsObj@peaks[, "rtmax"] - xcmsObj@peaks[, "rtmin"]) / 2, digits = 2)
  if (newXCMS) st <- round(median(xcmsObj@msFeatureData$chromPeaks[, "rtmax"] - xcmsObj@msFeatureData$chromPeaks[, "rtmin"]) / 2, digits = 2)

  if (is.null(MSMStag)) {
    msfiles <- 1:nfiles
    mslev <- 1
  }
  if (!is.null(MSMStag)) {
    if (length(msmsfiles) == 0) {
      stop("no MSMS files found", "\n")
    } else {
      msfiles <- (1:nfiles)[-msmsfiles]
      mslev <- 2
    }
  }


  ## check to make sure same number of MS and MSMS files
  if (mslev == 2) {
    if (length(msfiles) != length(msmsfiles)) {
      stop("detected ", length(msfiles), " ms files and ", length(msmsfiles), " msms files - ", "\n", "       number of MSMS files MUST be identical to number of MS files")
    }
  }


  if (!newXCMS) {
    data <- t(xcms::groupval(xcmsObj, value = "into"))
  }
  if (newXCMS) {
    data <- t(xcms::featureValues(xcmsObj, value = "into"))
  }

  if (length(msfiles) == 0) {
    stop("no msfiles recognized")
  }

  ## get phenotype file name associations for storage in new RC object
  if (newXCMS) {
    filepaths <- MSnbase::fileNames(xcmsObj)
    filenames <- basename(filepaths)
    phenotype <- getSampleMetadata(xcmsObj)
    phenotype <- data.frame(sample.names = phenotype, filenames, filepaths)
    if (mslev == 2) {
      phenotype <- phenotype[1:(nrow(phenotype) / 2), ]
    }
  } else {
    filepaths <- xcmsObj@filepaths
    filenames <- basename(filepaths)
    phenotype <- xcmsObj@phenoData[, 1]
    phenotype <- data.frame(sample.names = phenotype, filenames, filepaths)
    if (mslev == 2) {
      phenotype <- phenotype[1:(nrow(phenotype) / 2), ]
    }
  }

  history <- {
    paste0(
      "RAMClustR version ", utils::packageDescription("RAMClustR")$Version, " in ", R.Version()$version.string,
      ") was used to normalize, filter, and group features into spectra.",
      "XCMS (Smith 2006)(Tautenhahn 2008) output data was transferred to a ramclustR object using the rc.get.xcms.data function. ",
      "Feature data was extracted using the xcms ", if (newXCMS) {
        "featureValues"
      } else {
        "groupval"
      },
      " function."
    )
  }

  ## process data
  # get feature RTs
  if (!newXCMS) {
    times <- round(xcmsObj@groups[, "rtmed"], digits = 3)
  }
  if (newXCMS) {
    times <- round(xcmsObj@msFeatureData$featureDefinitions$rtmed, digits = 3)
  }
  if (any(is.na(times))) {
    do <- which(is.na(times))
    for (x in 1:length(do)) {
      if (!newXCMS) {
        times[do[x]] <- as.numeric((xcmsObj@groups[do[x], "rtmin"] + xcmsObj@groups[do[x], "rtmax"]) / 2)
      }
      if (newXCMS) {
        times[do[x]] <- as.numeric((xcmsObj@msFeatureData$featureDefinitions$rtmin[do[x]] + xcmsObj@msFeatureData$featureDefinitions$rtmax[do[x]]) / 2)
      }
    }
  }

  # get feature MZs
  if (!newXCMS) {
    mzs <- round(xcmsObj@groups[, "mzmed"], digits = mzdec)
  }
  if (newXCMS) {
    mzs <- round(xcmsObj@msFeatureData$featureDefinitions$mzmed, digits = mzdec)
  }

  featnames <- rownames(xcms::featureDefinitions(xcmsObj))
  
  # reorder feature data by RT, record original xcmsOrder
  xcmsOrd <- order(times)
  data <- data[, xcmsOrd]
  mzs <- mzs[xcmsOrd]
  times <- times[xcmsOrd]
  featnames <- featnames[xcmsOrd]
  dimnames(data)[[2]] <- featnames
  dimnames(data)[[1]] <- filenames

  if (mslev == 2) {
    MSMSdata <- data[msmsfiles, ]
  }

  ramclustObj <- create_ramclustObj(
    ExpDes = ExpDes,
    MSdata = data[msfiles, ],
    MSMSdata = MSMSdata,
    frt = times,
    fmz = mzs,
    st = st,
    input_history = history,
    phenoData = phenotype,
    feature_names = featnames,
    xcmsOrd = xcmsOrd,
    sample_names = filenames,
    ensure.no.na = ensure.no.na
  )

  if (is.null(ramclustObj$params)) {
    ramclustObj$params <- list()
  }
  ramclustObj$params$rc.get.xcms.data <- params

  return(ramclustObj)
}
