#' create_ramclustObj
#'
#' create ramclustr Object
#'
#' @param merger
#' @param height
#' @param order integer vector with length equal to number of injections
#' @param labels feature names
#' @param method method used e.g "RAMClustR"
#' @param call
#' @param dist.method
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param input_history input history
#' @param MSdata dataframe containing MS Data
#' @param MSMSdata dataframe containing MSMS Data
#' @param frt feature retention time, in whatever units were fed in
#' @param fmz feature retention time
#' @param st numeric: sigma t - time similarity decay value
#' @param MSdata_raw dataframe containing MS Data
#' @param MSMSdata_raw dataframe containing MSMS Data
#' @param phenoData dataframe containing phenoData
#' @param feature_names feature names extracted from the data
#' @param xcmsOrd original xcms order of features, for back-referencing when necessary
#' @param ensure.no.na logical: if TRUE, any 'NA' values in msint and/or msmsint are replaced with numerical values based on 10 percent of feature min plus noise.  Used to ensure that spectra are not written with NA values.
#' @return  an ramclustR object. this object is formatted as an hclust object with additional slots for holding feature and compound data.
#' @export

create_ramclustObj <- function(merger = vector(length = 0),
                               height = vector(length = 0),
                               order = vector(length = 0),
                               labels = vector(length = 0),
                               method = vector(length = 0),
                               call = vector(length = 0),
                               dist.method = vector(length = 0),
                               ExpDes = NULL,
                               input_history = NULL,
                               MSdata = NULL,
                               MSMSdata = NULL,
                               frt = NULL,
                               fmz = NULL,
                               st = NULL,
                               MSdata_raw = NULL,
                               MSMSdata_raw = NULL,
                               phenoData = NULL,
                               feature_names = NULL,
                               xcmsOrd = NULL,
                               ensure.no.na = TRUE) {
    ## create empty hclust object to ultimately hold clustering data
    ramclustObj <- list()
    class(ramclustObj) <- "hclust"
    ramclustObj$merger <- merger
    ramclustObj$height <- height
    ramclustObj$order <- order
    ramclustObj$labels <- labels
    ramclustObj$method <- method
    ramclustObj$call <- call
    ramclustObj$dist.method <- dist.method
    ramclustObj$ExpDes <- ExpDes
    ramclustObj$history <- list()
    ramclustObj$MSdata <- MSdata
    ramclustObj$MSMSdata <- MSMSdata
    ramclustObj$frt <- frt
    ramclustObj$fmz <- fmz
    ramclustObj$st <- st
    ramclustObj$history$input <- input_history
    ramclustObj$MSdata_raw <- MSdata_raw
    ramclustObj$MSMSdata_raw <- MSMSdata_raw
    ramclustObj$phenoData <- phenoData
    ramclustObj$featnames <- feature_names
    ramclustObj$xcmsOrd <- xcmsOrd

    global.min <- apply(cbind(ramclustObj$MSdata, ramclustObj$MSMSdata), 2, "min", na.rm = TRUE)

    ramclustObj$msint <- compute_wt_mean(ramclustObj$MSdata, global.min, ramclustObj$fmz, ensure.no.na)
    if (!is.null(MSMSdata)) {
        ramclustObj$msmsint <- compute_wt_mean(ramclustObj$MSMSdata, global.min, ramclustObj$fmz, ensure.no.na)
    }
    return(ramclustObj)
}

#' checks
#'
#' check if MS data contains mz and rt, and if MSMS data is present feature names and sample names are identical
#'
#' @param ms1_featureDefinitions dataframe with metadata with columns: mz, rt, feature names containing MS data
#' @param ms1_featureValues dataframe with rownames = sample names, colnames = feature names containing MS data
#' @param ms2_featureValues dataframe with rownames = sample names, colnames = feature names containing MSMS data
#' @param feature_names feature names extracted from the data
#' @export

checks <- function(ms1_featureDefinitions = NULL,
                   ms1_featureValues = NULL,
                   ms2_featureValues = NULL,
                   feature_names = NULL) {
    if (is.null(ms1_featureDefinitions$mz) && is.null(ms1_featureDefinitions$rt)) {
        stop("Please provide feature definition data frame which contain mz and rt")
    }

    if (!any(feature_names == colnames(ms1_featureValues))) {
        stop("ms1_featureValues column names are not equal to feature names of ms1_featureDefinitions")
    }

    if (!is.null(ms2_featureValues)) {
        if (!all(dimnames(ms1_featureValues)[[2]] == dimnames(ms2_featureValues)[[2]])) {
            stop("the feature names of your MS and idMSMS data are not identical")
        }

        if (!all(dimnames(ms1_featureValues)[[1]] == dimnames(ms2_featureValues)[[1]])) {
            stop("the order and names of your MS and idMSMS data sample names are not identical")
        }
    }
}

#' rc.get.df.data
#'
#' @description
#'
#' @param ms1_featureDefinitions dataframe with metadata with columns: mz, rt, feature names containing MS data
#' @param ms1_featureValues dataframe with rownames = sample names, colnames = feature names containing MS data
#' @param ms2_featureDefinitions dataframe with metadata with columns: mz, rt, feature names containing MSMS data
#' @param ms2_featureValues dataframe with rownames = sample names, colnames = feature names containing MSMS data
#' @param phenoData dataframe containing phenoData
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param featureNamesColumnIndex integer: which column in the dataframe contains feature names?
#' @param st numeric: sigma t - time similarity decay value
#' @param ensure.no.na logical: if TRUE, any 'NA' values in msint and/or msmsint are replaced with numerical values based on 10 percent of feature min plus noise.  Used to ensure that spectra are not written with NA values.
#' @details This function creates a ramclustObj which will be used as input for clustering.
#' @return  an empty ramclustR object.  this object is formatted as an hclust object with additional slots for holding feature and compound data. details on these found below.
#' @return   $frt: feature retention time, in whatever units were fed in
#' @return   $fmz: feature retention time, reported in number of decimal points selected in ramclustR function
#' @return   $ExpDes: the experimental design object used when running ramclustR.  List of two dataframes.
#' @return   $MSdata:  the MSdataset provided by either xcms or csv input
#' @return   $MSMSdata: the (optional) DIA(MSe, MSall, AIF etc) dataset
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
#' @author Zargham Ahmad, Helge Hecht, Corey Broeckling
#' @export

rc.get.df.data <- function(ms1_featureDefinitions = NULL,
                           ms1_featureValues = NULL,
                           ms2_featureDefinitions = NULL,
                           ms2_featureValues = NULL,
                           phenoData = NULL,
                           ExpDes = NULL,
                           featureNamesColumnIndex = 1,
                           st = NULL,
                           ensure.no.na = TRUE) {
    if (is.null(st)) {
        stop("please specify st:
      a recommended starting point is half the value of
      your average chromatographic peak width at half max (seconds)")
    }

    # If experimental design is NULL:
    if (is.null(ExpDes)) {
        ExpDes <- defineExperiment(force.skip = TRUE)
        warning(
            "\n", "you failed to define your experimental descriptor using 'defineExperiment()'", "\n",
            "RAMClustR must now guess at what you are trying to do ", "\n",
            "and your exported spectra will be labelled incorrectly"
        )
    }

    history <- paste("Raw mass spectrometry data were processed using an R based workflow for feature detection, retention time alignment, feature grouping, peak filling, feature clustering.")

    history <- paste(
        history,
        " Feature data was input as data frame"
    )

    feature_names <- ms1_featureDefinitions[[featureNamesColumnIndex]]

    checks(
        ms1_featureDefinitions,
        ms1_featureValues,
        ms2_featureValues,
        feature_names
    )

    mzs <- ms1_featureDefinitions$mz
    times <- ms1_featureDefinitions$rt
    sample_names <- rownames(ms1_featureValues)

    # data organization and parsing
    # sort data by retention time
    xcmsOrd <- order(times)
    ms1_featureValues <- ms1_featureValues[, xcmsOrd]
    mzs <- mzs[xcmsOrd]
    times <- times[xcmsOrd]

    if (!is.null(ms2_featureValues)) {
        ms2_featureValues <- ms2_featureValues[, xcmsOrd]
    }

    ramclustObj <- create_ramclustObj(
        ExpDes,
        MSdata = ms1_featureValues,
        MSMSdata = ms2_featureValues,
        frt = times,
        fmz = mzs,
        st = st,
        input_history = history,
        MSdata_raw = ms1_featureValues,
        MSMSdata_raw = ms2_featureValues,
        phenoData = phenoData,
        feature_names = feature_names,
        xcmsOrd = xcmsOrd
    )

    return(ramclustObj)
}
