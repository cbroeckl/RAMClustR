#' compute_wt_mean
#'
#' compute weighted.mean intensity of feature in ms/msms level data
#' @param data feature in ms/msms level data
#' @param global.min minimum intensity in ms/msms level data
#' @param fmz feature retention time
#' @param ensure.no.na logical: if TRUE, any 'NA' values in msint and/or msmsint are replaced with numerical values based on 10 percent of feature min plus noise.  Used to ensure that spectra are not written with NA values.
#' @return weighted.mean intensity of feature in ms/msms level data
#' @export
compute_wt_mean <- function(data, global.min, fmz, ensure.no.na) {
    wt_mean_int <- rep(0, length(fmz))
    for (i in 1:ncol(data)) {
        wt_mean_int[i] <- weighted.mean(data[, i], data[, i], na.rm = TRUE)
    }
    if (any(is.na(wt_mean_int)) & ensure.no.na) {
        rp <- which(is.na(wt_mean_int))
        r <- global.min[rp]
        r <- abs(jitter(r, factor = 0.01 * r))
        wt_mean_int[rp] <- r
    }

    return(wt_mean_int)
}

#' rc.get.csv.data
#'
#' extractor for csv objects in preparation for normalization and clustering
#'
#' @param csv filepath: csv input. Features as columns, rows as samples. Column header mz_rt
#' @param phenoData character: character string in 'taglocation' to designate files as either MS / DIA(MSe, MSall, AIF, etc) e.g. "01.mzML"
#' @param idmsms filepath: optional idMSMS / MSe csv data.  same dim and names as ms required
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param sampNameCol integer: which column from the csv file contains sample names?
#' @param st numeric: sigma t - time similarity decay value
#' @param timepos integer: which position in delimited column header represents the retention time
#' @param featdelim character: how feature mz and rt are delimited in csv import column header e.g. ="-"
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
#' @author Corey Broeckling
#' @export
#' @examples
#' ## Choose csv input file. Features as columns, rows as samples
#' ## Choose csv input file phenoData 
#' filename <- system.file("extdata", "peaks.csv", package = "RAMClustR", mustWork = TRUE)
#' phenoData <- system.file("extdata", "phenoData.csv", package = "RAMClustR", mustWork = TRUE)
#'
#' ramclustobj <- rc.get.csv.data(csv = filename, phenoData = phenoData, st = 5)
#'

rc.get.csv.data <- function(csv = NULL,
                            phenoData = NULL,
                            idmsms = NULL,
                            ExpDes = NULL,
                            sampNameCol = 1,
                            st = NULL,
                            timepos = 2,
                            featdelim = "_",
                            ensure.no.na = TRUE) {
    if (is.null(csv)) {
        stop("Please provide an MS dataset with features as columns
             (__)one column may contain sample names, defined by sampNameCol)")
    }
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
        " Feature data was input as .csv files"
    )

    MSdata <- read.csv(file = csv, header = TRUE, check.names = FALSE)
    data2 <- NULL

    if (!is.null(phenoData)) {
        phenoData <- read.csv(file = phenoData, header = TRUE, check.names = FALSE)
    }
    if (!is.null(idmsms)) {
        MSMSdata <- read.csv(file = idmsms, header = TRUE, check.names = FALSE)
    }
    if (is.null(sampNameCol)) {
        featcol <- 1:ncol(MSdata)
    } else {
        featcol <- setdiff(1:(ncol(MSdata)), sampNameCol)
    }
    if (is.null(sampNameCol)) {
        featcol <- 1:ncol(MSdata)
    } else {
        featcol <- setdiff(1:(ncol(MSdata)), sampNameCol)
    }
    sampnames <- MSdata[, sampNameCol]
    data1 <- as.matrix(MSdata[, featcol])
    rownames(data1) <- MSdata[, sampNameCol]
    colnames(data1) <- names(MSdata[, featcol])

    if (!is.null(idmsms)) {
        data2 <- as.matrix(MSMSdata[, featcol])
        rownames(data2) <- MSMSdata[, sampNameCol]
        colnames(data2) <- names(MSMSdata[, featcol])

        if (!all(colnames(data1) == colnames(data2))) {
            stop("the feature names of your MS and idMSMS data are not identical")
        }

        if (!all(rownames(data1) == rownames(data2))) {
            stop("the order and names of your MS and idMSMS data sample names are not identical")
        }
    }
    rtmz <- matrix(
        unlist(
            strsplit(colnames(data1), featdelim)
        ),
        byrow = TRUE, ncol = 2
    )
    times <- as.numeric(rtmz[, timepos])
    mzs <- as.numeric(rtmz[, which(c(1:2) != timepos)])

    if (any(is.na(times))) {
        stop("column(s) ", which(is.na(times)), " rt cannot be made numeric. ")
    }

    if (any(is.na(mzs))) {
        stop("column(s) ", which(is.na(times)), " mz cannot be made numeric. ")
    }

    rm(rtmz)

    # data organization and parsing
    # sort rt vector and data by retention time
    xcmsOrd <- order(times)
    data1 <- data1[, xcmsOrd]
    mzs <- mzs[xcmsOrd]
    times <- times[xcmsOrd]

    if (!is.null(data2)) {
        data2 <- data2[, xcmsOrd]
    }

    ramclustObj <- create_ramclustObj(
        ExpDes = ExpDes,
        MSdata = data1,
        MSMSdata = data2,
        frt = times,
        fmz = mzs,
        st = st,
        input_history = history,
        phenoData = phenoData,
        feature_names = colnames(data1),
        xcmsOrd = xcmsOrd,
        sample_names = rownames(data1)
    )

    return(ramclustObj)
}
