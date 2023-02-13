compute_wt_mean <- function(data, fmz, ensure.no.na) {
    global.min <- apply(data, 2, "min", na.rm = TRUE)

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

    if (!is.null(phenoData)) {
        phenoData <- read.csv(file = phenoData, header = TRUE, check.names = FALSE)
    }
    if (!is.null(idmsms)) {
        MSMSdata <- read.csv(file = idmsms, header = TRUE, check.names = FALSE)
    }
    if (is.null(idmsms)) {
        MSMSdata <- MSdata
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
    dimnames(data1)[[1]] <- MSdata[, sampNameCol]
    dimnames(data1)[[2]] <- names(MSdata[, featcol])
    data2 <- as.matrix(MSMSdata[, featcol])
    dimnames(data2)[[1]] <- MSMSdata[, sampNameCol]
    dimnames(data2)[[2]] <- names(MSMSdata[, featcol])
    if (!all(dimnames(data1)[[2]] == dimnames(data2)[[2]])) {
        stop("the feature names of your MS and idMSMS data are not identical")
    }

    if (!all(dimnames(data1)[[1]] == dimnames(data2)[[1]])) {
        stop("the order and names of your MS and idMSMS data sample names are not identical")
    }

    rtmz <- matrix(
        unlist(
            strsplit(dimnames(data1)[[2]], featdelim)
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
    data2 <- data2[, xcmsOrd]
    mzs <- mzs[xcmsOrd]
    times <- times[xcmsOrd]

    ## create empty hclust object to ultimately hold clustering data
    ramclustObj <- list()
    class(ramclustObj) <- "hclust"
    ramclustObj$merger <- vector(length = 0)
    ramclustObj$height <- vector(length = 0)
    ramclustObj$order <- vector(length = 0)
    ramclustObj$labels <- vector(length = 0)
    ramclustObj$method <- vector(length = 0)
    ramclustObj$call <- vector(length = 0)
    ramclustObj$dist.method <- vector(length = 0)
    ramclustObj$ExpDes <- ExpDes
    ramclustObj$history <- list()
    ramclustObj$MSdata <- data1
    ramclustObj$MSMSdata <- data2
    ramclustObj$frt <- times
    ramclustObj$fmz <- mzs
    ramclustObj$st <- st
    ramclustObj$history$input <- history
    ramclustObj$MSdata_raw <- ramclustObj$MSdata
    ramclustObj$MSMSdata_raw <- ramclustObj$MSMSdata
    ramclustObj$phenoData <- phenoData
    ramclustObj$featnames <- dimnames(data1)[[2]]
    ramclustObj$xcmsOrd <- xcmsOrd
    ramclustObj$msint <- compute_wt_mean(ramclustObj$MSdata, ramclustObj$fmz, ensure.no.na)
    ramclustObj$msmsint <- compute_wt_mean(ramclustObj$MSMSdata, ramclustObj$fmz, ensure.no.na)

    return(ramclustObj)
}
