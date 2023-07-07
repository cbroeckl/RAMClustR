normalize_quantiles <- function(x) {
    return(t(preprocessCore::normalize.quantiles(as.matrix(t(x)), keep.names = TRUE)))
}

#' rc.feature.normalize.quantile
#'
#' normalize data using quantile
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @return  ramclustR object with normalized data.
#' @export
rc.feature.normalize.quantile <- function(ramclustObj = NULL) {
    ramclustObj$MSdata <- normalize_quantiles(ramclustObj$MSdata)

    if (!is.null(ramclustObj$MSMSdata)) {
        ramclustObj$MSMSdata <- normalize_quantiles(ramclustObj$MSMSdata)
    }

    ## update msint and optionally msmsint
    global.min <- apply(cbind(ramclustObj$MSdata, ramclustObj$MSMSdata), 2, "min", na.rm = TRUE)

    ramclustObj$msint <- compute_wt_mean(
        ramclustObj$MSdata,
        global.min,
        ramclustObj$fmz,
        TRUE
    )

    if (!is.null(ramclustObj$MSMSdata)) {
        ramclustObj$msmsint <- compute_wt_mean(
            ramclustObj$MSMSdata,
            global.min,
            ramclustObj$fmz,
            TRUE
        )
    }

    ramclustObj$history$normalize.quantile <- paste(
        " Features were normalized using 'quantile' normalization."
    )

    return(ramclustObj)
}
