#' rc.feature.normalize.quantile
#'
#' normalize data using quantile
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @return  ramclustR object with normalized data.
#' @export

rc.feature.normalize.quantile <- function(ramclustObj = NULL) {
    tmpnames1 <- dimnames(ramclustObj$MSdata)
    ramclustObj$MSdata <- t(preprocessCore::normalize.quantiles(t(ramclustObj$MSdata)))
    dimnames(ramclustObj$MSdata) <- tmpnames1

    if (!is.null(ramclustObj$MSMSdata)) {
        tmpnames2 <- dimnames(ramclustObj$MSMSdata)
        ramclustObj$MSMSdata <- t(preprocessCore::normalize.quantiles(t(ramclustObj$MSMSdata)))
        dimnames(ramclustObj$MSMSdata) <- tmpnames2
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
