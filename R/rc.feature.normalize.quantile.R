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

    ramclustObj$history$normalize.quantile <- paste(
        " Features were normalized using 'quantile' normalization."
    )

    return(ramclustObj)
}
