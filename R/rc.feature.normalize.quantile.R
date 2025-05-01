quantile_normalisation <- function(df) {
    df_rank <- apply(df, 2, rank, ties.method = "min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)

    index_to_mean <- function(my_index, my_mean) {
        return(my_mean[my_index])
    }

    df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
}

#' rc.feature.normalize.quantile
#'
#' normalize data using quantile
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @return  ramclustR object with normalized data.
#' @export
rc.feature.normalize.quantile <- function(ramclustObj = NULL) {
    ramclustObj$MSdata <- t(quantile_normalisation(t(ramclustObj$MSdata)))

    if (!is.null(ramclustObj$MSMSdata)) {
        ramclustObj$MSMSdata <- t(quantile_normalisation(t(ramclustObj$MSMSdata)))
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
