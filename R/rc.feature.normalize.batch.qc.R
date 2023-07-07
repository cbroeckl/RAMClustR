#' normalized_data_batch_qc
#'
#' normalize data using batch.qc
#' @param data feature in ms/msms level data
#' @param batch integer vector with length equal to number of injections in xset or csv file or dataframe
#' @param order integer vector with length equal to number of injections in xset or csv file or dataframe
#' @param qc logical vector with length equal to number of injections in xset or csv file or dataframe
#' @param qc.inj.range integer: how many injections around each injection are to be scanned for presence of QC samples when using batch.qc normalization?  A good rule of thumb is between 1 and 3 times the typical injection span between QC injections.  i.e. if you inject QC ever 7 samples, set this to between 7 and 21.  smaller values provide more local precision but make normalization sensitive to individual poor outliers (though these are first removed using the boxplot function outlier detection), while wider values provide less local precision in normalization but better stability to individual peak areas.
#' @param output.plot logical set to TRUE to store plots
#' @return normalized data.

normalized_data_batch_qc <- function(data = NULL,
                                     batch = NULL,
                                     order = NULL,
                                     qc = NULL,
                                     qc.inj.range = 20,
                                     output.plot = FALSE) {
  max.ratio <- 4
  for (z in 1:ncol(data)) {
    tmp <- data[, z]
    featmed <- mean(tmp[qc])
    tmpn <- tmp

    for (i in unique(batch)) {
      do <- which(batch == i)
      doqc <- which(batch == i & qc)

      sds <- 1.96
      lcl <- mean(tmp[doqc]) - (sds * sd(tmp[doqc]))
      ucl <- mean(tmp[doqc]) + (sds * sd(tmp[doqc]))
      keep <- which(tmp[doqc] > lcl & tmp[doqc] < ucl)
      if (length(keep) > 0) doqc <- doqc[keep]

      batchmed <- mean(tmpn[doqc])
      f <- batchmed / featmed

      if (is.na(f)) next
      if (abs(log2(f)) > max.ratio) {
        if (f > 1) {
          f <- max.ratio
        }
        if (f < 1) {
          f <- 1 / max.ratio
        }
      }

      tmpn[do] <- tmp[do] / f
      tmpnqc <- tmpn

      for (x in do) {
        batchmed <- mean(tmpn[doqc])
        space <- abs(x - doqc)
        use <- which(space <= qc.inj.range)
        wts <- 1 / space[use]
        wts <- wts / sum(wts)
        localmed <- weighted.mean(tmpnqc[doqc[use]], weights = wts)

        if (is.na(localmed)) {
          for (y in 2:5) {
            use <- which(space <= qc.inj.range * y)
            wts <- 1 / space[use]
            wts <- wts / sum(wts)
            localmed <- weighted.mean(tmpnqc[doqc[use]], weights = wts)
            if (!is.na(localmed)) break
          }
        }

        if (is.na(localmed)) {
          localmed <- batchmed
        }
        f <- localmed / batchmed

        if (is.na(f)) next

        if (abs(log2(f)) > max.ratio) {
          if (f > 1) {
            f <- max.ratio
          }
          if (f < 1) {
            f <- 1 / max.ratio
          }
        }

        tmpn[x] <- tmpnqc[x] / f
        rm(localmed)
        rm(f)
      }
    }

    data[, z] <- tmpn
    if (output.plot) {
      par(mfrow = c(1, 2))
      plot(tmp,
        col = batch, cex = (qc + 1) / 2, ylim = c(0.9, 1.11) * range(tmp),
        main = paste(
          "all:", round(sd(tmp) / mean(tmp), digits = 2), "\n",
          "qc:", round(sd(tmp[qc]) / mean(tmp[qc]), digits = 2)
        )
      )
      plot(tmpn,
        col = batch, pch = 19, cex = (qc + 1) / 2, , ylim = c(0.9, 1.11) * range(tmp),
        main = paste(
          "all:", round(sd(tmpn) / mean(tmpn), digits = 2), "\n",
          "qc:", round(sd(tmpn[qc]) / mean(tmpn[qc]), digits = 2)
        )
      )
    }
    rm(tmp)
    rm(tmpn)
    if (exists("tmpnqc")) {
      rm(tmpnqc)
    }
    if (z %% 100 == 0) {
      gc()
    }
  }

  return(data)
}


#' order_datasets
#'
#' order the datasets first by batch and run order
#' @param data feature in ms/msms level data
#' @param batch integer vector with length equal to number of injections in xset or csv file or dataframe
#' @param order integer vector with length equal to number of injections in xset or csv file or dataframe
#' @param qc logical vector with length equal to number of injections in xset or csv file or dataframe
#' @return ordered feature in ms/msms level data, order, batch, qc

order_datasets <- function(order = NULL,
                           batch = NULL,
                           qc = NULL,
                           data = NULL) {
  ndf <- data.frame(batch, order, qc)
  new.ord <- order(ndf$order)
  ndf <- ndf[new.ord, ]
  data <- data[new.ord, ]
  new.ord <- order(ndf$batch)
  ndf <- ndf[new.ord, ]
  data <- data[new.ord, ]
  batch <- ndf[, "batch"]
  qc <- ndf[, "qc"]
  order <- ndf[, "order"]

  return(list(
    data = data,
    batch = batch,
    qc = qc,
    order = order
  ))
}

#' rc.feature.normalize.batch.qc
#'
#' normalize data using batch.qc
#' @param order integer vector with length equal to number of injections in xset or csv file or dataframe
#' @param batch integer vector with length equal to number of injections in xset or csv file or dataframe
#' @param qc logical vector with length equal to number of injections in xset or csv file or dataframe
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param qc.inj.range integer: how many injections around each injection are to be scanned for presence of QC samples when using batch.qc normalization?  A good rule of thumb is between 1 and 3 times the typical injection span between QC injections.  i.e. if you inject QC ever 7 samples, set this to between 7 and 21.  smaller values provide more local precision but make normalization sensitive to individual poor outliers (though these are first removed using the boxplot function outlier detection), while wider values provide less local precision in normalization but better stability to individual peak areas.
#' @param output.plot logical set to TRUE to store plots
#' @return  ramclustR object with normalized data.
#' @export

rc.feature.normalize.batch.qc <- function(order = NULL,
                                          batch = NULL,
                                          qc = NULL,
                                          ramclustObj = NULL,
                                          qc.inj.range = 20,
                                          output.plot = FALSE) {
  if (!all.equal(length(batch), length(qc), length(order), nrow(ramclustObj$MSdata))) {
    stop(
      "all lengths must be identical and are not: ", "\n",
      "  length(batch) = ", length(batch), "\n",
      "  length(order) = ", length(order), "\n",
      "  length(qc) = ", length(qc), "\n",
      "  number of injections = ", nrow(ramclustObj$MSdata), "\n"
    )
  }

  order2 <- order
  batch2 <- batch
  qc2 <- qc

  ## for batch.qc method, we will order the datasets first by batch and run order
  # backup <- data1

  ordered_data <- order_datasets(
    data = ramclustObj$MSdata,
    order = order,
    batch = batch,
    qc = qc
  )

  order <- ordered_data$order
  batch <- ordered_data$batch
  qc <- ordered_data$qc
  ramclustObj$MSdata <- ordered_data$data

  if (!is.null(ramclustObj$MSMSdata)) {
    ordered_data <- order_datasets(
      data = ramclustObj$MSMSdata,
      order = order2,
      batch = batch2,
      qc = qc2
    )

    order2 <- ordered_data$order
    batch2 <- ordered_data$batch
    qc2 <- ordered_data$qc
    ramclustObj$MSMSdata <- ordered_data$data
  }

  ramclustObj$history$normalize.batch.qc <- paste(
    " Features were normalized to nearby QC samples on a feature-by-feature basis using the 'batch.qc' option ",
    "with qc.inj.range = ", qc.inj.range, ".",
    " QC normalization was applied to ", length(order), " injections in ", length(unique(batch)),
    " bactches, and normization was based on ", length(which(qc)), " recognized QC samples.",
    sep = ""
  )

  if (output.plot) {
    pdf(file = "norm.plots.pdf", height = 4, width = 9)
  }

  ramclustObj$MSdata <- normalized_data_batch_qc(
    data = ramclustObj$MSdata,
    batch = batch,
    order = order,
    qc = qc,
    qc.inj.range = qc.inj.range,
    output.plot = output.plot
  )

  if (!is.null(ramclustObj$MSMSdata)) {
    ramclustObj$MSMSdata <- normalized_data_batch_qc(
      data = ramclustObj$MSMSdata,
      batch = batch2,
      order = order2,
      qc = qc2,
      qc.inj.range = qc.inj.range
    )
  }

  if (output.plot) {
    dev.off()
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

  return(ramclustObj)
}
