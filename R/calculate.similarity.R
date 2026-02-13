#' @exportS3Method NULL
compute.start <- function(position, blocksize, numcols) {
  return(min((1 + ((position - 1) * blocksize)), numcols))
}

#' @exportS3Method NULL
compute.stop <- function(position, blocksize, numcols) {
  if ((position) * blocksize > numcols) {
    stop_pos <- numcols
  } else {
    stop_pos <- (position) * blocksize
  }
  return(stop_pos)
}

#' @exportS3Method NULL
calculate.similarity <- function(numcols,
                                 data1,
                                 data2,
                                 times,
                                 blocksize,
                                 mult,
                                 maxt,
                                 st,
                                 sr,
                                 rt.only.low.n,
                                 cor.method,
                                 cor.use) {
  ########
  # establish some constants for downstream processing
  vlength <- (numcols * (numcols - 1)) / 2
  nblocks <- ceiling(numcols / blocksize)

  ramclustObj <- vector(mode = "integer", length = vlength)
  block = 1
  column <- NULL
  
  for (row in 1:(nblocks)) {
    for (col in row:(nblocks)) {
      block <- block + 1
      
      start_col <- compute.start(row, blocksize, numcols)
      stop_col <- compute.stop(row, blocksize, numcols)
      
      start_row <- compute.start(col, blocksize, numcols)
      stop_row <- compute.stop(col, blocksize, numcols)
      
      if (start_col <= start_row) {
        mint <- min(abs(outer(range(times[start_row:stop_row]), 
                              range(times[start_col:stop_col]), 
                              FUN = "-")))
        if (mint <= maxt) {
          # RT similarity
          RT_sim <- round(exp(-(((abs(outer(times[start_row:stop_row], 
                                            times[start_col:stop_col], 
                                            FUN = "-"))
          )) ^ 2) / (2 * (st ^ 2))), digits = 20)
          
          if (nrow(data1) < 5 & rt.only.low.n) {
            # correlational similarity
            corr_sim <- matrix(data = 1,
                               nrow = length(start_row:stop_row),
                               ncol = length(start_col:stop_col)
            )
          } else {
            suppressWarnings(
              max_value <- pmax(
                cor(
                  data1[, start_row:stop_row], data1[, start_col:stop_col], method = cor.method, use = cor.use),
                cor(
                  data1[, start_row:stop_row], data2[, start_col:stop_col], method = cor.method, use = cor.use),
                cor(
                  data2[, start_row:stop_row], data2[, start_col:stop_col], method = cor.method, use = cor.use)
                , na.rm = TRUE
              )

            )
            if(any(is.na(max_value))) {
              max_value[is.na(max_value)] <- 0
            }
            # correlational similarity
            corr_sim <- round(exp(-((1 - max_value) ^ 2) / (2 * (sr ^ 2))), digits = 20)
          }
          
          # overall similarity
          sim_matrix <- 1 - (RT_sim * corr_sim)
          # remove nans and infs
          sim_matrix[which(is.nan(sim_matrix)) | 
                       which(is.na(sim_matrix)) | 
                       which(is.infinite(sim_matrix))] <- 1
        } else {
          # overall similarity
          sim_matrix <- matrix(data = 1,
                               nrow = length(start_row:stop_row),
                               ncol = length(start_col:stop_col)
          )
        }
      }
      
      # merge computed similarity matrices to single matrix (extend rows)
      if (is.null(column)) {
        column <- sim_matrix
      } else {
        column <- rbind(column, sim_matrix)
      }
    }
    
    # remove values above the diagonal and convert to vector by columns
    column <- column[which(row(column) - col(column) > 0)]
    
    if (exists("startv") == FALSE)
      startv <- 1
    stopv <- startv + length(column) - 1

    # assign obtained vector to result
    ramclustObj[startv:stopv] <- column
    startv <- stopv + 1
    
    # remove column to start next iteration with clean slate
    column <- NULL
  }
  gc()
  
  return(ramclustObj)
}
