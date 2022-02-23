compute.start <- function(position, blocksize, numcols) {
  return(min((1 + (position * blocksize)), numcols))
}

compute.stop <- function(position, blocksize, numcols) {
  if ((position + 1) * blocksize > numcols) {
    stop_pos <- numcols
  } else {
    stop_pos <- (position + 1) * blocksize
  }
  return(stop_pos)
}

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
                                 cor.method) {
  ########
  # establish some constants for downstream processing
  vlength <- (numcols * (numcols - 1)) / 2
  nblocks <- floor(numcols / blocksize)
  
  cat(paste("Calculating ramclustR similarity using", sum((nblocks+1):1), "nblocks.\n"))
  
  ramclustObj <- vector(mode = "integer", length = vlength)
  block = 1
  
  for (row in 0:(nblocks)) {
    for (col in row:(nblocks)) {
      cat(block, " ")
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
            max_value <- pmax(cor(data1[, start_row:stop_row], data1[, start_col:stop_col], method = cor.method),
                              cor(data1[, start_row:stop_row], data2[, start_col:stop_col], method = cor.method),
                              cor(data2[, start_row:stop_row], data2[, start_col:stop_col], method = cor.method)
            )
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
      if (exists("column") == FALSE) {
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
    rm(column)
  }
  
  cat('\n RAMClust feature similarity matrix calculated and stored.\n')
  gc()
  
  return(ramclustObj)
}
