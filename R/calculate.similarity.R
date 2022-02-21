calculate.similarity <- function(n,
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
                                 fftempdir) {
  ########
  # establish some constants for downstream processing
  vlength <- (n * (n - 1)) / 2
  nblocks <- floor(n / blocksize)
  
  ########
  # set off ff matrix system for holding data.
  # manages RAM demands a bit.
  ffmat <- ff::ff(vmode = "double",
                  dim = c(n, n),
                  initdata = 0) ##reset to 1 if necessary
  gc()
  
  ########
  # make list of all row and column blocks to evaluate
  eval1 <- expand.grid(0:nblocks, 0:nblocks)
  names(eval1) <- c("j", "k") #j for cols, k for rows
  eval1 <- eval1[which(eval1[, "j"] <= eval1[, "k"]), ] #upper triangle only
  bl <- nrow(eval1)
  cat(paste("calculating ramclustR similarity: nblocks = ", bl, '\n'))
  
  ########
  # Define the RCsim function used to calculate feature similarities on selected blocks of data
  RCsim <- function(bl)  {
    cat(bl, ' ')
    j <- eval1[bl, "j"]  #columns
    k <- eval1[bl, "k"]  #rows
    startc <- min((1 + (j * blocksize)), n)
    if ((j + 1) * blocksize > n) {
      stopc <- n
    } else {
      stopc <- (j + 1) * blocksize
    }
    startr <- min((1 + (k * blocksize)), n)
    if ((k + 1) * blocksize > n) {
      stopr <- n
    } else {
      stopr <- (k + 1) * blocksize
    }
    if (startc <= startr) {
      mint <- min(abs(outer(range(times[startr:stopr]), range(times[startc:stopc]), FUN = "-")))
      if (mint <= maxt) {
        temp1 <- round(exp(-(((abs(outer(times[startr:stopr], times[startc:stopc], FUN = "-"))
                               )) ^ 2) / (2 * (st ^ 2))),
                       digits = 20)
        
        if (nrow(data1) < 5 & rt.only.low.n) {
          temp2 <- matrix(data = 1,
                          nrow = length(startr:stopr),
                          ncol = length(startc:stopc)
                         )
        } else {
          temp2 <-
            round (exp(-((
              1 - (pmax(
                cor(data1[, startr:stopr], data1[, startc:stopc], method = cor.method),
                cor(data1[, startr:stopr], data2[, startc:stopc], method =
                      cor.method),
                cor(data2[, startr:stopr], data2[, startc:stopc], method =
                      cor.method)
              ))
            ) ^ 2) / (2 * (sr ^ 2))),
            
            digits = 20)
        }
        #ffcor[startr:stopr, startc:stopc]<-temp
        temp <- 1 - (temp1 * temp2)
        temp[which(is.nan(temp))] <- 1
        temp[which(is.na(temp))] <- 1
        temp[which(is.infinite(temp))] <- 1
        ffmat[startr:stopr, startc:stopc] <- temp
        rm(temp1)
        rm(temp2)
        rm(temp)
        gc()
      }
      if (mint > maxt) {
        ffmat[startr:stopr, startc:stopc] <- 1
      }
    }
    gc()
  }
  
  ########
  # Call the similarity scoring function
  system.time(sapply(1:bl, RCsim))
  
  b <- Sys.time()
  
  ########
  # Report progress and timing
  cat('\n',
      "RAMClust feature similarity matrix calculated and stored:",
      '\n')
  gc()
  
  
  ########
  # extract lower diagonal of ffmat as vector
  blocksize <- mult * round(blocksize ^ 2 / n)
  nblocks <- floor(n / blocksize)
  remaind <- n - (nblocks * blocksize)
  
  ########
  # create vector for storing dissimilarities
  ramclustObj <- vector(mode = "integer", length = vlength)
  
  ########
  # fill vector with dissimilarities
  for (k in 0:(nblocks)) {
    startc <- 1 + (k * blocksize)
    if ((k + 1) * blocksize > n) {
      stopc <- n
    } else {
      stopc <- (k + 1) * blocksize
    }
    temp <- ffmat[startc:nrow(ffmat), startc:stopc]
    temp <- temp[which(row(temp) - col(temp) > 0)]
    if (exists("startv") == FALSE)
      startv <- 1
    stopv <- startv + length(temp) - 1
    ramclustObj[startv:stopv] <- temp
    gc()
    startv <- stopv + 1
    rm(temp)
    gc()
  }
  rm(startv)
  gc()
  
  ########
  # cleanup
  close(ffmat)
  rm(ffmat)
  gc()
  if (!is.null(fftempdir)) {
    options("fftempdir" = origffdir)
  }
  
  return(ramclustObj)
}
