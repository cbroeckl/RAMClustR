sonar.calibrate <- function(
    raw.sonar.data = d.spectra,
    nominal.mass.range = 10,
    ms.level = 1,
    mz.range = c(50,1200)
    
    
) {
  
  ## CHECKS
  # d.spectra <- Spectra::Spectra("C:/Users/cbroeckl/Documents/temp/20220124_sonarcaltest_025.mzML")
  ## check that object contains appropriate data
  ## should be a spectra object
  ## probably can also make an xcmsRaw object work, if need be.
  if(!class(raw.sonar.data)[1] == "Spectra") {
    stop("raw.sonar.data needs to be an object of class 'Spectra")
  }
  ## check that mz.range is numeric of length = 2
  if(!is.numeric(mz.range) | !length(mz.range == 2)) {
    stop("mz.range must be a numeric vector of length = 2, i.e. 'c(50,1200)'")
  }
  
  ## subset object to isolate MS2 data
  ms <- Spectra::filterMsLevel(d.spectra, ms.level)
  
  ## assign bins to MS2 data
  bins <- rep(1:200, length(ms)/200)
  unique.bins <- rep(1:200)
  
  
  ## aggregate spectra for all bins using density function on mz.values
  ## weighted by log(1+intensity)
  n.pts <- (10*(max(mz.range) - min(mz.range)))+1
  density.by.bin <- lapply(unique.bins, FUN = function(x) {
    ms.int <- unlist(Spectra::intensity(ms[which(bins == x)]))
    ms.mz  <- unlist(Spectra::mz(ms[which(bins == x)]))
    suppressWarnings(spec.dens <- density(ms.mz, weights = log10(1+ms.int), bw = 0.4, n = n.pts, from = min(mz.range), to = max(mz.range)))
    # plot(spec.dens, main = paste0("BIN: ", x), xlim = c(900,950))
    # abline(v = c(915,939))
    return(spec.dens)
  })
  
  # sig.matrix <- matrix(nrow = length(density.by.bin), ncol = n.pts)
  # dimnames(sig.matrix)[[1]] <- 1:length(density.by.bin)
  # dimnames(sig.matrix)[[2]] <- density.by.bin[[i]]$x
  # for(i in 1:nrow(sig.matrix)) {
  #   sig.matrix[i,] <- density.by.bin[[i]]$y/max(density.by.bin[[i]]$y)
  # }
  # 
  # sig.tall <- reshape2::melt(sig.matrix)
  # names(sig.tall) <- c("bin", "mz", "intensity")
  # dim(sig.tall)
  # library(ggplot2)
  # ggplot(sig.tall, aes(bin, mz, fill= intensity)) + 
  #   geom_tile()
  
  ## for each aggregate spectrum, determine bin boundaries
  bin.boundaries.index <- lapply(unique.bins, FUN = function(x) {
    
    # max.int <- which.max(density.by.bin[[x]]$y)
    ss <- smooth.spline(density.by.bin[[x]]$y, spar = 0.1)$y
    use <- range(which(ss >= 0.33*max(ss)))
    return(use)
    # plot(ss, type = "l")
    # abline(v = use, col = 2)         
    # Sys.sleep(0.5)
  }
  )

  plot(density.by.bin[[x]]$y, type = "l")  
points(max.int, 1, col = 2, cex = 2)
points(smooth.spline(density.by.bin[[x]]$y, spar = 0.2), type = "l", col = 2)
upper.pts <- density.by.bin[[x]]$y[max.int]
upper.der <- diff(density.by.bin[[x]]$y[max.int:length(density.by.bin[[x]]$y)], lag = 1)
plot(upper.der, type = "l")
upper.int.threshold <- which(density.by.bin[[x]]$y <= (0.4* max(density.by.bin[[x]]$y)))
upper.der.threshold <- which(upper.der > (-0.01*max(density.by.bin[[x]]$y)))

lower.infl <- max.int - which(
  c(
    FALSE, 
    diff(
      diff(
        density.by.bin[[x]]$y[max.int:1], lag = 20
      )>(
        0.02*max(
          density.by.bin[[x]]$y
        )
      )
    )
    !=0
  )
)[1]
upper.infl <- max.int + which(
  c(
    FALSE, 
    diff(
      diff(
        density.by.bin[[x]]$y[max.int:length(density.by.bin[[x]]$y)], lag = 20
      )>(
        0.01*max(
          density.by.bin[[x]]$y
        )
      )
    )
    !=0
  )
)[1]
plot(density.by.bin[[x]], type = "l")
abline(v = density.by.bin[[x]]$x[c(lower.infl, upper.infl)], lwd = 0.2, col = 2)
lower.cut <- which(under.cdf > 80)[1]
infl <- which(c(FALSE, diff(diff(under.cdf, lag = 2)>0.5)!=0))[1]
lower.cut <- which(under.cdf > 80)[1]
  
}




## regress/smooth upper and lower bin boundaries to reduce impact of outliers

## fix assignment of upper and lower mz bounds per bin

## return object with structured bin boundary data

