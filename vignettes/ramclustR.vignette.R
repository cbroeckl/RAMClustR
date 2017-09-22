## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----xcms faaKO, eval=FALSE, include=TRUE--------------------------------
#  library(xcms)
#  library(faahKO)
#  cdfpath <- system.file("cdf", package = "faahKO")
#  cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#  xset <- xcmsSet(cdffiles)  # detect features
#  xset <- group(xset)  # group features across samples by retention time and mass
#  xset <- retcor(xset, family = "symmetric", plottype = NULL)  # correct for drive in retention time
#  xset <- group(xset, bw = 10)  # regroup following rt correction
#  xset <- fillPeaks(xset)  # 'fillPeaks' to remove missing values in final dataset

## ----view xcms object summary, eval=FALSE, include=TRUE------------------
#  xset

## ----ramclustR installation, eval=FALSE, include=TRUE--------------------
#  install.packages("devtools", repos="http://cran.us.r-project.org", dependencies=TRUE)
#  library(devtools)
#  install_github("cbroeckl/RAMClustR")
#  library(RAMClustR)

## ----ramclustR of xcms processed faaKO, eval=FALSE, include=TRUE---------
#  experiment <- defineExperiment(csv = TRUE)
#  RC <- ramclustR(xcmsObj = xset, ExpDes=experiment)

## ----export csv, eval=FALSE, include=TRUE--------------------------------
#  write.csv(RC$SpecAbund, file="SpecAbund.csv", row.names=TRUE)

## ----csv input, eval=FALSE, include=TRUE---------------------------------
#  experiment <- defineExperiment(csv = TRUE)
#  RC <- ramclustR(ms = "mymsdata.csv", featdelim - "_", timepos = 2, st = 5, ExpDes=experiment)

## ----cars----------------------------------------------------------------
summary(cars)

## ----pressure, echo=FALSE------------------------------------------------
plot(pressure)

