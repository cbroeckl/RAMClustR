library(xcms, quietly=TRUE)
load("inst/exampledata/xset4.Rdata")

source("R/Params.R")
source("R/ramclustR.R")
RC<-ramclustR(xcmsObj=xset4, mspout=TRUE, MStag="01.cdf", idMSMStag="02.cdf")
gc()

##test functionality on xcms object
#load("inst/exampledata/MSdata.Rdata")
#load("inst/exampledata/MSMSdata.Rdata")
#RC<-ramclustR(xcmsObj=xset4, MStag="01.cdf", idMSMStag="02.cdf")