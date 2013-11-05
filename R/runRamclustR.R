#small dataset for devel
load("inst/exampledata/MSdata.Rdata")
load("inst/exampledata/MSMSdata.Rdata")

#bigger datasets for testing
#load("CSFMSdata.Rdata")
#load("CSFMSMSdata.Rdata")

source("R/ramclustR.R")
#source("UPLC_C18params.R")
#source("mspout.R")

RC<-ramclustR(ms=MSdata, idmsms=MSMSdata, sr=0.2, st=2, hmax=1.05, mspout=TRUE)
gc()