##
##  first set your working directory to directory containing chosen ramclustR directory
##  i.e. C:/R/libraries/

library(devtools)
install(pkg=paste(getwd(), "/RAMClustR", sep=""), dependencies=FALSE)
library(xcms)

library(RAMClustR)
nonxcms_MSonly<- ramclustR (xcmsObj = NULL, ms = paste(getwd(),"/RAMClustR/inst/exampledata/MSdata.csv", sep=""), 
                   idmsms = NULL, 
                   taglocation = "filepaths", 
                   MStag = NULL, idMSMStag = NULL, featdelim = "_", timepos = 2, 
                   st = 5, sr = NULL, maxt = 20, deepSplit = FALSE, 
                   blocksize = 2000, mult = 5, hmax = 0.5, sampNameCol = 1, 
                   collapse = TRUE, mspout = FALSE, mslev = 1, ExpDes = NULL, 
                   normalize = "TIC", minModuleSize = 2, linkage="average")

nonxcms_MSMS<- ramclustR (xcmsObj = NULL, ms = paste(getwd(),"/RAMClustR/inst/exampledata/MSdata.csv", sep=""), 
                             idmsms = paste(getwd(),"/RAMClustR/inst/exampledata/MSMSdata.csv", sep=""), 
                             taglocation = "filepaths", 
                             MStag = NULL, idMSMStag = NULL, featdelim = "_", timepos = 2, 
                             st = 5, sr = NULL, maxt = 20, deepSplit = FALSE, 
                             blocksize = 2000, mult = 5, hmax = 0.5, sampNameCol = 1, 
                             collapse = TRUE, mspout = FALSE, mslev = 2, ExpDes = NULL, 
                             normalize = "TIC", minModuleSize = 2, linkage="average")

## type 'newGC' into the cell next to 'platform',
## on windows, then close window, close window again
load(paste(getwd(),"/RAMClustR/inst/exampledata/xset_GCMS.Rdata", sep=""))
xcms_GCMS<- ramclustR (xcmsObj = xset4,   
                            ms = NULL, idmsms = NULL, 
                            taglocation = "filepaths", 
                            MStag = NULL, idMSMStag = NULL, featdelim = "_", timepos = 2, 
                            st = NULL, sr = NULL, maxt = 20, deepSplit = FALSE, 
                            blocksize = 2000, mult = 5, hmax = 0.5, sampNameCol = 1, 
                            collapse = TRUE, mspout = TRUE, mslev = 1, ExpDes = NULL, 
                            normalize = "TIC", minModuleSize = 2, linkage="average")
rm(xset4)

load(paste(getwd(),"/RAMClustR/inst/exampledata/xset_LC_idMSMS.Rdata", sep=""))
## type 'newLC' into the cell next to 'platform', 
## on windows, then close window, close window again
xcms_LCMS_idMSMS<- ramclustR (xcmsObj = xset, 
                       ms = NULL, idmsms = NULL, 
                       taglocation = "filepaths", 
                       MStag = "01.cdf", idMSMStag = "02.cdf", featdelim = "_", timepos = 2, 
                       st = NULL, sr = NULL, maxt = 20, deepSplit = FALSE, 
                       blocksize = 2000, mult = 5, hmax = 0.5, sampNameCol = 1, 
                       collapse = TRUE, mspout = TRUE, mslev = 1, ExpDes = NULL, 
                       normalize = "TIC", minModuleSize = 2, linkage="average")
rm(xset)

