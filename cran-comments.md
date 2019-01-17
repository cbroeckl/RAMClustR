---
title: "cran-comments"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## This document enables CRAN checks and describes results for record
```{cran check}
devtools::check('C:/Users/cbroeckl/Documents/GitHub/RAMClustR')
devtools::build_win('C:/Users/cbroeckl/Documents/GitHub/RAMClustR', version ="R-release")

```

## pasted output from console
Updating RAMClustR documentation
Loading RAMClustR
Writing NAMESPACE
Writing NAMESPACE
Setting env vars -------------------------------------------------------------------------------------
CFLAGS  : -Wall -pedantic
CXXFLAGS: -Wall -pedantic
Building RAMClustR -----------------------------------------------------------------------------------
"C:/PROGRA~1/R/R-35~1.1/bin/x64/R" --no-site-file --no-environ --no-save --no-restore --quiet CMD  \
  build "C:\Users\cbroeckl\Documents\GitHub\RAMClustR" --no-resave-data --no-manual 

* checking for file 'C:\Users\cbroeckl\Documents\GitHub\RAMClustR/DESCRIPTION' ... OK
* preparing 'RAMClustR':
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
* building 'RAMClustR_0.99.001.tar.gz'

Setting env vars -------------------------------------------------------------------------------------
_R_CHECK_CRAN_INCOMING_USE_ASPELL_: TRUE
_R_CHECK_CRAN_INCOMING_           : FALSE
_R_CHECK_FORCE_SUGGESTS_          : FALSE
Checking RAMClustR -----------------------------------------------------------------------------------
"C:/PROGRA~1/R/R-35~1.1/bin/x64/R" --no-site-file --no-environ --no-save --no-restore --quiet CMD  \
  check "C:\Users\cbroeckl\AppData\Local\Temp\Rtmp00ZA6R/RAMClustR_0.99.001.tar.gz" --as-cran  \
  --timings --no-manual 

* using log directory 'C:/Users/cbroeckl/AppData/Local/Temp/Rtmp00ZA6R/RAMClustR.Rcheck'
* using R version 3.5.1 (2018-07-02)
* using platform: x86_64-w64-mingw32 (64-bit)
* using session charset: ISO8859-1
* using options '--no-manual --as-cran'
* checking for file 'RAMClustR/DESCRIPTION' ... OK
* checking extension type ... Package
* this is package 'RAMClustR' version '0.99.001'
* package encoding: UTF-8
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking serialization versions ... OK
* checking whether package 'RAMClustR' can be installed ... OK
* checking installed package size ... NOTE
  installed size is 48.6Mb
  sub-directories of 1Mb or more:
    exampledata  48.3Mb
* checking package directory ... OK
* checking 'build' directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking examples ... NONE
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in 'inst/doc' ... OK
* checking re-building of vignette outputs ... OK
* DONE

Status: 1 NOTE
See
  'C:/Users/cbroeckl/AppData/Local/Temp/Rtmp00ZA6R/RAMClustR.Rcheck/00check.log'
for details.


R CMD check results
0 errors | 0 warnings | 1 note 
checking installed package size ... NOTE
  installed size is 48.6Mb
  sub-directories of 1Mb or more:
    exampledata  48.3Mb

## Notes on results: 
One note reflecting the presence of xcms objects that were in the past used as test data.  Consider deleting them to reduce package size.  
