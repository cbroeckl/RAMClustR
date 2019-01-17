---
title: "cran-comments"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## This document enables CRAN checks and describes results for record
```{cran check}
devtools::check_win_release('C:/Users/cbroeckl/Documents/GitHub/RAMClustR')

```

## pasted output
* using log directory 'd:/RCompile/CRANguest/R-release/RAMClustR.Rcheck'
* using R version 3.5.2 (2018-12-20)
* using platform: x86_64-w64-mingw32 (64-bit)
* using session charset: ISO8859-1
* checking for file 'RAMClustR/DESCRIPTION' ... OK
* checking extension type ... Package
* this is package 'RAMClustR' version '1.00.00'
* package encoding: UTF-8
* checking CRAN incoming feasibility ... NOTE
Maintainer: '"Broeckling,Corey" <Corey.Broeckling@ColoState.EDU>'

New submission

Version contains leading zeroes (1.00.00)

Possibly mis-spelled words in DESCRIPTION:
  GC (13:131)
  MSall (13:193)
  MSe (13:188)
  Metabolomics (4:48)
  RAMClust (3:8)
  indiscriminant (13:166)
  metabolomics (13:81)
  spectrometric (13:67)

The Title field should be in title case, current version then in title case:
'RAMClust: A Novel Feature Clustering Method Enables Spectral-Matching Based Annotation For Metabolomics Data'
'RAMClust: A Novel Feature Clustering Method Enables Spectral-Matching Based Annotation for Metabolomics Data'
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking serialization versions ... OK
* checking whether package 'RAMClustR' can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking 'build' directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* loading checks for arch 'i386'
** checking whether the package can be loaded ... OK
** checking whether the package can be loaded with stated dependencies ... OK
** checking whether the package can be unloaded cleanly ... OK
** checking whether the namespace can be loaded with stated dependencies ... OK
** checking whether the namespace can be unloaded cleanly ... OK
** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* loading checks for arch 'x64'
** checking whether the package can be loaded ... OK
** checking whether the package can be loaded with stated dependencies ... OK
** checking whether the package can be unloaded cleanly ... OK
** checking whether the namespace can be loaded with stated dependencies ... OK
** checking whether the namespace can be unloaded cleanly ... OK
** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... [26s] OK
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
* checking re-building of vignette outputs ... [2s] OK
* checking PDF version of manual ... OK
* DONE
Status: 1 NOTE
