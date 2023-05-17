# RAMClustR: Mass Spectrometry Metabolomics Feature Clustering and Interpretation
A feature clustering algorithm for non-targeted mass spectrometric metabolomics data. This method is compatible with gas and liquid chromatography coupled mass spectrometry, including indiscriminant tandem mass spectrometry data.

## Documentation for users

### Installation
The newest version of the package can be installed through conda from the [bioconda](https://anaconda.org/bioconda/r-ramclustr) channel:

```bash
conda install -c bioconda r-ramclustr
```

Or you can alternatively Install from R console:

install.packages("devtools", repos="http://cran.us.r-project.org", dependencies=TRUE)

library(devtools)

install_github("cbroeckl/RAMClustR", build_vignettes = TRUE, dependencies = TRUE)

library(RAMClustR)

vignette("RAMClustR")

### Introduction
Main clustering function output - see citation for algorithm description or vignette('RAMClustR') for a walk through. batch.qc. normalization requires input of three vectors (1) batch (2) order (3) qc. This is a feature centric normalization approach which adjusts signal intensities first by comparing batch median intensity of each feature (one feature at a time) QC signal intensity to full dataset median to correct for systematic batch effects and then secondly to apply a local QC median vs global median sample correction to correct for run order effects.

There are two pathways for using RAMClustR; You can use either use the main ramclustR function or the individual stepwise workflow. 

Below is a small example of using main ramclustR function.
```R
## Choose input file with feature column names `mz_rt` (expected by default).
## Column with sample name is expected to be first (by default).
## These can be adjusted with the `featdelim` and `sampNameCol` parameters.
wd <- getwd()
filename <- file.path(wd, "testdata/peaks.csv")
pheno <- file.path(wd, "testdata/phenoData.csv") 
print(filename)
head(data.frame(read.csv(filename)), c(6L, 5L))

## If the file contains features from MS1, assign those to the `ms` parameter.
## If the file contains features from MS2, assign those to the `idmsms` parameter.
## If you ran `xcms` for the feature detection, the assign the output to the `xcmsObj` parameter.
## In this example we use a MS1 feature table stored in a `csv` file.
setwd(tempdir())
ramclustobj <- ramclustR(
    ms = filename,
    pheno_csv = pheno,
    st = 5,
    maxt = 1,
    blocksize = 1000
  )

## Investigate the deconvoluted features in the `spectra` folder in MSP format
## or inspect the `ramclustobj` for feature retention times, annotations etc.
print(ramclustobj$ann)
print(ramclustobj$nfeat)
print(ramclustobj$SpecAbund[,1:6])
setwd(wd)
```

#### Individual stepwise workflow
![alt text](https://github.com/zargham-ahmad/RAMClustR/blob/issue_14/docs/ramclustR.png)

Below is a small example of using Individual stepwise workflow.
```R
set.seed(123) # to get reproducible results with jitters
wd <- getwd()
tmp <- tempdir()
load(file.path("testdata", "test.rc.ramclustr.fillpeaks"))

setwd(tmp)

ramclustObj <- rc.get.xcms.data(xcmsObj = xdata)
ramclustObj <- rc.expand.sample.names(ramclustObj = ramclustObj)
ramclustObj <- rc.feature.replace.na(ramclustObj = ramclustObj)
ramclustObj <- rc.feature.filter.blanks(ramclustObj = ramclustObj, blank.tag = "Blanc")
ramclustObj <- rc.feature.normalize.qc(ramclustObj = ramclustObj, qc.tag = "QC")
ramclustObj <- rc.feature.filter.cv(ramclustObj = ramclustObj)
ramclustObj <- rc.ramclustr(ramclustObj = ramclustObj)
ramclustObj <- rc.qc(ramclustObj = ramclustObj)
ramclustObj <- do.findmain(ramclustObj = ramclustObj)

## Investigate the deconvoluted features in the `spectra` folder in MSP format
## or inspect the `ramclustobj` for feature retention times, annotations etc.
print(ramclustobj$ann)
print(ramclustobj$nfeat)
print(ramclustobj$SpecAbund[,1:6])
setwd(wd)
```

## Documentation for developers

### Installation
```bash
git clone https://github.com/cbroeckl/RAMClustR.git
cd RAMClustR
conda env create -n ramclustr-dev -f=conda/environment-dev.yaml
conda activate ramclustr-dev
```

### Testing
```R
# Activate the ramclustr-dev environment
# Run the below command on R console
devtools::test()
```

## References
Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d. Epub 2014 Jun 26. PubMed PMID: 24927477.

Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.

