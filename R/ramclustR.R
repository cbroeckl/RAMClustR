#' compute_SpecAbundAve
#'
#' further aggregate by sample names for 'SpecAbundAve' dataset
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @return  ramclustR object with aggregate by sample names for 'SpecAbundAve' dataset
#' @export
compute_SpecAbundAve <- function(ramclustObj = NULL) {
  if (length(dimnames(ramclustObj$SpecAbund)[[1]]) > length(unique(dimnames(ramclustObj$SpecAbund)[[1]]))) {
    ramclustObj$SpecAbundAve <- aggregate(ramclustObj$SpecAbund[, 1:ncol(ramclustObj$SpecAbund)],
      by = list(dimnames(ramclustObj$SpecAbund)[[1]]),
      FUN = "mean", simplify = TRUE
    )
    dimnames(ramclustObj$SpecAbundAve)[[1]] <- ramclustObj$SpecAbundAve[, 1]
    ramclustObj$SpecAbundAve <- as.matrix(ramclustObj$SpecAbundAve[, 2:ncol(ramclustObj$SpecAbundAve)])
    dimnames(ramclustObj$SpecAbundAve)[[2]] <- dimnames(ramclustObj$SpecAbund)[[2]]
    gc()
  }

  return(ramclustObj)
}

#' ramclustR
#'
#' Main clustering function for grouping features based on their analytical behavior.
#'
#' @param xcmsObj xcmsObject: containing grouped feature data for clustering by ramclustR
#' @param ms filepath: optional csv input. Features as columns, rows as samples. Column header mz_rt
#' @param pheno_csv filepath: optional csv input containing phenoData
#' @param idmsms filepath: optional idMSMS / MSe csv data.  same dim and names as ms required
#' @param MStag character: character string in 'taglocation' to designat MS / MSe files e.g. "01.cdf"
#' @param idMSMStag character: character string in 'taglocation' to designat idMSMS / MSe files e.g. "02.cdf"
#' @param taglocation character: "filepaths" by default, "phenoData[,1]" is another option. refers to xcms slot
#' @param featdelim character: how feature mz and rt are delimited in csv import column header e.g. ="-"
#' @param timepos integer: which position in delimited column header represents the retention time (csv only)
#' @param st numeric: sigma t - time similarity decay value
#' @param sr numeric: sigma r - correlational similarity decay value
#' @param maxt numeric: maximum time difference to calculate retention similarity for - all values beyond this are assigned similarity of zero
#' @param deepSplit logical: controls how agressively the HCA tree is cut - see ?cutreeDynamicTree
#' @param blocksize integer: number of features (scans?) processed in one block  =1000,
#' @param mult numeric: internal value, can be used to influence processing speed/ram usage
#' @param hmax numeric: precut the tree at this height, default 0.3 - see ?cutreeDynamicTree
#' @param sampNameCol integer: which column from the csv file contains sample names?
#' @param collapse logical: reduce feature intensities to spectrum intensities?
#' @param usePheno logical: transfer phenotype data from XCMS object to SpecAbund dataset?
#' @param mspout logical: write msp formatted spectra to file?
#' @param ExpDes either an R object created by R ExpDes object: data used for record keeping and labelling msp spectral output
#' @param normalize character: either "none", "TIC", "quantile", or "batch.qc" normalization of feature intensities.  see batch.qc overview in details.
#' @param qc.inj.range integer: how many injections around each injection are to be scanned for presence of QC samples when using batch.qc normalization?  A good rule of thumb is between 1 and 3 times the typical injection span between QC injections.  i.e. if you inject QC ever 7 samples, set this to between 7 and 21.  smaller values provide more local precision but make normalization sensitive to individual poor outliers (though these are first removed using the boxplot function outlier detection), while wider values provide less local precision in normalization but better stability to individual peak areas.
#' @param batch integer vector with length equal to number of injections in xset or csv file
#' @param order integer vector with length equal to number of injections in xset or csv file
#' @param qc logical vector with length equal to number of injections in xset or csv file.
#' @param minModuleSize integer: how many features must be part of a cluster to be returned? default = 2
#' @param linkage character: heirarchical clustering linkage method - see ?hclust
#' @param mzdec integer: number of decimal places used in printing m/z values
#' @param cor.method character: which correlational method used to calculate 'r' - see ?cor
#' @param rt.only.low.n logical: default = TRUE  At low injection numbers, correlational relationships of peak intensities may be unreliable.  by defualt ramclustR will simply ignore the correlational r value and cluster on retention time alone.  if you wish to use correlation with at n < 5, set this value to FALSE.
#' @param fftempdir valid path: if there are file size limitations on the default ff package temp directory  - getOptions('fftempdir') - you can change the directory used as the fftempdir with this option.
#' @param replace.zeros logical: TRUE by default.  NA, NaN, and Inf values are replaced with zero, and zero values are sometimes returned from peak peaking.  When TRUE, zero values will be replaced with a small amount of noise, with noise level set based on the detected signal intensities for that feature.
#' @details Main clustering function output - see citation for algorithm description or vignette('RAMClustR') for a walk through.  batch.qc. normalization requires input of three vectors (1) batch (2) order (3) qc.   This is a feature centric normalization approach which adjusts signal intensities first by comparing batch median intensity of each feature (one feature at a time) QC signal intensity to full dataset median to correct for systematic batch effects and then secondly to apply a local QC median vs global median sample correction to correct for run order effects.
#' @return   $featclus: integer vector of cluster membership for each feature
#' @return   $frt: feature retention time, in whatever units were fed in (xcms uses seconds, by default)
#' @return   $fmz: feature retention time, reported in number of decimal points selected in ramclustR function
#' @return   $xcmsOrd: the original XCMS (or csv) feature order for cross referencing, if need be
#' @return   $clrt: cluster retention time
#' @return   $clrtsd: retention time standard deviation of all the features that comprise that cluster
#' @return   $nfeat: number of features in the cluster
#' @return   $nsing: number of 'singletons' - that is the number of features which clustered with no other feature
#' @return   $ExpDes: the experimental design object used when running ramclustR.  List of two dataframes.
#' @return   $cmpd: compound name.  C#### are assigned in order of output by dynamicTreeCut.  Compound with the most features is classified as C0001...
#' @return   $ann: annotation.  By default, annotation names are identical to 'cmpd' names.  This slot is a placeholder for when annotations are provided
#' @return   $MSdata:  the MSdataset provided by either xcms or csv input
#' @return   $MSMSdata: the (optional) MSe/idMSMS dataset provided be either xcms or csv input
#' @return   $SpecAbund: the cluster intensities after collapsing features to clusters
#' @return   $SpecAbundAve: the cluster intensities after averaging all samples with identical sample names
#' @return   - 'spectra' directory is created in the working directory.  In this directory a .msp is (optionally) created, which contains the spectra for all compounds in the dataset following clustering.  if MSe/idMSMS data are provided, they are listed width he same compound name as the MS spectrum, with the collision energy provided in the ExpDes object provided to distinguish low from high CE spectra.
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis boxplot hist "legend" "par" "plot" "points" "title"
#' @importFrom stats aggregate cor fitted lm loess median predict quantile sd weighted.mean
#' @importFrom utils edit read.csv read.delim read.delim2 write.csv packageVersion
#' @importFrom ff ff
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom e1071 skewness
#' @importFrom gplots heatmap.2
#' @importFrom pcaMethods pca
#' @importFrom jsonlite fromJSON
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom utils citation packageVersion
#'
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export
#' @examples
#' ## Choose input file with feature column names `mz_rt` (expected by default).
#' ## Column with sample name is expected to be first (by default).
#' ## These can be adjusted with the `featdelim` and `sampNameCol` parameters.
#' wd <- getwd()
#' filename <- system.file("extdata", "peaks.csv", package = "RAMClustR", mustWork = TRUE)
#' print(filename)
#' head(data.frame(read.csv(filename)), c(6L, 5L))
#'
#' ## If the file contains features from MS1, assign those to the `ms` parameter.
#' ## If the file contains features from MS2, assign those to the `idmsms` parameter.
#' ## If you ran `xcms` for the feature detection, the assign the output to the `xcmsObj` parameter.
#' ## In this example we use a MS1 feature table stored in a `csv` file.
#' setwd(tempdir())
#' ramclustobj <- ramclustR(ms = filename, st = 5, maxt = 1, blocksize = 1000)
#'
#' ## Investigate the deconvoluted features in the `spectra` folder in MSP format
#' ## or inspect the `ramclustobj` for feature retention times, annotations etc.
#' print(ramclustobj$ann)
#' print(ramclustobj$nfeat)
#' print(ramclustobj$SpecAbund[, 1:6])
#' setwd(wd)
#'
ramclustR <- function(xcmsObj = NULL,
                      ms = NULL,
                      pheno_csv = NULL,
                      idmsms = NULL,
                      taglocation = "filepaths",
                      MStag = NULL,
                      idMSMStag = NULL,
                      featdelim = "_",
                      timepos = 2,
                      st = NULL,
                      sr = NULL,
                      maxt = NULL,
                      deepSplit = FALSE,
                      blocksize = 2000,
                      mult = 5,
                      hmax = NULL,
                      sampNameCol = 1,
                      collapse = TRUE,
                      usePheno = TRUE,
                      mspout = TRUE,
                      ExpDes = NULL,
                      normalize = "TIC",
                      qc.inj.range = 20,
                      order = NULL,
                      batch = NULL,
                      qc = NULL,
                      minModuleSize = 2,
                      linkage = "average",
                      mzdec = 3,
                      cor.method = "pearson",
                      rt.only.low.n = TRUE,
                      fftempdir = NULL,
                      replace.zeros = TRUE) {
  ########
  # If experimental design is NULL:
  if (is.null(ExpDes)) {
    ExpDes <- defineExperiment(force.skip = TRUE)
    warning(
      "\n", "you failed to define your experimental descriptor using 'defineExperiment()'", "\n",
      "RAMClustR must now guess at what you are trying to do ", "\n",
      "and your exported spectra will be labelled incorrectly"
    )
    if (!is.null(idmsms)) {
      ExpDes[[2]][which(row.names(ExpDes[[2]]) == "MSlevs"), 1] <- 2
    }
  }

  if (normalize == "batch.qc") {
    if (is.null(order) | is.null(batch) | is.null(qc)) {
      stop("to use batch.qc normalization you must provide vectors for batch, order (run order) and qc information as vectors.  see help ?ramclustR")
    }
  }

  if (!is.null(fftempdir)) {
    origffdir <- getOption("fftempdir")
    options("fftempdir" = fftempdir)
  }

  ########
  # define ms levels, used several times below
  mslev <- as.integer(as.numeric(as.character(ExpDes[[2]][which(row.names(ExpDes[[2]]) == "MSlevs"), 1])))

  history <- paste("Raw mass spectrometry data were processed using an R based workflow for feature detection, retention time alignment, feature grouping, peak filling, feature clustering.")

  ########
  # do some checks to make sure we have everything we need before proceeding
  if (is.null(xcmsObj) & is.null(ms)) {
    stop("you must select either
          1: an MS dataset with features as columns
             (__)one column may contain sample names, defined by sampNameCol)
          2: an xcmsObj. If you choose an xcms object, set taglocation: 'filepaths' by default
            and both MStag and idMSMStag")
  }

  if (!is.null(xcmsObj) & mslev == 2 & any(is.null(MStag), is.null(idMSMStag), is.null(taglocation))) {
    stop("you must specify the the MStag, idMSMStag, and the taglocations")
  }

  if (normalize != "none" & normalize != "TIC" & normalize != "quantile" & normalize != "batch.qc") {
    stop("please selected either 'none', 'TIC', or 'quantile' for the normalize setting")
  }


  a <- Sys.time()

  ########
  # set default parameters of not defined in function call

  if (is.null(hmax)) {
    hmax <- 0.3
  }

  cat(paste("organizing dataset", "\n"))

  # read csv through rc.get.csv.data
  if (!is.null(ms)) {
    ramclustObj <- rc.get.csv.data(
      csv = ms,
      phenoData = pheno_csv,
      idmsms = idmsms,
      ExpDes = ExpDes,
      sampNameCol = sampNameCol,
      st = st,
      timepos = timepos,
      featdelim = featdelim,
      ensure.no.na = replace.zeros
    )
  }

  ########
  # if xcms object as input, do this:
  if (!is.null(xcmsObj)) {
    ramclustObj <- rc.get.xcms.data(
      xcmsObj = xcmsObj,
      taglocation = taglocation,
      MStag = MStag,
      MSMStag = idMSMStag,
      ExpDes = ExpDes,
      mzdec = mzdec,
      ensure.no.na = replace.zeros
    )

    history <- paste(history,
      " XCMS", paste0("(v.", packageVersion("xcms"), ")"), "was used for feature detection and retention time alighment. ",
      "Processing was performed using R", paste0("(v.", citation()$author, " ", citation()$year, ")."),
      " Feature data was input as an xcms object with ramclustR parameter settings of ",
      " st = ", st,
      " sr = ", sr,
      " and maxt = ", maxt, ".",
      sep = ""
    )
  }

  history <- paste(history, "RAMClustR (version ",
    packageVersion("RAMClustR"),
    ") was utilized to cluster features into spectra (Broeckling 2014).",
    sep = ""
  )
  
  if (mslev == 2) {
   history <- paste(
     history,
     "Feature data included both MS and indiscriminant MS/MS data (Broeckling 2012)."
   )
  }

  ########
  # ensure that we have all numeric values,
  # then optionally ensure we have all non-zero values in the dataset.
  # uses a noise addition 'jitter' around minimum values with missing data points.
  ramclustObj <- rc.feature.replace.na(
    ramclustObj = ramclustObj,
    replace.zero = replace.zeros
  )

  ########
  # Optional normalization of data, either Total ion signal or quantile

  if (normalize != "none") {
    cat(paste("  normalizing dataset", "\n"))
  }

  if (normalize == "TIC") {
    ramclustObj <- rc.feature.normalize.tic(ramclustObj = ramclustObj)

    history <- paste(
      history,
      " Features were normalized to total signal using 'tic' normalization."
    )
  }

  if (normalize == "quantile") {
    ramclustObj <- rc.feature.normalize.quantile(ramclustObj)

    history <- paste(
      history,
      " Features were normalized using 'quantile' normalization."
    )
  }

  if (normalize == "batch.qc") {
    ramclustObj <- rc.feature.normalize.batch.qc(
      ramclustObj = ramclustObj,
      qc.inj.range = qc.inj.range
    )

    history <- paste(history,
     " Features were normalized to nearby QC samples on a feature-by-feature basis using the 'batch.qc' option "
    )
  }

  ramclustObj <- rc.ramclustr(
    ramclustObj = ramclustObj,
    st = st,
    sr = sr,
    maxt = maxt,
    deepSplit = deepSplit,
    blocksize = blocksize,
    mult = mult,
    hmax = hmax,
    collapse = collapse,
    minModuleSize = minModuleSize,
    linkage = linkage,
    cor.method = cor.method,
    rt.only.low.n = rt.only.low.n,
    fftempdir = fftempdir
  )

  history <- paste(
    history,
    "The feature similarity matrix was clustered using fastcluster package heirarchical clustering method using the",
    linkage, "method."
  )

  if (!is.null(ms) & normalize == "quantile") {
    dimnames(ramclustObj$SpecAbund)[[1]] <- dimnames(ramclustObj$MSdata)[[1]]
  }

  ########
  # further aggregate by sample names for 'SpecAbundAve' dataset
  if (!is.null(ramclustObj$SpecAbund)) {
    ramclustObj <- compute_SpecAbundAve(ramclustObj)
  }
  gc()

  ########
  # write msp formatted spectra
  if (mspout == TRUE) {
    cat(paste("writing msp formatted spectra", "\n"))
    write.msp(ramclustObj, one.file = TRUE)
    cat(paste("msp file complete", "\n"))
  }

  ramclustObj$history <- history
  if (nrow(ramclustObj$MSdata) < 5 & rt.only.low.n) {
    warning("\n", "too few samples to use correlational similarity, clustering by retention time only", "\n")
    ramclustObj$history <- paste(
      ramclustObj$history,
      "Since there were fewer than five injections, clustering was performed only using retention time simiilarity."
    )
  }
  return(ramclustObj)
}
