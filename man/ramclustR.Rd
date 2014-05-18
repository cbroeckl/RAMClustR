\name{ramclustR}
\alias{RAMClustR}
\title{ramclustR}
\usage{
  ramclustR (xcmsObj = NULL, ms = NULL, idmsms = NULL, taglocation = "filepaths", 
    MStag = NULL, idMSMStag = NULL, featdelim = "_", timepos = 2, 
    mzpos = 1, st = NULL, sr = NULL, maxt = NULL, deepSplit = FALSE, 
    blocksize = 2000, mult = 5, hmax = NULL, sampNameCol = NULL, 
    collapse = TRUE, mspout = TRUE, mslev = 2, ExpDes = NULL, 
    normalize = "TIC", minModuleSize = 2, linkage="average")
}
\arguments{
  \item{xcmsObj}{the xcms object to be used as input for clustering}

  \item{mslev}{the number of full scan collision energies used for the data in the xcms object. If only low collision energy is used, set this parameter equal to 1. If indiscriminant MS/MS is used, set to 2}

  \item{normalize}{Method used to normalize the xcms peak intensities.  Method options currently include 'TIC' (default), 'quantile', or 'none'}

  \item{taglocation}{The location for the tag in the xcms object in which to look for the idMSMStag (see below) default is "filepaths", either "filepaths" or "phenoData".}

  \item{MStag}{character e.g. "01.cdf".  THe character string in the taglocation by which to tag files as indiscriminant MSMS data}

  \item{idMSMStag}{character e.g. "02.cdf".  The character string in the taglocation by which to tag files as indiscriminant MSMS data}

  \item{ms}{alternate input format for low collision energy data.  If peak detection is performed by other software, csv formatted data is requred, and the file name (with full path, if not in working directory) is specified here}

  \item{idmsms}{alternate input format for high collision energy data (MS^e, idMS/MS).  If peak detection is perfomred by other software, csv formatted data is requred, and the file name (with full path, if not in working directory) is specified here. feature and sample names must be identical to those in the ms data format.  If using only high collision energy data (or GC-MS, electron impact), use the ms data option instead.}

  \item{featdelim}{alternate input format option: if ms is specified instead of xcmsObj, feature names must contain mass and retention time.  The featdelim parameter sets the delimiter used to separate mass and time components of the feature name.  character e.g. "_".}

  \item{sampNameCol}{The numeric (integer) position of the column in the ms and or idmsms datasets that contains the sample names.}

  \item{timepos}{alternate input format option: If importing data from another peak detection tool, feature names must contain mass and retention time. The position of the time component of the feature name (numeric, either 1 or 2) is defined here.}

  \item{st}{Sigma t value controlling the rate at which feature similarity decays with respect to differences in retention times between feature pairs. When an xcms object is used as input, this parameter is by default set to half the peak width of all the peaks in the xcmsSet.  Can be manually overridden to tune performance when necessary.  Larger values are less selective, smaller values are more selective.}

  \item{sr}{Sigma r value controlling the rate at which feature similarity decays with respect to correlational simiarity between a pair of features.  numeric e.g. 0.5, Sigma r: the rate of similarity decay in response to correlational r .  Can be manually overridden to tune performance.  Larger values are less selective, smaller values are more selective.}

  \item{hmax}{a value passed to dynamicTreeCut.  The maximum allowable branch height.  The dendrogram is first cut at branches containing longer than this value, before the shape-based cutting occurs.  Tends to result in more reliable cutting when used, default = 0.3.}

  \item{maxt}{numeric e.g 20, the maximum time allowed for similarity between features to be greater than zero.  For example, if two features are separated by more than 20 seconds, the similarity score is automatically set to zero, ensuring that they end up in different clusters.  If you are importing data from another vendor, ensure that either you change the feature times to seconds (default in ramclustR) or you change this value to a value representing minutes, rather than seconds (e.g. 0.33).  Used primarily to speed up compuation time by not calculating values above this threshold.}

  \item{deepSplit}{boolean e.g. =FALSE,  access to the deepSplit function in dynamicTreeCut.  False generates fewer spectra, TRUE splits spectra more readily and thereby generates more spectra. For more detailed description, see documentation for the dynamicTreeCut pacakage.}

  \item{blocksize}{integer number of features processed in one block =2000, ramclustR makes attempts to minimize memory demands to maintain processing speed by avoiding memory constraints.  It does so by processing the full similarity matrix in blocks of blocksize dimensions.  This has two benefits: one is that a given PC is less likely to run out of memory during processing and the second is that we need not process blocks that contain not feature pairs with a delta rt of greater than maxt, thereby avoiding needless calculations and processing time.}

  \item{minModuleSize}{the number of features required for a cluster to be reported. That is, if a feature is clustered in a group with less than minModuleSize features, it will not be exported as a SpecAbund or SpecAbundAve column. If you wish to analyze singletons, set minModuleSize=1}

  \item{mult}{constant for reshaping blocksize during custom similarity matrix generation, generally unecessary to set manually.  If memory demands become apparent later in processing, reduce the value.}

  \item{collapse}{logical, default=TRUE; controls whether RAMClustR should collapse feature intensities into spectral intensities after cluster membership for all features has been determined.}

  \item{mspout}{logical, default=TRUE; controls whether RAMClustR writes an msp formatted spectral file containing spectra for all compounds with more than minModuleSize features}

  \item{normalize}{One of "TIC", "quantile", or "none" currently.  Controls how the xcms data is normilized before clustering into spectra.  Normalization is important when MS sensitivity drift over the experiment might artificially induce correlation between features derived from two compounds that coelute.}

  \item{linkage}{The linkage method used to perform heirarchical fastcluster-based heirarchical clustering.  default="average"}

  \item{ExpDes}{The experimental design information for the experiment.  Will be filled by the user using the edit function.  THe editor used by this function will vary based on the platform.  On a windows platform, a popup spreadsheet-style editor will appear, allowing users to enter the descritions for the experimental parameters.  If you wish to save these parameters, enter a new parameter set name in the 'saveAs' cell of the second sheet.  The values input will be used during msp spectrum export.} 

\author{
  Corey Broeckling
}

