\name{ramclustR}
\alias{ramclustR}
\title{ramclustR}
\usage{
  ramclustR (xcmsObj = NULL, ms = NULL, idmsms = NULL, taglocation = "filepaths", 
    MStag = NULL, idMSMStag = NULL, featdelim = "_", mzpos = 1, timepos = 2, 
    st = NULL, sr = NULL, maxt = NULL, deepSplit = FALSE, 
    blocksize = 2000, mult = 5, hmax = NULL, sampNameCol = NULL, 
    collapse = TRUE, mspout = TRUE, mslev = 2, ExpDes = NULL, 
    normalize = "TIC", minModuleSize = 2, linkage="average")
}
\arguments{
  \item{xcmsObj}{the xcms object to be used as input for clustering}

  \item{mslev}{the number of full scan collision energies used for the data in the xcms object. If only low collision energy is used, set this parameter equal to 1. If indiscriminant MS/MS is used, set to 2}

  \item{taglocation}{The location for the tag in the xcms object in which to look for the idMSMStag (see below) default is "filepaths", either "filepaths" or "phenoData".}

  \item{MStag}{character e.g. "01.cdf".  THe character string in the taglocation by which to tag files as indiscriminant MSMS data}

  \item{idMSMStag}{character e.g. "02.cdf".  The character string in the taglocation by which to tag files as indiscriminant MSMS data}

  \item{ms}{alternate input format for low collision energy data.  If peak detection is performed by other software, csv formatted data is requred, and the file name (with full path, if not in working directory) is specified here}

  \item{idmsms}{alternate input format for high collision energy data (MS^e, idMS/MS).  If peak detection is perfomred by other software, csv formatted data is requred, and the file name (with full path, if not in working directory) is specified here. feature and sample names must be identical to those in the ms data format.  If using only high collision energy data (or GC-MS, electron impact), use the ms data option instead.}

  \item{featdelim}{alternate input format option: if ms is specified instead of xcmsObj, feature names must contain mass and retention time.  The featdelim parameter sets the delimiter used to separate mass and time components of the feature name.  character e.g. "_".}

  \item{sampNameCol}{The numeric (integer) position of the column in the ms and or idmsms datasets that contains the sample names.}

  \item{timepos}{alternate input format option: If importing data from another peak detection tool, feature names must contain mass and retention time. The position of the time component of the feature name (numeric, either 1 or 2) is defined here.}

  \item{mzpos}{alternate input format option: If importing data from another peak detection tool, feature names must contain mass and retention time. The position of the mz component of the feature name (numeric, either 1 or 2) is defined here.}

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

  \item{ExpDes}{The experimental design information for the experiment.  Will be filled by the user using the edit function.  The editor used by this function will vary based on the platform.  On a windows platform, a popup spreadsheet-style editor will appear, allowing users to enter the descritions for the experimental parameters.  If you wish to save these parameters, enter a new parameter set name in the 'saveAs' cell of the second sheet.  The values input will be used during msp spectrum export.} 
}
\value{
an S3 RAMClustR object. THe base object structure is that of an hclust object, thus the plot methods available for other R objects of hlcust format are available.  In addition, several slots are added:
  featclus: the cluster to which each feature belongs. zero values indicate singletons when minModuleSize is set to a value > 1. 
  frt: feature retention time used in clustering
  fmz: feature m/z value used to draw spectrum
  nfeat: the number of features contained in each cluster
  nsing: the number of singletons following clustering (unclustered features, only when minModuleSize is set to a value > 1. 
  SpecAbund: the reduced dataset.  Features within a cluster are aggregated using a weighted mean function into a single value for each file and cluster.  The weight used is the intensity of the peak, thereby allowing more intense signals to contribute more to the cluster value, reducing the influence of noise in the quantitative value.  
  SpecAbundAve: an (optionally) further reduced dataset in which the injection replicates have been averaged. Injection replicates are taken from the xcms 'phenoData' slot, and samples with identical names (and taken from the same collision energy raw data) will be assumed to be replicate injections.  If ms is set instead of xcmsObj, the sample names will be taken from the optional sampNameCol argument.
  ExpDes: the instrument parameters used for generation of raw data, the values of which are (optionally) written to the output msp file. 
  MSdata: a matrix version of the low collision energy feature intensities
  MSMSdata: a matrix version of the high collision energy intensities
  cmpd: a cluster name: by default, all clusters are assigned the label C1... Cn, where n is the number of clusters with at least minModuleSize features
  ann: a place holder for compound annotations.  by default, all clusters are annotated with the label C1... Cn, where n is the number of clusters with at least minModuleSize features
  annconf: a placeholder for a numerical annotation confidence for each annotation (see Metabolomics Standards initiative for more information)
  annnotes: a placeholder for a vector of annotation notes: in the event of an ambiguous annotation, notes can be used to document the evidence and/or caveats associated with the annotation.
Spectra groups are then used for two purposes.  1. to guide the condenstation of the low collision energy dataset from features to spectral abundances, whereby the weighted mean of each feature in the group is taken for each sample, such that a spectral intensity is calculated for each compound (cluster, spectra) in each sample. 2. to redraw spectra.  Spectra can be exported in .msp format, allowing for further spectral searching using NIST MSsearch or PepSearch, or to be uploaded to the MassBank batch search tool.  Upon annotation of compounds, users are encouraged to store their annotations, Metabolomics Standards Initiative confidence scores, and annotation notes or observations in the ramclust object, slots $ann, $anncof, and $annnotes, respectively. 
}

\description{
  RAMClustR: a feature clustering algorithm designed to groupd features into compounds in an unsupervised manner. 
}

\details{
  The RAMClustR package implements a feature grouping method which determines retention time and correlational similarity of all features in a metabolomics dataset.  These scores are induced to decay using independently tuned sigmoid functions (st and sr for time and correlation, respectively).  A total similarity score is then calculated as the product of the time and correlational similarities.  This numerical similarity score is calculated for all features in the dataset within maxt of each other (20 seconds by default), thereby generating a full similarity matrix reflecting high similarity scores between features only when the two features (approximately) coelute and are correlated with each other across the dataset.  for datasets with many features, memory is managed by processed in blocksize^2 'chunks', reducing the risk of memory limitation induced failures.  
This similarity matrix is used as input to heirarchical clustering via the fastcluster package, using (by default) average linkage.  The output of the clustering is then cut using the dynamicTreeCut package, which prunes complex trees, rather than by branch height, but the relative 'shape' of the clusters, as determined by local branch height patterns.  The cutting process generates discrete spectral grouping assignments for each feature of the dataset.  A RAMClustR object is returned.
There are a few practical considerations to keep in mind while processing your data.  
1. In-source spectra will contain signals above the expected molecular ions, and many compounds generate significan in-source fragments.  Do not be alarmed.  This is observation is what drove us to build RAMClustR.  If you are quite sure that the clusters contain signals from multiple compounds, first reduce sr (from the default value of 0.5), then change deepSpilt to TRUE, and if all else fails, reduce the st value (by default, have the median peak width (when taken from xcms, the median of the peak stop time minus peak start time of all detected peaks).    
2. Since we are using in-source spectra, RAMClustR spectra will often contain masses above the precursor ions typically targeted by generators of MS/MS libraries.  This lack of a perfect match reflects differences in acquisition methods, more than evidence of a poor match. 
3. The RAMClustR spectra can only be as good as the peak detection.  Ensure that you have set XCMS (or any other peak picking software) appropriately.  We generally recommend setting signal to noise levels relatively low (3 or 5) so as to enable detection (and hence clustering) of as many potentially diagnostic signals as possible for more accurate and rich spectra.  In addition, consider that if you have performed peak detection on MS and idMS/MS raw data, and a given feature is present in all the samples, but only at the MS level, the resulting XCMS object will only observe that feature in 50% of the files.  To recover more features in the final dataset, set minfrac rather low, i.e. 10%.  
4. Since good quantitative peak detection is dependent on good separation, minimally by one of mass or retention time, features with a common mass and retention time that are derived from two compounds will naturally not represent the quantitative abundance of either of the two compounds well.  Features such as this may either cluster with the contributing compound that has higher intensity, or may not cluster with either of the two compounds, thereby generating a singleton feature.  A good example of this is the phosphotidyl cholines by reversed phase LC-MS, which generate a common 184 ion, or monosaccharides by GC-MS which generate common 73 and 147 ions.  
}
\author{
  Corey Broeckling
}

