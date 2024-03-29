% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc.calibrate.ri.R
\name{rc.calibrate.ri}
\alias{rc.calibrate.ri}
\title{rc.calibrate.ri}
\usage{
rc.calibrate.ri(ramclustObj = NULL, calibrant.data = "", poly.order = 3)
}
\arguments{
\item{ramclustObj}{ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)}

\item{calibrant.data}{character vector defining the file path/name to a csv file containing columns including 'rt', and 'ri'.  Alternatively, a data.frame with those columnn names (case sensitive)}

\item{poly.order}{integer default = 3.  polynomical order used to fit rt vs ri data, and calculate ri for all feature and metabolite rt values. poly.order should be apprciably smaller than the number of calibrant points.}
}
\value{
ramclustR object with retention index assigned for features ($fri) and compounds ($clri).
}
\description{
extractor for xcms objects in preparation for clustering
}
\details{
This function generates a new slot in the ramclustR object for retention index.  Calibration is performed using a polynomial fit of order poly.order.  It is the user's responsibility to ensure that the number and span of calibrant points is sufficient to calibrate the full range of feature and compound retention times.  i.e. if the last calibration point is at 1000 seconds, but the last eluting peak is at 1300 seconds, the calibration will be very poor for the late eluting compound.
}
\references{
Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
}
\author{
Corey Broeckling
}
\concept{MSFinder}
\concept{RAMClustR}
\concept{clustering}
\concept{feature}
\concept{mass spectrometry}
\concept{metabolomics}
\concept{ramclustR}
\concept{xcms}
