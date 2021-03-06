% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc.feature.normalize.tic.R
\name{rc.feature.normalize.tic}
\alias{rc.feature.normalize.tic}
\title{rc.feature.normalize.tic}
\usage{
rc.feature.normalize.tic(ramclustObj = NULL)
}
\arguments{
\item{ramclustObj}{ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)}
}
\value{
ramclustR object with total extracted ion normalized data.
}
\description{
extractor for xcms objects in preparation for clustering
}
\details{
This function offers normalization by total extracted ion signal.  it is recommended to first run 'rc.feature.filter.blanks' to remove non-sample derived signal.
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
