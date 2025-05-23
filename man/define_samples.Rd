% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc.feature.filter.cv.R
\name{define_samples}
\alias{define_samples}
\title{define_samples}
\usage{
define_samples(ramclustObj, tag, return.logical = FALSE)
}
\arguments{
\item{ramclustObj}{ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)}

\item{tag}{character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.}

\item{return.logical}{optionally convert numeric vector with length equal to the number of matched samples to a logical vector of length equal to number of samples, with TRUE representing matching samples.}
}
\value{
samples found using the tag
}
\description{
define samples in each set
}
