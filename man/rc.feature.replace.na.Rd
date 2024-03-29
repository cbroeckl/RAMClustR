% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc.feature.replace.na.R
\name{rc.feature.replace.na}
\alias{rc.feature.replace.na}
\title{rc.feature.replace.na}
\usage{
rc.feature.replace.na(
  ramclustObj = NULL,
  replace.int = 0.1,
  replace.noise = 0.1,
  replace.zero = TRUE,
  which.data = c("MSdata", "MSMSdata")
)
}
\arguments{
\item{ramclustObj}{ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)}

\item{replace.int}{default = 0.1.  proportion of minimum feature value to replace NA (or zero) values with}

\item{replace.noise}{default = 0.1.  proportion ofreplace.int value by which noise is added via 'jitter'}

\item{replace.zero}{logical if TRUE, any zero values are replaced with noise as if they were NA values}

\item{which.data}{name of dataset}
}
\value{
ramclustR object with NA and zero values removed.
}
\description{
replaces any NA (and optionally zero) values with small signal (20% of minimum feature signal value + 20% random noise)
}
\details{
noise is added by finding for each feature the minimum detected value, multiplying that value by replace.int, then adding (replace.int*replace.noise) noise.  abs() is used to ensure no negative values result.
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
