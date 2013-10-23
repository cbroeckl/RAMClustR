\name{ramclustR}
\alias{ramclustR}
\title{ramclustR}
\usage{
  ramclustR(ms = MSdata, idmsms = ms, idMSMStag = "02.cdf",
    featdelim = "_", timepos = 2, st = 5, sr = 5,
    maxt = 20, deepSplit = FALSE, blocksize = 1000,
    mult = 10)
}
\arguments{
  \item{filename}{character Filename of the nmrML to check}

  \item{ms}{MS1 intensities =MSdata,}

  \item{idmsms}{=ms,}

  \item{idMSMStag}{character e.g. "02.cdf"}

  \item{featdelim}{character e.g. ="_"}

  \item{timepos}{numeric 2}

  \item{st}{numeric no clue e.g. =5,}

  \item{sr}{numeric also no clue yet =5,}

  \item{maxt}{numeric again no clue =20,}

  \item{deepSplit}{boolean e.g. =FALSE,}

  \item{blocksize}{integer number of features (scans?)
  processed in one block =1000,}

  \item{mult}{numeric =10}
}
\value{
  A vector with the numeric values of the processed data
}
\description{
  Main clustering function
}
\details{
  This is the Details section
}
\author{
  Corey Broeckling
}

