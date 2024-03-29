% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rectGateFlowFrame.R
\name{rectGateFlowFrame}
\alias{rectGateFlowFrame}
\title{rectGateFlowFrame}
\usage{
rectGateFlowFrame(
  rawDir = NA,
  flowName,
  xVariable = "FL1-A",
  yVariable = "SSC-A",
  xMinValue = 10000,
  xMaxValue = 9e+05,
  yMinValue = 10000,
  yMaxValue = 9e+05,
  savePlot = TRUE
)
}
\arguments{
\item{rawDir}{The directory of the raw .fcs data}

\item{flowName}{The name of the .fcs file}

\item{xVariable}{The fluorescence channel on the x axis}

\item{yVariable}{The fluorescence channel on the y axis}

\item{xMinValue}{The lower bound x value for the gate}

\item{xMaxValue}{The upper bound x value for the gate}

\item{yMinValue}{The lower bound y value for the gate}

\item{yMaxValue}{The upper bound y value for the gate}

\item{savePlot}{A side by side graph comparison the raw data and the gated
data - by default TRUE}
}
\value{
A .fcs of the gated data and plots of gated data
}
\description{
This function will gate your raw .fcs file and save your gated data, a plot
of your gated data, and a plot of the original data
}
\examples{
rectGateFlowFrame(
 rawDir = paste0(system.file(package = "PloidyPeaks"), "/raw_data/"),
 flowName = "A01-A01",
 xVariable = "FITC-A",
 yVariable = "SSC-A",
 xMinValue = 50,
 xMaxValue = 850,
 yMinValue = 50,
 yMaxValue = 850,
 savePlot = TRUE
)
}
