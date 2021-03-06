% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowPeakDetection.R
\name{flowPeakDetection}
\alias{flowPeakDetection}
\title{flowPeakDetection}
\usage{
flowPeakDetection(
  xVariable,
  flowDir = NA,
  doublet = FALSE,
  saveGraph = TRUE,
  singleThreshold = 8,
  usedCellsThreshold = 86,
  MaxDoubletHeight = NA,
  subsetDs = NA
)
}
\arguments{
\item{xVariable}{The fluorescence channel on the x axis}

\item{flowDir}{The directory of the gated .fcs data}

\item{doublet}{T/F for having the information about doublets in the output
dataset - by default FALSE}

\item{saveGraph}{T/F for saving the graphs as an output of the NLS}

\item{singleThreshold}{threshold for classifying single populations}

\item{usedCellsThreshold}{threshold for classifying multiple populations}

\item{MaxDoubletHeight}{The maximum height a doublet can be. If left as NA
the algorithm will find a value based on the other peaks}

\item{subsetDs}{A vector of samples that the user wants to analyze if they
do not want to analyse the full dataset. Default value is NA}
}
\description{
The wrapper function flowPeakDetection is the main function that consists of
five main peak algorithm functions, as well as five helper functions used within
the main functions. The output of flowPeakDetection is a single .csv file with
information about the location of each peak and their height, as well as a
“for review” and doublet indicator. In addition a statistical measure of confidence
for the number of sub-populations (RSE).
}
\examples{

flowPeakDetection(
 xVariable = "FL1-A",
 flowDir = here("vignettes/data/gated_data"),
 doublet = FALSE,
 saveGraph = TRUE,
 singleThreshold = 8,
 usedCellsThreshold = 86,
 MaxDoubletHeight = 50,
 subsetDs = c("A01-A01", "A01-A02")
 )

}
