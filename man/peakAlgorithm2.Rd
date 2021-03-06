% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakAlgorithm2.R
\name{peakAlgorithm2}
\alias{peakAlgorithm2}
\title{peakAlgorithm2}
\usage{
peakAlgorithm2(
  flowDir,
  flaggedData_,
  xVariable,
  usedCellsThreshold,
  MaxDoubletHeight
)
}
\arguments{
\item{flowDir}{The directory of the gated .fcs data}

\item{flaggedData_}{List of names of the flow frames that got flagged by PeakAlgorithm1}

\item{xVariable}{The fluorescence channel on the x axis}

\item{usedCellsThreshold}{Threshold for classifying multiple populations}

\item{MaxDoubletHeight}{The maximum height a doublet can be. If left as NA
the algorithm will find a value based on the other peaks}
}
\description{
The second branching point is in peakAlgorithm2, which seeks to identify all
possible sub populations. The algorithm looks at the proportion of cells used
in the subpopulations (CellProp). In other words, once the populations have been
identified for a given flow frame, peakAlgorithm2 will find the cells associated with the G1/G2 pairs.
}
\examples{
peakAlgorithm2(
 flowDir = "FlowData/T10_FLC/gated_data",
 flaggedData_ = flaggedData ,
 xVariable = "FITC-A",
 usedCellsThreshold = 86,
 MaxDoubletHeight = 50
 )

}
