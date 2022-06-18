#' flowPeakDetection
#'
#' The wrapper function flowPeakDetection is the main function that consists of
#' five main peak algorithm functions, as well as five helper functions used within
#' the main functions. The output of flowPeakDetection is a single .csv file with
#' information about the location of each peak and their height, as well as a
#' “messy” and doublet indicator. As well as a statistical measure of confidence
#' for the number of sub-populations.
#'
#' @param xVariable The fluorescence channel on the x axis
#' @param flowDir The directory of the gated .fcs data
#' @param doublet T/F for having the information about doublets in the output
#'  dataset - by default FALSE
#' @param saveGraph T/F for saving the graphs as an output of the NLS
#' @param singleThreshold threshold for classifying single populations
#' @param usedCellsThreshold threshold for classifying multiple populations
#' 
#' @import magrittr
#' @import scorepeak
#' @import data.table
#' @import tcltk
#' @import ggplot2
#' 
#' @export
#'
#' @examples
#' flowPeakDetection(
#'  xVariable = "FITC-A",
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  doublet = FALSE,
#'  saveGraph = TRUE,
#'  singleThreshold = 8,
#'  usedCellsThreshold = 86
#'  )
#'
flowPeakDetection = function(xVariable, flowDir = NA, doublet = FALSE, saveGraph = TRUE, singleThreshold = 8, usedCellsThreshold = 86){

  ##Progress bar iterations
  total <- 7
  ##Create progress bar
  pb <- txtProgressBar(min=0, max=total, style=3)

  ##Directory: If user do not hard codes their directory in the function
  ##a window will open and ask the user to pick their directory.
  ##This should be where the data is located
  if(is.na(flowDir)){
    getwd()
    flowDir <- tclvalue(tkchooseDirectory())
  }

  flowSet <- list.files(flowDir)

  ##Creating a log that will be saved in the global environment
  ##This will tell you which algorithm ran and which flow frames were included as
  ##well as an indicator of success/failure. If the algorithm fails, the user
  ##will know which algorithm and flow framed caused it to fail
  .GlobalEnv$logDs <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logDs) <- c("Algorithm", "Data", "Success")

  ##Running first peak algorithm
  peakAlg1 <- peakAlgorithm1(flowDir, flowSet,  xVariable, singleThreshold)

  ##Getting the flagged flow frames that will be passed to the next algorithm
  flaggedData <- as.data.frame(peakAlg1[[1]])
  names(flaggedData)[1] <- "data"
  ##The data that is not flagged
  singleData <- peakAlg1[2]

  ##update progress bar
  setTxtProgressBar(pb, 1)

  ##Checking if there are any flagged flow frames from the first peak algorithm
  ##If so, the second peak algorithm will run
  if( nrow(peakAlg1[[1]]) != 0 ){
    peakAlg2 <- peakAlgorithm2(
      flowDir,
      flaggedData,
      xVariable,
      usedCellsThreshold
      )
    flaggedData <- as.data.frame(peakAlg2[[1]])
    names(flaggedData)[1] <- "data"
    finishedData <- peakAlg2[2]
    ##update progress bar
    setTxtProgressBar(pb, 2)
  }else{
    ##update progress bar
    setTxtProgressBar(pb, 7)
    finishedData <- NA
    messyData <- NA
    outputData(
      flowDir,
      singleData,
      finishedData,
      messyData,
      xVariable,
      doublet,
      saveGraph
      )
    print("Done! - Check 'analysis' folder for results")
    break
  }


  ##Checking if there are any flagged flow frames from the second peak algorithm
  ##If so, the third peak algorithm will run
  if( nrow(peakAlg2[[1]]) != 0 ){
    peakAlg3 <- peakAlgorithm3(
      flowDir,
      flaggedData,
      xVariable,
      finishedData,
      usedCellsThreshold
      )
    flaggedData <- as.data.frame(peakAlg3[[1]])
    names(flaggedData)[1] <- "data"
    finishedData <- peakAlg3[2]
    ##update progress bar
    setTxtProgressBar(pb, 3)
  }else{
    messyData <- NA
    ##update progress bar
    setTxtProgressBar(pb, 7)
    outputData(
      flowDir,
      singleData,
      finishedData,
      messyData,
      xVariable,
      doublet,
      saveGraph
      )
    print("Done! - Check 'analysis' folder for results")
    break
  }

  ##Checking if there are any flagged flow frames from the third peak algorithm
  ##If so, the forth peak algorithm will run
  if( nrow(peakAlg3[[1]]) != 0 ){
    peakAlg4 <- peakAlgorithm4(
      flowDir,
      flaggedData,
      xVariable,
      finishedData,
      usedCellsThreshold
      )
    flaggedData <- as.data.frame(peakAlg4[[1]])
    messyDataNoNA <- flaggedData %>% dplyr::filter(!is.na(x))
    messyDataNoNA <- messyDataNoNA %>% dplyr::select(-propCellsUsed)
    finishedData <- peakAlg4[2]
    ##update progress bar
    setTxtProgressBar(pb, 4)
  }else{
    messyData <- NA
    ##update progress bar
    setTxtProgressBar(pb, 7)
    outputData(
      flowDir,
      singleData,
      finishedData,
      messyData,
      xVariable,
      doublet,
      saveGraph
      )
    print("Done! - Check 'analysis' folder for results")
    break
  }

  ##Checking if there are any flow frames that have not been analyzed
  ##If so, the last peak algorithm will run
  if(TRUE %in% is.na(peakAlg4[[1]]$x)){
    peakAlg5 <- peakAlgorithm5(flowDir, flaggedData, xVariable, messyDataNoNA)
    messyData <- peakAlg5
  }else{
    messyData <- messyDataNoNA
    outputData(
      flowDir,
      singleData,
      finishedData,
      messyData,
      xVariable,
      doublet,
      saveGraph
      )
    setTxtProgressBar(pb, 7)
  }
  setTxtProgressBar(pb, 5)
  ##Making the dataset
  outputData(
    flowDir,
    singleData,
    finishedData,
    messyData,
    xVariable,
    doublet,
    saveGraph
    )
  ##update progress bar
  setTxtProgressBar(pb, 7)
  print("Done! - Check 'analysis' folder for results")

  ##Closing progress bar
  close(pb)

}

