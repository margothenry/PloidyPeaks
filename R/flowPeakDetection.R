#' flowPeakDetection
#'
#' The wrapper function flowPeakDetection is the main function that consists of
#' five main peak algorithm functions, as well as five helper functions used within
#' the main functions. The output of flowPeakDetection is a single .csv file with
#' information about the location of each peak and their height, as well as a
#' “for review” and doublet indicator. In addition a statistical measure of confidence
#' for the number of sub-populations (RSE).
#'
#' @param xVariable The fluorescence channel on the x axis
#' @param flowDir The directory of the gated .fcs data
#' @param doublet T/F for having the information about doublets in the output
#'  dataset - by default FALSE
#' @param saveGraph T/F for saving the graphs as an output of the NLS
#' @param singleThreshold threshold for classifying single populations
#' @param usedCellsThreshold threshold for classifying multiple populations
#' @param MaxDoubletHeight The maximum height a doublet can be. If left as NA
#'  the algorithm will find a value based on the other peaks
#' @param subsetDs A vector of samples that the user wants to analyze if they
#' do not want to analyse the full dataset. Default value is NA
#' 
#' @import magrittr
#' @import scorepeak
#' @import data.table
#' @import tcltk
#' @import ggplot2
#' 
#' @importFrom grDevices dev.off png
#' @importFrom graphics hist
#' @importFrom stats nls nls.control predict quantile
#' @importFrom utils data setTxtProgressBar txtProgressBar write.csv
#' 
#' @export
#'
#' @examples
#' 
#' flowPeakDetection(
#'  xVariable = "FL1-A",
#'  flowDir = here("vignettes/data/gated_data"),
#'  doublet = FALSE,
#'  saveGraph = TRUE,
#'  singleThreshold = 8,
#'  usedCellsThreshold = 86,
#'  MaxDoubletHeight = 50,
#'  subsetDs = c("A01-A01", "A01-A02")
#'  )
#'
flowPeakDetection = function(
  xVariable,
  flowDir = NA,
  doublet = FALSE,
  saveGraph = TRUE,
  singleThreshold = 8,
  usedCellsThreshold = 86,
  MaxDoubletHeight = NA,
  subsetDs = NA
){
  ##Removing NOTE 'no visible binding for global variable'
  x<-propCellsUsed<-NULL
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
  
  if(TRUE %in% is.na(subsetDs)){
    flowSet <- list.files(flowDir)
  }else{
    flowSet <- subset(list.files(flowDir), list.files(flowDir) %in% subsetDs)
  }
  
  
  ##Creating a log that will be saved in the global environment
  ##This will tell you which algorithm ran and which flow frames were included as
  ##well as an indicator of success/failure. If the algorithm fails, the user
  ##will know which algorithm and flow framed caused it to fail
  .GlobalEnv$logDs <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logDs) <- c("Algorithm", "Data", "Success")
  
  ##Running first peak algorithm
  peakAlg1 <- peakAlgorithm1(
    flowDir,
    flowSet,
    xVariable,
    singleThreshold
  )
  
  ##Getting the flagged flow frames that will be passed to the next algorithm
  if(length(peakAlg1) == 2){
    singleData <- peakAlg1[2]
    flaggedData <- as.data.frame(peakAlg1[[1]])
    names(flaggedData)[1] <- "data"
  }else if( length(as.data.frame(peakAlg1[[1]])) == 1){
    singleData <- NULL
    flaggedData <- as.data.frame(peakAlg1[[1]])
    names(flaggedData)[1] <- "data"
  }else{
    flaggedData <- NULL
    singleData <- as.data.frame(peakAlg1[[1]])
  }
  
  
  ##update progress bar
  setTxtProgressBar(pb, 1)
  
  ##Checking if there are any flagged flow frames from the first peak algorithm
  ##If so, the second peak algorithm will run
  if( !is.null(flaggedData) ){
    peakAlg2 <- peakAlgorithm2(
      flowDir,
      flaggedData,
      xVariable,
      usedCellsThreshold,
      MaxDoubletHeight
    )
    
    if(length(peakAlg2) == 2){
      finishedData <- peakAlg2[2]
      flaggedData <- as.data.frame(peakAlg2[[1]])
      names(flaggedData)[1] <- "data"
    }else if( length(as.data.frame(peakAlg2[[1]])) == 1){
      finishedData <- NULL
      flaggedData<-as.data.frame(peakAlg2[[1]])
      names(flaggedData)[1] <- "data"
    }else{
      flaggedData <- NULL
      finishedData<-peakAlg2[[1]]
    }
    
    ##update progress bar
    setTxtProgressBar(pb, 2)
  }else{
    ##update progress bar
    finishedData <- NULL
    messyData <- NULL
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
    print("Done! - Check 'analysis' folder for results")
  }
  
  
  ##Checking if there are any flagged flow frames from the second peak algorithm
  ##If so, the third peak algorithm will run
  if( !is.null(flaggedData) ){
    peakAlg3 <- peakAlgorithm3(
      flowDir,
      flaggedData,
      xVariable,
      finishedData,
      usedCellsThreshold,
      MaxDoubletHeight
    )
    
    if(length(peakAlg3) == 2){
      finishedData <- peakAlg3[2]
      flaggedData <- as.data.frame(peakAlg3[[1]])
      names(flaggedData)[1] <- "data"
    }else if( length(as.data.frame(peakAlg3[[1]])) == 1){
      finishedData <- finishedData
      flaggedData<-as.data.frame(peakAlg3[[1]])
      names(flaggedData)[1] <- "data"
    }else{
      flaggedData <- NULL
      finishedData<-peakAlg3[[1]]
    }
    ##update progress bar
    setTxtProgressBar(pb, 3)
  }else{
    messyData <- NULL
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
  }
  
  ##Checking if there are any flagged flow frames from the third peak algorithm
  ##If so, the forth peak algorithm will run
  if( !is.null(flaggedData) ){
    peakAlg4 <- peakAlgorithm4(
      flowDir,
      flaggedData,
      xVariable,
      finishedData,
      usedCellsThreshold,
      MaxDoubletHeight
    )
    if(length(peakAlg4) == 2){
      flaggedData <- as.data.frame(peakAlg4[[1]])
      messyDataNoNA <- flaggedData %>% dplyr::filter(!is.na(x))
      messyDataNoNA <- messyDataNoNA %>% dplyr::select(-propCellsUsed)
      finishedData <- peakAlg4[2]
    }else if( length(as.data.frame(peakAlg4[[1]])) == 1){
      finishedData <- finishedData
      flaggedData<-as.data.frame(peakAlg4[[1]])
      names(flaggedData)[1] <- "data"
    }else{
      flaggedData <- NULL
      messyDataNoNA <- NULL
      finishedData <- peakAlg4[[1]]
    }
    
    ##update progress bar
    setTxtProgressBar(pb, 4)
  }else{
    messyData <- NULL
    ##update progress bar
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
    print("Done! - Check 'analysis' folder for results")
  }
  
  ##Checking if there are any flow frames that have not been analyzed
  ##If so, the last peak algorithm will run
  if(TRUE %in% is.na(flaggedData$x)){
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
  print("Done! - Check 'analysis' folder for results")
  
  ##Closing progress bar
  close(pb)
  
}

