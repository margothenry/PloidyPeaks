#' flowPeakDetection
#'
#' The wrapper function flowPeakDetection is the main function that consists of
#' five main peak algorithm functions, as well as five helper functions used
#' within the main functions. The output of flowPeakDetection is a single .csv 
#' file with information about the location of each peak and their height, 
#' as well as a “for review” and doublet indicator. In addition a statistical
#' measure of confidence for the number of sub-populations (RSE).
#'
#' @param xVariable The fluorescence channel on the x axis
#' @param flowDir The directory of the gated .fcs data
#' @param doublet T/F for having the information about doublets in the output
#'  dataset - by default FALSE
#' @param saveGraph T/F for saving the graphs as an output of the NLS
#' @param singleThreshold threshold for classifying single populations
#' @param usedCellsThreshold threshold for classifying multiple populations
#' @param maxDoubletHeight The maximum height a doublet can be. If left as NA
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
#' @return a .csv with information about each sample and nls graphs
#' @export
#'
#' @examples
#' flowPeakDetection(
#'  xVariable = "FITC-A",
#'  flowDir = paste0(system.file(package = "PloidyPeaks"), "/gated_data/"),
#'  doublet = FALSE,
#'  saveGraph = FALSE,
#'  singleThreshold = 8,
#'  usedCellsThreshold = 86,
#'  maxDoubletHeight = NA,
#'  subsetDs = c("A01-A02", "A04-D12", "A07-G12", "A09-A02", "T1-D08")
#'  )
#'  
flowPeakDetection = function(
  xVariable = "FL1-A",
  flowDir = NA,
  doublet = FALSE,
  saveGraph = TRUE,
  singleThreshold = 8,
  usedCellsThreshold = 86,
  maxDoubletHeight = NA,
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
    flowSet <- subset(
      list.files(flowDir),
      list.files(flowDir) %in% subsetDs
    )
  }
  
  
  ##Creating a log that will be saved in the global environment
  ##This will tell you which algorithm ran and which flow frames were included
  ##as well as an indicator of success/failure. If the algorithm fails, the
  ##user will know which algorithm and flow framed caused it to fail
  .GlobalEnv$logDs <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logDs) <- c("Algorithm", "Data", "Success")
  
  ##Running first peak algorithm
  peakAlg1 <- .peakAlgorithm1(
    flowDir,
    flowSet,
    xVariable,
    singleThreshold
  )
  
  ##Getting the flagged flow frames that will be passed to the next algorithm
  if(length(peakAlg1) == 3){
    singleData <- peakAlg1[[2]]
    flaggedData <- as.data.frame(peakAlg1[[1]])
    names(flaggedData)[1] <- "data"
    PACheck <- as.data.frame(peakAlg1[[3]])
  }else if(length(peakAlg1) == 2){
    if(length(as.data.frame(peakAlg1[[1]])) == 1){
      flaggedData <- as.data.frame(peakAlg1[[1]])
      names(flaggedData)[1] <- "data"
      PACheck <- NULL
    }else if(length(as.data.frame(peakAlg1[[1]])) == 2){
      PACheck <- as.data.frame(peakAlg1[[1]])
    }else{
      singleData <- peakAlg1[[1]]
    }
    
    if(length(as.data.frame(peakAlg1[[2]])) == 1){
      flaggedData <- as.data.frame(peakAlg1[[2]])
      names(flaggedData)[1] <- "data"
      PACheck <- NULL
    }else if(length(as.data.frame(peakAlg1[[2]])) == 2){
      PACheck <- as.data.frame(peakAlg1[[2]])
    }else{
      singleData <- peakAlg1[[2]]
    }
    
  }else if(length(peakAlg1) == 1){
    if(length(as.data.frame(peakAlg1[[1]])) == 1){
      singleData <- NULL
      flaggedData <- as.data.frame(peakAlg1[[1]])
      names(flaggedData)[1] <- "data"
      PACheck <- NULL
    }else if(length(as.data.frame(peakAlg1[[1]])) == 2){
      singleData <- NULL
      flaggedData <- NULL
      PACheck <- as.data.frame(peakAlg1[[1]])
    }else{
    flaggedData <- NULL
    singleData <- as.data.frame(peakAlg1[[1]])
    PACheck <- NULL
    }
  }
  
  ##update progress bar
  setTxtProgressBar(pb, 1)
  
  ##Checking if there are any flagged flow frames from the first peak 
  ##algorithm. If so, the second peak algorithm will run
  if( !is.null(flaggedData) ){
    peakAlg2 <- .peakAlgorithm2(
      flowDir,
      flaggedData,
      xVariable,
      usedCellsThreshold,
      maxDoubletHeight
    )
    
    if(length(peakAlg2) == 2){
      finishedData <- peakAlg2[[2]]
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
    investigateData <- NULL
    .outputData(
      flowDir,
      singleData,
      finishedData,
      investigateData,
      PACheck,
      xVariable,
      doublet,
      saveGraph
    )
    setTxtProgressBar(pb, 7)
    print("Done! - Check 'analysis' folder for results")
  }
  
  #updating the PA flags
  if(!purrr::is_empty(PACheck) & !purrr::is_empty(finishedData)){
    subsetCheck <- finishedData %>% dplyr::filter(
      data %in% PACheck$data
    )
    if(nrow(subsetCheck) > 0){
      PACheck$flag[subsetCheck$data %in% PACheck$data] = 0
    }
  }
  
  ##Checking if there are any flagged flow frames from the second peak 
  ##algorithm. If so, the third peak algorithm will run
  if( !is.null(flaggedData) ){
    peakAlg3 <- .peakAlgorithm3(
      flowDir,
      flaggedData,
      xVariable,
      finishedData,
      usedCellsThreshold,
      maxDoubletHeight
    )
    
    if(length(peakAlg3) == 2){
      finishedData <- peakAlg3[[2]]
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
    investigateData <- NULL
    ##update progress bar
    setTxtProgressBar(pb, 7)
    .outputData(
      flowDir,
      singleData,
      finishedData,
      investigateData,
      PACheck,
      xVariable,
      doublet,
      saveGraph
    )
    print("Done! - Check 'analysis' folder for results")
  }
  
  ##Checking if there are any flagged flow frames from the third peak 
  ##algorithm. If so, the forth peak algorithm will run
  if( !is.null(flaggedData) ){
    peakAlg4 <- .peakAlgorithm4(
      flowDir,
      flaggedData,
      xVariable,
      finishedData,
      usedCellsThreshold,
      maxDoubletHeight
    )
    if(length(peakAlg4) == 2){
      flaggedData <- as.data.frame(peakAlg4[[1]])
      investigateDataNoNA <- flaggedData %>% dplyr::filter(!is.na(x))
      investigateDataNoNA <- investigateDataNoNA %>%
        dplyr::select(-propCellsUsed)
      finishedData <- peakAlg4[[2]]
    }else if( length(as.data.frame(peakAlg4[[1]])) == 1){
      finishedData <- finishedData
      flaggedData<-as.data.frame(peakAlg4[[1]])
      names(flaggedData)[1] <- "data"
    }else{
      flaggedData <- NULL
      investigateDataNoNA <- NULL
      finishedData <- peakAlg4[[1]]
    }
    
    ##update progress bar
    setTxtProgressBar(pb, 4)
  }else{
    investigateData <- NULL
    ##update progress bar
    R.utils::doCall(
      .outputData(
        flowDir,
        singleData,
        finishedData,
        investigateData,
        PACheck,
        xVariable,
        doublet,
        saveGraph
      )
    ) 
    setTxtProgressBar(pb, 7)
    print("Done! - Check 'analysis' folder for results")
  }
  
  ##Checking if there are any flow frames that have not been analyzed
  ##If so, the last peak algorithm will run
  if(TRUE %in% is.na(flaggedData$x)){
    peakAlg5 <- .peakAlgorithm5(
      flowDir,
      flaggedData,
      xVariable,
      investigateDataNoNA
    )
    investigateData <- peakAlg5
  }else{
    investigateData <- investigateDataNoNA
    .outputData(
      flowDir,
      singleData,
      finishedData,
      investigateData,
      PACheck,
      xVariable,
      doublet,
      saveGraph
    )
    setTxtProgressBar(pb, 7)
  }
  setTxtProgressBar(pb, 5)
  ##Making the dataset
  .outputData(
    flowDir,
    singleData,
    finishedData,
    investigateData,
    PACheck,
    xVariable,
    doublet,
    saveGraph
  )
  ##update progress bar
  print("Done! - Check 'analysis' folder for results")
  
  ##Closing progress bar
  close(pb)
  
}
## Helper functions within the wrapper function

## peakAlgorithm1
.peakAlgorithm1 = function(
  flowDir,
  flowSet,
  xVariable,
  singleThreshold = 8,
  missedThreshold = 6
){
  
  singleData <- c()
  flaggedData <- c()
  checkFlag <- c()
  ##Adding flow frame to log
  logFlow <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  algorithmNum <- 1
  
  ##Looping through each flow frame
  for(k in seq_len(length(flowSet))){
    ##Checking the data in the folder is in proper format
    if(
      tools::file_ext(flowSet[k]) %in% c(
        "csv", "xls", "xlsx","html", "ppt",
        "pptx" ,"pdf", "doc", "docx"
      )
    ){
      errorMsg<-paste0(
        "The flow frame ",flowSet[k]," is does not seem 
                                to be in the valid flow format. Consider 
                                changing the format or removing from the folder"
      )
      stop(errorMsg)
      rm(errorMsg)
    }
    
    ##Reading in and smoothing data
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowSet[k]), transformation=FALSE
    )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    flowData <- .smoothData( flowName, xVariable, 11)
    
    logFlow[1, ] <- c(
      algorithmNum,
      flowName@description[["GUID"]],
      0
    )
    
    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )
    
    ##Finding local peaks
    localPeaks <- scorepeak::detect_localmaxima(flowData$y, 3)
    possiblePeaks <- flowData[localPeaks, ]
    
    ##Removing the peaks that are identified at the base of the
    ##histogram
    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > quantile(flowData$y)[4]),
    ]
    xVarMax <- max(flowData$x)
    ##Removing the peaks that are identified at the extreme left 
    ##side of the histogram that could be caused by debris/improper
    ##gating
    possiblePeaks2 <- possiblePeaks2[
      which(possiblePeaks2$x > xVarMax/9.5), ]
    
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > quantile(flowData$y)[4]),
      ]
    }
    
    possiblePeaks3 <- .findTruePeaks(possiblePeaks2, 40, xVarMax)
    
    ##Ordering the peaks
    possiblePeaks4 <- possiblePeaks3[
      order(possiblePeaks3$x, decreasing=FALSE),]
    ##Removing the peaks that are identified at the extreme
    ## right side of the histogram
    possiblePeaks4 <- possiblePeaks4[
      which(possiblePeaks4$x < quantile(flowData$x)[4]), ]
    peaksFix <- possiblePeaks4
    ##Identifying the first peak in the dataset
    firstPeak <- possiblePeaks4[1, seq_len(2)]
    ##If there is only one peak identified, try to identify the pair
    ##for this peak
    if(nrow(possiblePeaks4) == 1){
      
      initialPeak <- firstPeak
      ##The peak is G1 and G2 is missing
      g2PeakRadiusUL <- initialPeak$x*2.2 + 1
      g2PeakRadiusLL <- initialPeak$x*1.75 - 1
      g2ToTestDs2 <- which(
        g2PeakRadiusLL < possiblePeaks$x & 
          possiblePeaks$x < g2PeakRadiusUL &
          possiblePeaks$x < quantile(flowData$x)[4]
      )
      g2ToTestDs3 <- possiblePeaks[g2ToTestDs2,]
      if(nrow(g2ToTestDs3) == 0){
        g2ToTestDs3 <- data.frame(
          x=0,
          y=0
        )
      }
      #The peak is G2 and G1 is missing
      g1PeakRadiusLL <- initialPeak$x/2.2 + 1
      g1PeakRadiusUL <- initialPeak$x/1.75 - 1
      g1ToTestDs2 <- which(
        g1PeakRadiusLL < possiblePeaks$x & 
          possiblePeaks$x < g1PeakRadiusUL
      )
      g1ToTestDs3 <- possiblePeaks[g1ToTestDs2,]
      if(nrow(g1ToTestDs3) == 0){
        g1ToTestDs3 <- data.frame(
          x=0,
          y=0
        )
      }
      
      if(g2ToTestDs3$y[1] < g1ToTestDs3$y[1]){
        peak2 <- g1ToTestDs3[1,]
        peak1 <- initialPeak
        peaksFix <- rbind(
          peak1,
          peak2
        )
      }else{
        peak2 <- g2ToTestDs3[1,]
        peak1 <- initialPeak
        peaksFix <- rbind(
          peak1,
          peak2
        )
      }
    }
    
    ##checking if we have the proper G1 peak
    peaksFix <- peaksFix[order(peaksFix$x, decreasing=FALSE), ]
    firstPeak <- peaksFix[1, seq_len(2)]
    if(firstPeak$x > xVarMax/3.1 | firstPeak$x < xVarMax/6){
      
      initialPeak <- firstPeak
      otherPeak <- peaksFix[-1, seq_len(2)]
      
      ##the peak is G2 missing G1
      g1PeakRadiusLL <- initialPeak$x/2.3 - 1
      g1PeakRadiusUL <- initialPeak$x/1.7 + 1
      g1ToTestDs2 <- which(
        g1PeakRadiusLL <= possiblePeaks$x & 
          possiblePeaks$x <= g1PeakRadiusUL
      )
      g1ToTestDs3 <- possiblePeaks[g1ToTestDs2, ]
      
      if(nrow(g1ToTestDs3) == 0){
        g1ToTestDs2 <- which(
          g1PeakRadiusLL <= flowData$x & 
            flowData$x <= g1PeakRadiusUL
        )
        g1ToTestDs3 <- flowData[g1ToTestDs2, ]
      }
      g1ToTestDs3 <- g1ToTestDs3[
        order(g1ToTestDs3$y, decreasing=TRUE), ]
      g1ToTestDs3 <- g1ToTestDs3[
        which(g1ToTestDs3$y > quantile(flowData$y)[2]),]
      peak2 <- g1ToTestDs3[1, ]
      peak1 <- initialPeak
      peaksFix <- rbind(
        peak1,
        peak2,
        otherPeak
      )
    }
    
    ##Finding the are to the right of the G2 peak
    peaksFix <- peaksFix[order(peaksFix$x, decreasing=FALSE), ]
    peaksFix <- peaksFix %>% dplyr::distinct()
    
    ##Finding the distance between the G1 and G2 peak
    
    epsilon <- abs((peaksFix$x[2] - peaksFix$x[1])/2)
    epsilonLeft <- peaksFix$x[2] - epsilon
    epsilonRight <- peaksFix$x[2] + epsilon
    
    xEpsilon <- which(
      abs(
        flowData$x-epsilonRight) == 
        min(abs(flowData$x-epsilonRight)
        )
    )
    if(length(xEpsilon)>1){
      xEpsilon <- xEpsilon[1]
    }
    
    xEpsilon2 <- which(
      abs(
        flowData$x-epsilonLeft) == 
        min(abs(flowData$x-epsilonLeft)
        )
    )
    if(length(xEpsilon2)>1){
      xEpsilon2 <- xEpsilon2[1]
    }
    
    postEpsilon <- flowData$y[c(xEpsilon:length(flowData$x))]
    preEpsilon <- flowData$y[c(0:xEpsilon2)]
    cellCountPostEpsilon <- sum(postEpsilon)
    cellCountPreEpsilon <- sum(preEpsilon)
    totalCellCount <- sum(flowData$y)
    ##Proportion of cells that are to the right of the G2 peak
    detectR <- round((cellCountPostEpsilon/totalCellCount)*100, 4)
    detectL <- round((cellCountPreEpsilon/totalCellCount)*100, 4)
    
    if(detectR <= singleThreshold & detectL > missedThreshold){
      
      if(peaksFix$x[1] == 0){
        peaksFix[1, ] <- peaksFix[2, ]
        peaksFix$x[2] <- 0
        peaksFix$y[2] <- 0
      }
      
      singleDs <- data.frame(
        data=flowName@description[["GUID"]],
        x=peaksFix$x[1],
        y=peaksFix$y[1],
        possiblePairX=peaksFix$x[2],
        possiblePairY=peaksFix$y[2]
      )
      
      rangeLength <- nchar(format(xVarMax, scientific=FALSE))
      multiplier <- 10^(rangeLength-3)
      singleDs2 <- .doubletCheck(
        singleDs,
        possiblePeaks,
        10*multiplier,
        15*multiplier
      )
      
      singleData <- rbind(
        singleData,
        singleDs2
      )
      singleDataUpdated <- .updatedMeans(
        singleData,
        flowDir,
        xVariable
      )
    }else{
      
      flaggedData <- rbind(
        flaggedData,
        flowName@description[["GUID"]]
      )
    }
    
    if(detectL <= missedThreshold | detectL >= 85){
      checkPeakFlag <- data.frame(
        data=flowName@description[["GUID"]],
        flag=1
      )
      
      checkFlag <- rbind(
        checkFlag,
        checkPeakFlag
      )
    }
    
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
    
  }
  
  if(!purrr::is_empty(singleData) & !purrr::is_empty(flaggedData)){
    if(!purrr::is_empty(checkFlag)){
      returnedList <- list(flaggedData, singleDataUpdated, checkFlag)
    }else{
      returnedList <- list(flaggedData, singleDataUpdated)
    }
  }else if(!purrr::is_empty(singleData) & purrr::is_empty(flaggedData)){
    if(!purrr::is_empty(checkFlag)){
      returnedList <- list(singleData, checkFlag) 
    }else{
      returnedList <- list(singleData)      
    }
  }else{
    if(!purrr::is_empty(checkFlag)){
      returnedList <- list(flaggedData, checkFlag) 
    }else{
      returnedList <- list(flaggedData) 
    }
    
  }
  
  return(returnedList)
}

## peakAlgorithm2
.peakAlgorithm2 = function(
  flowDir,
  flaggedData_,
  xVariable,
  usedCellsThreshold,
  maxDoubletHeight
){
  ##Removing NOTE 'no visible binding for global variable'
  y<-possiblePairY<-keep<-truePeak<-NULL
  
  flowNameDs <- flaggedData_$data
  finishedData <- c()
  flaggedData <- c()
  
  logFlow <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  algorithmNum <- 2
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    flowData <- .smoothData( flowName, xVariable, 5)
    
    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)
    
    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )
    
    localPeaks <- scorepeak::detect_localmaxima(flowData$y, 3)
    possiblePeaks <- flowData[localPeaks, ]
    
    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/3),
    ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/10), ]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > max(possiblePeaks$y)/4),
      ]
    }
    possiblePeaks3 <- .findTruePeaks(possiblePeaks2, 40, xVarMax)
    
    possiblePeaks4 <- .findPairs(possiblePeaks3, possiblePeaks3, 1.75, 2.2)
    
    possiblePeaks4 <- possiblePeaks4 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        possiblePairY > quantile(flowData$y)[4]
      )
    possiblePeaks5 <- possiblePeaks4[
      !duplicated(possiblePeaks4$possiblePairY),]
    
    if(nrow(possiblePeaks5) > 1){
      tempDs<-possiblePeaks5
      tempDs2<-c()
      for(i in seq_len(nrow(tempDs))){
        popInQuestion<-tempDs[i,]
        if(is.na(maxDoubletHeight)){
          maxDoubletHeight_ <- round((popInQuestion$y)/2,2)
        }else{
          maxDoubletHeight_<-maxDoubletHeight
        }
        if(popInQuestion$possiblePairY < maxDoubletHeight_){
          popInQuestion$truePeak = FALSE
        }else{
          popInQuestion$truePeak = TRUE
        }
        tempDs2<-rbind(
          tempDs2,
          popInQuestion
        )
      }
      possiblePeaks5<-tempDs2 %>% 
        dplyr::filter(truePeak == TRUE) %>% 
        dplyr::select(-truePeak)
      rm(tempDs, tempDs2)  
    }
    
    ## Finding Doublets
    rangeLength <- nchar(format(xVarMax, scientific=FALSE))
    multiplier <- 10^(rangeLength-3)
    possiblePeaks6 <- .doubletCheck(
      possiblePeaks5,
      possiblePeaks,
      10*multiplier,
      15*multiplier
    )
    
    ##Checking if G2+G2 doublet is a true doublet and not a G2 peak for
    ##Another subpopulation
    tempDs<-possiblePeaks6
    if(nrow(tempDs) > 1 & !is.na(tempDs$g1G2DoubletCount[1])){
      pop1<-tempDs[1,]
      for(i in 2:nrow(tempDs)){
        popInQuestion<-tempDs[i,]
        if(!is.na(pop1$g2G2Doublet)){
          if(pop1$g2G2Doublet == popInQuestion$possiblePairX & 
             popInQuestion$possiblePairY < maxDoubletHeight_){
            tempDs <- tempDs[-i,]
          }else{
            tempDs$g1G2Doublet<-NA
            tempDs$g2G2Doublet<-NA
            tempDs$g1G2DoubletCount<-NA
            tempDs$g2G2DoubletCount<-NA
          } 
        }
      }
    }
    possiblePeaks6<-tempDs
    rm(tempDs) 
    
    if(nrow(possiblePeaks6) != 0){
      
      if(nrow(possiblePeaks6) == 1){
        peakRow <- data.frame(
          peaks = c(
            possiblePeaks6$x,
            possiblePeaks6$possiblePairX,
            possiblePeaks6$g1G2Doublet,
            possiblePeaks6$g2G2Doublet
          )
        )
      }else{
        peakRow <- data.frame(
          peaks=c(
            possiblePeaks6$x,
            possiblePeaks6$possiblePairX
          )
        )
      }
      
      peakRow <- peakRow %>%
        dplyr::distinct() %>% tidyr::drop_na()
      cellsUsed <- c()
      for(j in seq_len(nrow(peakRow))){
        
        xPeak <- which(flowData$x == peakRow$peaks[j])
        if(length(xPeak)>1){
          xPeak <- xPeak[1]
        }
        
        flowData$keep <- FALSE
        flowData$keep[xPeak] <- TRUE
        
        ##Forwards
        for(i in xPeak:nrow(flowData)){
          if(i == nrow(flowData)){
            flowData$keep[nrow(flowData)] <- TRUE
            break
          }else{
            if(flowData$y[i+1] > flowData$y[i]){
              flowData$keep[i+1] <- TRUE
              break
            }
          }
        }
        
        ##Backwards
        for(i in xPeak:1){
          if(i == 1){
            flowData$keep[1] <- TRUE
            break
          }else{
            if(flowData$y[i-1] > flowData$y[i]){
              flowData$keep[i-1] <- TRUE
              break
            }
          }
        }
        
        intPeak <- flowData %>% dplyr::filter(
          keep == TRUE
        )
        
        xEpsilonRight <- which(flowData$x == intPeak$x[3])
        xEpsilonLeft <- which(flowData$x == intPeak$x[1])
        peakCount <- flowData$y[c(xEpsilonLeft:xEpsilonRight)]
        peakCount01 <- sum(peakCount)
        peakCount01[j] <- peakCount01
        cellsUsed <- c(
          cellsUsed,
          peakCount01[j]
        )
        
      }
      
      cellsUsed01 <- sum(cellsUsed)
      totalCellCount <- sum(flowData$y)
      propCellsUsed <- round((cellsUsed01/totalCellCount)*100, 2)
      possiblePeaks6$propCellsUsed <- propCellsUsed
    }else{
      possiblePeaks6 <- data.frame(
        x=NA,
        y=NA,
        cluster=NA,
        LL=NA,
        UL=NA,
        possiblePairX=NA,
        possiblePairY=NA,
        data=flowName@description[["GUID"]],
        propCellsUsed=NA
      )
    }
    
    if(
      possiblePeaks6$propCellsUsed[1] >= usedCellsThreshold &
      !is.na(possiblePeaks6$propCellsUsed[1])
    ){
      
      possiblePeaks7 <- possiblePeaks6 %>% dplyr::mutate(
        data=flowName@description[["GUID"]]
      )
      
      finishedData <- rbind(
        finishedData,
        possiblePeaks7
      )
      
    }else{
      flaggedData <- rbind(
        flaggedData,
        flowName@description[["GUID"]]
      )
      
    }
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
  }
  
  if(!purrr::is_empty(finishedData) & !purrr::is_empty(flaggedData)){
    returnedList <- list(flaggedData, finishedData)
  }else if(!purrr::is_empty(finishedData) & purrr::is_empty(flaggedData)){
    returnedList <- list(finishedData) 
  }else{
    returnedList <- list(flaggedData) 
  }
  
  return(returnedList)
  
}

## peakAlgorithm3
.peakAlgorithm3 = function(
  flowDir,
  flaggedData_,
  xVariable,
  appendData,
  usedCellsThreshold,
  maxDoubletHeight
){
  
  ##Removing NOTE 'no visible binding for global variable'
  y<-possiblePairY<-keep<-truePeak<-NULL
  
  flowNameDs <- flaggedData_$data
  flaggedData <- c()
  
  logFlow <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  
  algorithmNum <- 3
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    
    flowData <- .smoothData( flowName, xVariable, 5)
    
    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)
    
    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )
    
    localPeaks <- scorepeak::detect_localmaxima(flowData$y, 5)
    possiblePeaks <- flowData[localPeaks, ]
    
    possiblePeaks1 <- possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/5),
    ]
    
    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/3.5),
    ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/10), ]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > max(possiblePeaks$y)/3.5),
      ]
    }
    possiblePeaks3 <- .findTruePeaks(possiblePeaks2, 40, xVarMax)
    
    possiblePeaks4 <- .findPairs(possiblePeaks3, possiblePeaks1, 1.7, 2.3)
    
    possiblePeaks5 <- possiblePeaks4 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        possiblePairY > quantile(flowData$y)[3]+10
      )
    possiblePeaks6 = possiblePeaks5[
      !duplicated(possiblePeaks5$possiblePairY), ]
    
    
    if(nrow(possiblePeaks6) > 1){
      possiblePeaks7 <- .findTruePeaks(possiblePeaks6, 40, xVarMax)
    }else{
      possiblePeaks7 <- possiblePeaks6
    }
    
    ##Checking peak sizes to see if it's a true peak
    if(nrow(possiblePeaks7) > 1){
      tempDs<-possiblePeaks7
      tempDs2<-c()
      for(i in seq_len(nrow(tempDs))){
        popInQuestion<-tempDs[i,]
        if(is.na(maxDoubletHeight)){
          maxDoubletHeight_ <- round((popInQuestion$y)/3,2)
        }else{
          maxDoubletHeight_<-maxDoubletHeight
        }
        if(popInQuestion$possiblePairY < maxDoubletHeight_){
          popInQuestion$truePeak = FALSE
        }else{
          popInQuestion$truePeak = TRUE
        }
        tempDs2<-rbind(
          tempDs2,
          popInQuestion
        )
      }
      possiblePeaks7<-tempDs2 %>% 
        dplyr::filter(truePeak == TRUE) %>% 
        dplyr::select(-truePeak)
      rm(tempDs, tempDs2)  
    }
    
    ## Finding Doublets
    rangeLength <- nchar(format(xVarMax, scientific=FALSE))
    multiplier <- 10^(rangeLength-3)
    possiblePeaks8 <- .doubletCheck(
      possiblePeaks7,
      possiblePeaks,
      10*multiplier,
      15*multiplier
    )
    
    ##Checking if G2+G2 doublet is a true doublet and not a G2 peak for
    ##Another subpopulation
    tempDs<-possiblePeaks8
    if(nrow(tempDs) > 1 & !is.na(tempDs$g1G2DoubletCount[1])){
      pop1<-tempDs[1,]
      for(i in 2:nrow(tempDs)){
        popInQuestion = tempDs[i,]
        if(pop1$g2G2Doublet == popInQuestion$possiblePairX & 
           popInQuestion$possiblePairY < maxDoubletHeight_){
          tempDs <- tempDs[-i,]
        }else{
          tempDs$g1G2Doublet<-NA
          tempDs$g2G2Doublet<-NA
          tempDs$g1G2DoubletCount<-NA
          tempDs$g2G2DoubletCount<-NA
        }
      }
    }
    possiblePeaks8<-tempDs
    rm(tempDs) 
    
    
    if(nrow(possiblePeaks8) != 0){
      if(nrow(possiblePeaks8) == 1){
        peakRow <- data.frame(
          peaks=c(
            possiblePeaks8$x,
            possiblePeaks8$possiblePairX,
            possiblePeaks8$g1G2Doublet,
            possiblePeaks8$g2G2Doublet
          )
        )
      }else{
        peakRow <- data.frame(
          peaks = c(possiblePeaks8$x, possiblePeaks8$possiblePairX)
        )
      }
      
      peakRow <- peakRow %>%
        dplyr::distinct() %>% tidyr::drop_na()
      
      cellsUsed <- c()
      
      smoothedData <- .smoothData(flowName, xVariable, 10)
      
      for(j in seq_len(nrow(peakRow))){
        xPeak <- which(smoothedData$x == peakRow$peaks[j])
        if(length(xPeak)>1){
          xPeak <- xPeak[1]
        }
        
        if(xPeak+10 > nrow(smoothedData) & xPeak-10 < 1){
          peakRowSmoothedRange <- smoothedData[
            seq(1, nrow(smoothedData), 1), ]
        }else if(xPeak+10 > nrow(smoothedData) & xPeak-10 >= 1){
          peakRowSmoothedRange <- smoothedData[
            seq(xPeak-10, nrow(smoothedData), 1),
          ]
        }else if(xPeak+10 <= nrow(smoothedData) & xPeak-10 < 1){
          peakRowSmoothedRange <- smoothedData[seq(1, xPeak+10, 1), ]
        }else{
          peakRowSmoothedRange <- smoothedData[
            seq(xPeak-10, xPeak+10, 1), ]
        }
        
        xPeakSmoothed <- which(
          max(peakRowSmoothedRange$y) == smoothedData$y
        )
        if(length(xPeakSmoothed) > 1){
          xPeakSmoothed <- xPeakSmoothed[1]
        }
        smoothedData$keep <- FALSE
        smoothedData$keep[xPeakSmoothed] <- TRUE
        
        #Forwards
        for(i in xPeakSmoothed:nrow(smoothedData)){
          if(i == nrow(smoothedData)){
            smoothedData$keep[nrow(smoothedData)] <- TRUE
            break
          }else{
            if(smoothedData$y[i+1] > smoothedData$y[i]){
              smoothedData$keep[i+1] <- TRUE
              break
            }
          }
        }
        
        #Backwards
        for(i in xPeakSmoothed:1){
          if(i == 1){
            smoothedData$keep[1] <- TRUE
            break
          }else{
            if(smoothedData$y[i-1] > smoothedData$y[i]){
              smoothedData$keep[i-1] <- TRUE
              break
            }
          }
        }
        
        intPeak <- smoothedData %>% dplyr::filter(
          keep == TRUE
        )
        
        xEpsilonRight <- which(flowData$x == intPeak$x[3])
        xEpsilonLeft <- which(flowData$x == intPeak$x[1])
        peakCount <- flowData$y[c(xEpsilonLeft:xEpsilonRight)]
        peakCount01 <- sum(peakCount)
        peakCount01[j] <- peakCount01
        cellsUsed <- c(
          cellsUsed,
          peakCount01[j]
        )
      }
      
      cellsUsed01 <- sum(cellsUsed)
      totalCellCount <- sum(smoothedData$y)
      propCellsUsed <- round((cellsUsed01/totalCellCount)*100, 2)
      possiblePeaks8$propCellsUsed <- propCellsUsed
      
    }else{
      possiblePeaks8 <- data.frame(
        x=NA,
        y=NA,
        cluster=NA,
        LL=NA,
        UL=NA,
        possiblePairX=NA,
        possiblePairY=NA,
        data=flowName@description[["GUID"]],
        propCellsUsed=NA
      )
    }
    
    if(
      possiblePeaks8$propCellsUsed[1] >= usedCellsThreshold &
      !is.na(possiblePeaks8$propCellsUsed[1])
    ){
      possiblePeaks9 <- possiblePeaks8 %>% dplyr::mutate(
        data=flowName@description[["GUID"]]
      )
      
      appendData <- rbind(
        appendData,
        possiblePeaks9
      )
    }else{
      flaggedData <- rbind(
        flaggedData,
        flowName@description[["GUID"]]
      )
      
    }
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
  }
  
  if(!purrr::is_empty(appendData) & !purrr::is_empty(flaggedData)){
    returnedList <- list(flaggedData, appendData)
  }else if(!purrr::is_empty(appendData) & purrr::is_empty(flaggedData)){
    returnedList <- list(appendData) 
  }else{
    returnedList <- list(flaggedData) 
  }
  
  
  return(returnedList)
  
}

## peakAlgorithm4
.peakAlgorithm4 = function(
  flowDir,
  flaggedData_, 
  xVariable, 
  appendData, 
  usedCellsThreshold,
  maxDoubletHeight
){
  
  ##Removing NOTE 'no visible binding for global variable'
  y<-possiblePairY<-keep<-truePeak<-NULL
  
  flowNameDs <- flaggedData_$data
  flaggedData <- c()
  
  logFlow <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  
  algorithmNum <- 4
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    
    flowData <- .smoothData( flowName, xVariable, 4)
    
    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)
    
    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )
    
    localPeaks <- scorepeak::detect_localmaxima(flowData$y, 3)
    possiblePeaks <- flowData[localPeaks, ]
    
    possiblePeaks1 <- possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/7),
    ]
    
    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/3.5),
    ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/9.5),]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > max(possiblePeaks$y)/3.5),
      ]
    }
    possiblePeaks3 <- .findTruePeaks(possiblePeaks2, 40, xVarMax)
    
    possiblePeaks4 <- .findPairs(possiblePeaks3, possiblePeaks, 1.7, 2.3)
    
    possiblePeaks5 = possiblePeaks4 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        possiblePairY > quantile(flowData$y)[3]+10
      )
    possiblePeaks6 <- possiblePeaks5[
      !duplicated(possiblePeaks5$possiblePairY),]
    
    
    if(nrow(possiblePeaks6) > 1){
      possiblePeaks7 <- .findTruePeaks(possiblePeaks6, 40, xVarMax)
    }else{
      possiblePeaks7 <- possiblePeaks6
    }
    
    ##Checking peak sizes to see if it's a true peak
    if(nrow(possiblePeaks7) > 1){
      tempDs<-possiblePeaks7
      tempDs2<-c()
      for(i in seq_len(nrow(tempDs))){
        popInQuestion<-tempDs[i,]
        if(is.na(maxDoubletHeight)){
          maxDoubletHeight_ <- round((popInQuestion$y)/3.5,2)
        }else{
          maxDoubletHeight_<-maxDoubletHeight
        }
        if(popInQuestion$possiblePairY < maxDoubletHeight_){
          popInQuestion$truePeak = FALSE
        }else{
          popInQuestion$truePeak = TRUE
        }
        tempDs2<-rbind(
          tempDs2,
          popInQuestion
        )
      }
      possiblePeaks7<-tempDs2 %>% 
        dplyr::filter(truePeak == TRUE) %>% 
        dplyr::select(-truePeak)
      rm(tempDs, tempDs2)  
    }
    
    ## Finding Doublets
    rangeLength <- nchar(format(xVarMax, scientific=FALSE))
    multiplier <- 10^(rangeLength-3)
    possiblePeaks8 <- .doubletCheck(
      possiblePeaks7,
      possiblePeaks,
      10*multiplier,
      15*multiplier
    )
    
    ##Checking if G2+G2 doublet is a true doublet and not a G2 peak for
    ##Another subpopulation
    tempDs<-possiblePeaks8
    if(nrow(tempDs) > 1 & !is.na(tempDs$g1G2DoubletCount[1])){
      pop1<-tempDs[1,]
      for(i in 2:nrow(tempDs)){
        popInQuestion = tempDs[i,]
        if(pop1$g2G2Doublet == popInQuestion$possiblePairX & 
           popInQuestion$possiblePairY < maxDoubletHeight_){
          tempDs <- tempDs[-i,]
        }else{
          tempDs$g1G2Doublet<-NA
          tempDs$g2G2Doublet<-NA
          tempDs$g1G2DoubletCount<-NA
          tempDs$g2G2DoubletCount<-NA
        }
      }
    }
    possiblePeaks8<-tempDs
    rm(tempDs) 
    
    if(nrow(possiblePeaks8) != 0){
      if(nrow(possiblePeaks8) == 1){
        peakRow <- data.frame(
          peaks=c(
            possiblePeaks8$x,
            possiblePeaks8$possiblePairX,
            possiblePeaks8$g1G2Doublet,
            possiblePeaks8$g2G2Doublet
          )
        )
      }else{
        peakRow <- data.frame(
          peaks=c(possiblePeaks8$x, possiblePeaks8$possiblePairX)
        )
      }
      
      peakRow <- peakRow %>%
        dplyr::distinct() %>% tidyr::drop_na()
      
      cellsUsed <- c()
      smoothedData <- .smoothData(flowName, xVariable, 15)
      
      for(j in seq_len(nrow(peakRow))){
        xPeak <- which(smoothedData$x == peakRow$peaks[j])
        if(length(xPeak)>1){
          xPeak <-xPeak[1]
        }
        
        if(xPeak+10 > nrow(smoothedData) & xPeak-10 < 1){
          peakRowSmoothedRange <- smoothedData[
            seq(1, nrow(smoothedData), 1), ]
        }else if(xPeak+10 > nrow(smoothedData) & xPeak-10 >= 1){
          peakRowSmoothedRange <- smoothedData[
            seq(xPeak-10, nrow(smoothedData), 1),
          ]
        }else if(xPeak+10 <= nrow(smoothedData) & xPeak-10 < 1){
          peakRowSmoothedRange <- smoothedData[seq(1, xPeak+10, 1), ]
        }else{
          peakRowSmoothedRange <- smoothedData[
            seq(xPeak-10, xPeak+10, 1), ]
        }
        
        xPeakSmoothed <- which(
          max(peakRowSmoothedRange$y) == smoothedData$y
        )
        if(length(xPeakSmoothed) > 1){
          xPeakSmoothed <- xPeakSmoothed[1]
        }
        smoothedData$keep <- FALSE
        smoothedData$keep[xPeakSmoothed] <- TRUE
        
        #Forwards
        for(i in xPeakSmoothed:nrow(smoothedData)){
          if(i == nrow(smoothedData)){
            smoothedData$keep[nrow(smoothedData)] <- TRUE
            break
          }else{
            if(smoothedData$y[i+1] > smoothedData$y[i]){
              smoothedData$keep[i+1] <- TRUE
              break
            }
          }
        }
        
        #Backwards
        for(i in xPeakSmoothed:1){
          if(i == 1){
            smoothedData$keep[1] <- TRUE
            break
          }else{
            if(smoothedData$y[i-1] > smoothedData$y[i]){
              smoothedData$keep[i-1] <- TRUE
              break
            }
          }
        }
        
        intPeak <- smoothedData %>% dplyr::filter(
          keep == TRUE
        )
        
        xEpsilonRight <- which(flowData$x == intPeak$x[3])
        xEpsilonLeft <- which(flowData$x == intPeak$x[1])
        peakCount <- flowData$y[c(xEpsilonLeft:xEpsilonRight)]
        peakCount01 <- sum(peakCount)
        peakCount01[j] <- peakCount01
        cellsUsed <- c(
          cellsUsed,
          peakCount01[j]
        )
        
      }
      
      cellsUsed01 <- sum(cellsUsed)
      totalCellCount <- sum(smoothedData$y)
      propCellsUsed <- round((cellsUsed01/totalCellCount)*100, 2)
      possiblePeaks8$propCellsUsed <- propCellsUsed
      
    }else{
      possiblePeaks8 <- data.frame(
        x=NA,
        y=NA,
        cluster=NA,
        distToNext=NA,
        LL=NA,
        UL=NA,
        possiblePairX=NA,
        possiblePairY=NA,
        g3LL=NA,
        g3UL=NA,
        g4LL=NA,
        g4UL=NA,
        g1G2Doublet=NA,
        g1G2DoubletCount=NA,
        g2G2Doublet=NA,
        g2G2DoubletCount=NA,
        data=flowName@description[["GUID"]],
        propCellsUsed=NA
      )
    }
    
    if(
      possiblePeaks8$propCellsUsed[1] >= 75 &
      !is.na(possiblePeaks8$propCellsUsed[1])
    ){
      possiblePeaks9 <- possiblePeaks8 %>% dplyr::mutate(
        data=flowName@description[["GUID"]]
      )
      
      appendData <- rbind(
        appendData,
        possiblePeaks9
      )
      
    }else{
      possiblePeaks9 <- possiblePeaks8 %>% dplyr::mutate(
        data=flowName@description[["GUID"]],
        investigate=1
      )
      
      flaggedData <- rbind(
        flaggedData,
        possiblePeaks9
      )
      
    }
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
  }
  
  if(!purrr::is_empty(appendData) & !purrr::is_empty(flaggedData)){
    returnedList <- list(flaggedData, appendData)
  }else if(!purrr::is_empty(appendData) & purrr::is_empty(flaggedData)){
    returnedList <- list(appendData) 
  }else{
    returnedList <- list(flaggedData) 
  }
  
  return(returnedList)
  
}

## peakAlgorithm5
.peakAlgorithm5 = function(
  flowDir,
  flaggedData_, 
  xVariable, 
  appendData
){
  ##Removing NOTE 'no visible binding for global variable'
  x<-y<-possiblePairY<-keep<-NULL
  
  naDs <- flaggedData_ %>% dplyr::filter(is.na(x))
  flowNameDs <- unique(naDs$data)
  
  logFlow <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  
  algorithmNum <- 5
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    
    flowData <- .smoothData( flowName, xVariable, 5)
    
    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)
    
    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )
    
    localPeaks <- scorepeak::detect_localmaxima(flowData$y, 3)
    possiblePeaks <- flowData[localPeaks, ]
    
    possiblePeaks1 <- possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/9),
    ]
    
    possiblePeaks2 = possiblePeaks[
      which(possiblePeaks$y > max(possiblePeaks$y)/4),
    ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/9.5),]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > max(possiblePeaks$y)/4.5),
      ]
    }
    possiblePeaks3 <- .findTruePeaks(possiblePeaks2, 30, xVarMax)
    
    possiblePeaks4 <- possiblePeaks3[
      order(possiblePeaks3$x, decreasing=FALSE),
    ]
    
    peaksFix <- possiblePeaks4
    firstPeak <- possiblePeaks4[1, seq_len(2)]
    
    if(nrow(possiblePeaks4) == 1){
      
      initialPeak <- firstPeak
      #the peak is G1 missing G2
      g2PeakRadiusUL <- initialPeak$x*2.2 + 1
      g2PeakRadiusLL <- initialPeak$x*1.75 - 1
      g2ToTestDs2 <- which(
        g2PeakRadiusLL < possiblePeaks$x &
          possiblePeaks$x < g2PeakRadiusUL &
          possiblePeaks$x < quantile(flowData$x)[4]
      )
      g2ToTestDs3 <- possiblePeaks[g2ToTestDs2, ]
      if(nrow(g2ToTestDs3) == 0){
        g2ToTestDs3 <- data.frame(
          x=0,
          y=0
        )
      }
      #the peak is G2 missing G1
      g1PeakRadiusLL <- initialPeak$x/2.2 + 1
      g1PeakRadiusUL <- initialPeak$x/1.75 - 1
      g1ToTestDs2 <- which(
        g1PeakRadiusLL < possiblePeaks$x & 
          possiblePeaks$x < g1PeakRadiusUL
      )
      g1ToTestDs3 <- possiblePeaks[g1ToTestDs2, ]
      if(nrow(g1ToTestDs3) == 0){
        g1ToTestDs3 <- data.frame(
          x=0,
          y=0
        )
      }
      
      if(g2ToTestDs3$y[1] < g1ToTestDs3$y[1]){
        peak2 <- g1ToTestDs3[1, ]
        peak1 <- initialPeak
        peaksFix <- rbind(
          peak1,
          peak2
        )
      }else{
        peak2 <- g2ToTestDs3[1, ]
        peak1 <- initialPeak
        peaksFix <- rbind(
          peak1,
          peak2
        )
      }
      
    }
    peaksFix<-peaksFix[order(peaksFix$x, decreasing=FALSE), ]
    firstPeak<-peaksFix[1,]
    if(firstPeak$x > xVarMax/3.1 | firstPeak$x < xVarMax/6){
      
      initialPeak <- firstPeak
      otherPeak <- peaksFix[-1, seq_len(2)]
      
      ##the peak is G2 missing G1
      g1PeakRadiusLL <- initialPeak$x/2.3 - 1
      g1PeakRadiusUL <- initialPeak$x/1.7 + 1
      g1ToTestDs2 <- which(
        g1PeakRadiusLL <= possiblePeaks$x & 
          possiblePeaks$x <= g1PeakRadiusUL
      )
      g1ToTestDs3 <- possiblePeaks[g1ToTestDs2, ]
      
      if(nrow(g1ToTestDs3) == 0){
        g1ToTestDs2 <- which(
          g1PeakRadiusLL <= flowData$x & flowData$x <= g1PeakRadiusUL
        )
        g1ToTestDs3 <- flowData[g1ToTestDs2, ]
      }
      g1ToTestDs3 <- g1ToTestDs3[order(g1ToTestDs3$y, decreasing=TRUE), ]
      g1ToTestDs3 <- g1ToTestDs3[
        which(g1ToTestDs3$y > quantile(flowData$y)[2]), ]
      peak2 <- g1ToTestDs3[1, ]
      peak1 <- initialPeak
      peaksFix <- rbind(
        peak1,
        peak2,
        otherPeak
      )
    }
    possiblePeaks5 <- .findPairs(peaksFix, peaksFix, 1.4, 2.5)
    
    possiblePeaks5 <- possiblePeaks5 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        possiblePairY > quantile(flowData$y)[3]
      )
    possiblePeaks6 <- possiblePeaks5[
      !duplicated(possiblePeaks5$possiblePairY),]
    
    if(nrow(possiblePeaks6) > 1){
      possiblePeaks7 <- .findTruePeaks(possiblePeaks6, 30, xVarMax)
    }else{
      possiblePeaks7 <- possiblePeaks6 %>% dplyr::mutate(
        cluster=1,
        distToNext=0
      )
    }
    
    if(nrow(possiblePeaks7) == 0){
      peaksFix<-peaksFix %>% tidyr::drop_na()
      maxPeaksFix <- peaksFix[which(max(peaksFix$y) == peaksFix$y), ]
      possiblePeaks7 <- data.frame(
        x=maxPeaksFix$x,
        y=maxPeaksFix$y,
        cluster=1,
        distToNext=0,
        LL=NA,
        UL=NA,
        possiblePairX=0,
        possiblePairY=0
      )
      
    }
    
    ## Finding Doublets
    rangeLength <- nchar(format(xVarMax, scientific=FALSE))
    multiplier <- 10^(rangeLength-3)
    possiblePeaks8 <- .doubletCheck(
      possiblePeaks7,
      possiblePeaks,
      10*multiplier,
      15*multiplier
    )
    
    possiblePeaks9 <- possiblePeaks8 %>% dplyr::mutate(
      data=flowName@description[["GUID"]],
      investigate=1
    )
    
    appendData <- rbind(
      appendData,
      possiblePeaks9
    )
    
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
  }
  
  returnedList <- list(appendData)
  
  return(returnedList)
}

## smoothData
.smoothData = function(flowDs, xVariable, smoothLevel){
  
  ##Get counts and breaks from histogram
  histData <- hist(flowDs@exprs[, xVariable], breaks=256, plot=FALSE)
  
  ##Data that will be smoothed
  data <- histData$counts
  ##Apply smoothing to the counts with 'smoothLevel'
  smoothedDs <- zoo::rollmean(data, k=smoothLevel, fill=0)
  
  ##Create data frame with smoothed data
  ds <- data.frame(
    x=histData$breaks,
    y=c(0, smoothedDs)
  )
  return(ds)
}

## findTruePeaks
.findTruePeaks = function(ds, clusterDist, maxXValue){
  ##Removing NOTE 'no visible binding for global variable'
  cluster<-y<-NULL
  ##Finding initial distance between all identified peaks
  tempDs <- ds
  tempDs$cluster <- 1
  tempDs$distToNext <- c(diff(tempDs$x), 0)
  
  ##Creating an distance requirement between two peaks
  maxDist <- maxXValue/clusterDist
  clusterNum <- 1
  ##Going through each peak and checking if their distance to the next peak is
  ##greater than the distance requirement 'maxDist'. If not, the peaks will be
  ##grouped together in the same cluster.
  for (i in seq_len(nrow(tempDs))) {
    if(tempDs$distToNext[i] > maxDist){
      tempDs$cluster[i] <- clusterNum
      clusterNum <- clusterNum + 1
    } else {
      tempDs$cluster[i] <- clusterNum
    }
  }
  ##If clusters have more than one peak identified, the tallest peak will be
  ##selected as the peak in that cluster
  tempDs <- tempDs %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice(which.max(y))
  
  return(tempDs)
  
}

## findPairs
.findPairs = function(ds, peaks, LL, UL){
  ##Removing NOTE 'no visible binding for global variable'
  x<-NULL
  ##creating upper and lowed bounds that are read into the function
  findingPairsDs <- ds %>% dplyr::mutate(
    LL = x*LL,
    UL = x*UL
  )
  ##Finding the peaks that are in the in (LL, UL)
  ##These are the possible peaks that will be considered for their G1/G2 
  ##pairing if there is more than one peak identified, we pick the tallest
  ##of those peaks
  findingPairsDs$possiblePairX <- NA
  findingPairsDs$possiblePairY <- NA
  for( i in seq_len(nrow(findingPairsDs))){
    possible_ <- which(
      peaks$x >= findingPairsDs$LL[i] &
        peaks$x <= findingPairsDs$UL[i]
    )
    
    if(length(possible_) > 1){
      possible_ <- which(
        peaks$y == max(peaks[possible_,]$y) &
          peaks$x >= findingPairsDs$LL[i]
      )
      possible_<-possible_[1]
    }
    
    if(!purrr::is_empty(possible_)){
      if(possible_ != i) {
        maxPossible_ <- peaks[possible_, ]
        maxPossibleRows <- maxPossible_[
          order(maxPossible_$y, decreasing = TRUE),][1,]
        
        findingPairsDs$possiblePairX[i] <- maxPossibleRows$x
        findingPairsDs$possiblePairY[i] <- maxPossibleRows$y
      }else{
        possible_ <- NA
      }
    }
  }
  
  return(findingPairsDs)
}

## doubletCheck
.doubletCheck = function(doubletCheckDs, peaks, g1G2Range, g2G2Range){
  ##Removing NOTE 'no visible binding for global variable'
  x<-possiblePairX<-NULL
  ##creating lower bounds and upper bounds for peaks that
  ##could be classified as doublets
  doubletCheckDs2 <- doubletCheckDs %>% dplyr::mutate(
    g3LL=x + possiblePairX - g1G2Range,
    g3UL=x + possiblePairX + g1G2Range,
    g4LL=possiblePairX + possiblePairX - g2G2Range,
    g4UL=possiblePairX + possiblePairX + g2G2Range
  )
  
  ##finding the peaks that are in the G1+G2 range
  doubletCheckDs2$g1G2Doublet <- NA
  doubletCheckDs2$g1G2DoubletCount <- NA
  for( i in seq_len(nrow(doubletCheckDs2))){
    possible_ <- which(
      peaks$x >= doubletCheckDs2$g3LL[i] &
        peaks$x <= doubletCheckDs2$g3UL[i]
    )
    
    if(length(possible_) > 1){
      possible_ <- which(
        peaks$y == max(peaks[possible_,]$y) &
          peaks$x >= doubletCheckDs2$g3LL[i]   
      )
      possible_ <-possible_[1]
    }
    
    if(!purrr::is_empty(possible_)){
      if(possible_ != i){
        ##If there are more than one peak identified, we pick the 
        ##tallest of those peaks
        maxPossible_ <- peaks[possible_, ]
        maxPossibleRows <- maxPossible_[
          order(maxPossible_$y, decreasing = TRUE),
        ][1, ]
        
        doubletCheckDs2$g1G2Doublet[i] <- maxPossibleRows$x
        doubletCheckDs2$g1G2DoubletCount[i] <- maxPossibleRows$y
      }else{
        possible_ <- NA
      }
    }
  }
  
  ##similarly for G2+G2
  ##finding the peaks that are in the G2+G2 range
  doubletCheckDs2$g2G2Doublet <- NA
  doubletCheckDs2$g2G2DoubletCount <- NA
  for( i in seq_len(nrow(doubletCheckDs2))){
    possible_ <- which(
      peaks$x >= doubletCheckDs2$g4LL[i] &
        peaks$x <= doubletCheckDs2$g4UL[i]
    )
    if(length(possible_) > 1){
      possible_ <- which(
        peaks$y == max(peaks[possible_,]$y) &
          peaks$x >= doubletCheckDs2$g4LL[i]  
      )
      possible_ <-possible_[1]
    }
    
    
    if(!purrr::is_empty(possible_)){
      if(possible_ != i) {
        maxPossible_ <- peaks[possible_,]
        maxPossibleRows <- maxPossible_[
          order(maxPossible_$y, decreasing = TRUE), ][1, ]
        
        doubletCheckDs2$g2G2Doublet[i] <- maxPossibleRows$x
        doubletCheckDs2$g2G2DoubletCount[i] <- maxPossibleRows$y
      }else{
        possible_ <- NA
      }
    }
  }
  
  return(doubletCheckDs2)
}


## updatedMeans
.updatedMeans = function(ds, flowDir, xVariable){
  ##Removing NOTE 'no visible binding for global variable'
  x<-y<-possiblePairX<-possiblePairY<-g3LL<-g3UL<-NULL
  g4LL<-g4UL<-g1G2Doublet<-g1G2DoubletCount<-g2G2Doublet<-NULL
  g2G2DoubletCount<-possiblePairXOld<-possiblePairYOld<-NULL
  
  flowNameDs <- ds$data
  singlePopUpdated <- c()
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    flowData <- .smoothData( flowName, xVariable, 5)
    singlePop01 <- ds %>% dplyr::filter(
      data == flowNameDs[k]
    ) %>% dplyr::rename(
      xOld=x,
      yOld=y,
      possiblePairXOld=possiblePairX,
      possiblePairYOld=possiblePairY
    )
    
    midPoint <- flowData[
      which(
        abs(
          flowData$x-mean(c(singlePop01$xOld,
                            singlePop01$possiblePairXOld))
        ) == min(
          abs(
            flowData$x-mean(c(singlePop01$xOld,
                              singlePop01$possiblePairXOld))
          )
        )
      ),
    ]
    
    if(nrow(midPoint)>1){
      midPoint <- midPoint[nrow(midPoint), ]
    }
    
    xPeak1 <- which(flowData$x == singlePop01$xOld)
    if(length(xPeak1)>1){
      xPeak1 <- xPeak1[1]
    }
    
    if(xPeak1+10 > midPoint$x & xPeak1-10 < 1){
      peakRowSmoothedRange <- flowData[seq(1, midPoint$x, 1), ]
    }else if(xPeak1+10 > midPoint$x & xPeak1-10 >= 1){
      peakRowSmoothedRange <- flowData[
        seq(xPeak1-10, midPoint$x, 1),
      ]
    }else if(xPeak1+10 <= midPoint$x & xPeak1-10 < 1){
      peakRowSmoothedRange <- flowData[seq(1, xPeak1+10, 1), ]
    }else{ 
      peakRowSmoothedRange <- flowData[seq(xPeak1-10, xPeak1+10, 1), ]
    }
    
    xPeak1Smoothed <- which(max(peakRowSmoothedRange$y) == flowData$y)
    if(length(xPeak1Smoothed) > 1){
      xPeak1Smoothed <- xPeak1Smoothed[1]
    }
    
    if(singlePop01$possiblePairXOld != 0){
      xPeak2 <- which(flowData$x == singlePop01$possiblePairXOld)
      if(length(xPeak2)>1){
        xPeak2 <- xPeak2[1]
      }
      midPointX = which(midPoint$x == flowData$x)
      if(xPeak2+10 > nrow(flowData) & xPeak2-10 < midPointX){
        peakRowSmoothedRange <- flowData[
          seq(midPointX, nrow(flowData), 1), ]
      }else if(xPeak2+10 > nrow(flowData) & xPeak2-10 >= midPointX){
        peakRowSmoothedRange <- flowData[
          seq(xPeak2-10, nrow(flowData), 1),
        ]
      }else if(xPeak2+10 <= nrow(flowData) & xPeak2-10 < midPointX){
        peakRowSmoothedRange <- flowData[seq(midPointX, xPeak2+10, 1), ]
      }else{
        peakRowSmoothedRange <- flowData[seq(xPeak2-10, xPeak2+10, 1), ]
      }
      
      xPeak2Smoothed <- which(max(peakRowSmoothedRange$y) == flowData$y)
      if(length(xPeak2Smoothed) > 1){
        xPeak2Smoothed <- xPeak2Smoothed[
          which(
            xPeak2Smoothed > midPointX
          )
        ]
        xPeak2Smoothed <- xPeak2Smoothed[1]
      }
      
      singlePop02 <- singlePop01 %>% dplyr::mutate(
        x=flowData[xPeak1Smoothed, ]$x,
        y=flowData[xPeak1Smoothed, ]$y,
        possiblePairX=flowData[xPeak2Smoothed, ]$x,
        possiblePairY=flowData[xPeak2Smoothed, ]$y
      ) %>% dplyr::select(
        data,
        x,
        y,
        possiblePairX,
        possiblePairY,
        g3LL,
        g3UL,
        g4LL,
        g4UL,
        g1G2Doublet,
        g1G2DoubletCount,
        g2G2Doublet,
        g2G2DoubletCount
      )
    }else{
      singlePop02 <- singlePop01 %>% dplyr::mutate(
        x=flowData[xPeak1Smoothed, ]$x,
        y=flowData[xPeak1Smoothed, ]$y,
        possiblePairX=possiblePairXOld,
        possiblePairY=possiblePairYOld
      ) %>% dplyr::select(
        data,
        x,
        y,
        possiblePairX,
        possiblePairY,
        g3LL,
        g3UL,
        g4LL,
        g4UL,
        g1G2Doublet,
        g1G2DoubletCount,
        g2G2Doublet,
        g2G2DoubletCount
      )
    }
    
    
    singlePopUpdated <- rbind(
      singlePopUpdated,
      singlePop02
    )
    
  }
  returnedList <- list(singlePopUpdated)
}


## outputData
.outputData = function(
  flowDir,
  singleDs,
  finishedDs,
  investigateDs,
  peak1Check,
  xVariable,
  doubletFlag,
  saveGraph = TRUE
){
  ##Removing NOTE 'no visible binding for global variable'
  x<-y<-.<-possiblePairX<-possiblePairY<-G1<-G1Count<-G2<-G2Count<-id<-NULL
  g1G2Doublet<-g1G2DoubletCount<-g2G2Doublet<-g2G2DoubletCount<-NULL
  residual<-residualDoublet<-Success<-Algorithm<-Data<-`doublet G1+G2`<-NULL
  `doublet G1+G2 count`<-`doublet G2+G2`<- `doublet G2+G2 count`<-NULL
  residual3Pop<-residual2Pop<-residualMultiple<-NULL
  ##If the algorithm classified each sub population in the first algorithm
  if(
    purrr::is_empty(finishedDs) &
    purrr::is_empty(investigateDs) &
    !purrr::is_empty(singleDs)
  ){
    ##Formatting the single pop data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
      dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    finalData=finalPart1 %>% dplyr::select(
      c(
        "data",
        "G1",
        "G2",
        "g1G2Doublet",
        "g2G2Doublet",
        "G1Count",
        "G2Count",
        "g1G2DoubletCount",
        "g2G2DoubletCount",
        "investigate"
      )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    if(!purrr::is_empty(peak1Check)){
      finalData$investigate[finalData$data %in% peak1Check$data] = 1
    }
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
      dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
      dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
      dplyr::mutate(id=rowid(data)) %>%
      tidyr::pivot_wider(
        names_from=id,
        values_from=c(G1, G1Count, G2, G2Count)
      )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    
    finalData3=finalData3 %>%
      dplyr::rename(
        `doublet G1+G2`=g1G2Doublet,
        `doublet G1+G2 count`=g1G2DoubletCount,
        `doublet G2+G2`=g2G2Doublet,
        `doublet G2+G2 count`=g2G2DoubletCount
      )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
      doubletRSE<-.popConfidenceDoublet(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      doubletRSE<-finalData3 %>% dplyr::mutate(
        residualDoublet=NA
      )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
      if(TRUE %in% grepl("_3", names(finalData3))){
        ds1<-finalData3 %>% dplyr::filter(is.na(G1_3))    
      }else{
        ds1<-finalData3
      }
      
      if(TRUE %in% !is.na(ds1$G1_2)){
        twoPopRSE<-.popConfidence2Pop(
          flowDir, ds=ds1, xVariable, saveGraph
        )
      }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
          residual2Pop=NA
        )
      }
    }else{
      twoPopRSE<-finalData3 %>% dplyr::mutate(
        residual2Pop=NA
      )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
      threePopRSE<-.popConfidence3Pop(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      threePopRSE<-finalData3 %>% dplyr::mutate(
        residual3Pop=NA
      )
    }
    
    finalData4<-sqldf::sqldf(
      "select ds.*,
                ds2.residual,
                ds3.residualDoublet,
                ds4.residual2Pop,
                ds5.residual3Pop
                from finalData3 ds
                left join initialRSE ds2
                on ds.data = ds2.data
                left join doubletRSE ds3
                on ds.data = ds3.data
                left join twoPopRSE ds4
                on ds.data = ds4.data
                left join threePopRSE ds5
                on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
      residualMultiple=ifelse(
        is.na(residual3Pop), residual2Pop, residual3Pop)
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
      finalRSE=do.call(
        pmin, 
        c(
          subset(., select=c(
            residual,
            residualDoublet,
            residualMultiple)
          ), na.rm=TRUE
        )
      )
    ) %>% dplyr::rename(
      singleRSE=residual,
      doubletRSE=residualDoublet,
      multipleRSE=residualMultiple
    ) %>% dplyr::select(
      -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
      if(
        !is.na(finalData5[i,]$multipleRSE) &
        finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
      ){
        finalData5[i,]$investigate = 1
      }
      
    }
    
    ##Flagging to investigate if the RSE is considered an outlier       
    outlierRSE = as.numeric(
      quantile(finalData5$finalRSE)[4]+
        1.5*stats::IQR(finalData5$finalRSE)
    )
    
    for(i in seq_len(nrow(finalData5))){
      if( finalData5$finalRSE[i] >= outlierRSE ){
        finalData5[i,]$investigate = 1
      }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
      "Sample" = "data"
    )
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
      finalData6=finalData6 %>%
        dplyr::select(
          -c(
            `doublet G1+G2`,
            `doublet G1+G2 count`,
            `doublet G2+G2`,
            `doublet G2+G2 count`
          )
        )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings=FALSE)
    experimentName <- basename(dirname(getwd()))
    experimentName <- sub(" ", "_", experimentName)   # replaces spaces with an underscore
    write.csv(
      finalData6,
      paste0(
        file.path(dirname(getwd()), subDir),
        "/", experimentName, "_ploidyPeaksOutput.csv"
      ),
      row.names = FALSE
    )
    
  }else if(
    purrr::is_empty(investigateDs) &
    !purrr::is_empty(finishedDs) &
    !purrr::is_empty(singleDs)
  ){
    ##If the algorithm classified each sub population before prior the 4th
    ##algorithm
    
    ##Formatting the diploid data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
      dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    ##Formatting the data from the other peak algorithms
    finalPart2=finishedDs %>% data.frame() %>%
      dplyr::select(
        -c(
          "g3LL",
          "g3UL",
          "g4LL",
          "g4UL",
          "propCellsUsed",
          "cluster",
          "distToNext",
          "LL",
          "UL"
        )
      ) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    ##Merging the two datasets
    finalData=rbind(
      finalPart1,
      finalPart2
    ) %>% dplyr::select(
      c(
        "data",
        "G1",
        "G2",
        "g1G2Doublet",
        "g2G2Doublet",
        "G1Count",
        "G2Count",
        "g1G2DoubletCount",
        "g2G2DoubletCount",
        "investigate"
      )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    if(!purrr::is_empty(peak1Check)){
      finalData$investigate[finalData$data %in% peak1Check$data] = 1
    }
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
      dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
      dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
      dplyr::mutate(id=rowid(data)) %>%
      tidyr::pivot_wider(
        names_from=id,
        values_from=c(G1, G1Count, G2, G2Count))
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    finalData3=finalData3 %>%
      dplyr::rename(
        `doublet G1+G2`=g1G2Doublet,
        `doublet G1+G2 count`=g1G2DoubletCount,
        `doublet G2+G2`=g2G2Doublet,
        `doublet G2+G2 count`=g2G2DoubletCount
      )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
      doubletRSE<-.popConfidenceDoublet(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      doubletRSE<-finalData3 %>% dplyr::mutate(
        residualDoublet=NA
      )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
      if(TRUE %in% grepl("_3", names(finalData3))){
        ds1<-finalData3 %>% dplyr::filter(is.na(G1_3))    
      }else{
        ds1<-finalData3
      }
      
      if(TRUE %in% !is.na(ds1$G1_2)){
        twoPopRSE<-.popConfidence2Pop(
          flowDir, ds=ds1, xVariable, saveGraph
        )
      }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
          residual2Pop=NA
        )
      }
    }else{
      twoPopRSE<-finalData3 %>% dplyr::mutate(
        residual2Pop=NA
      )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
      threePopRSE<-.popConfidence3Pop(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      threePopRSE<-finalData3 %>% dplyr::mutate(
        residual3Pop=NA
      )
    }
    
    finalData4<-sqldf::sqldf(
      "select ds.*,
                ds2.residual,
                ds3.residualDoublet,
                ds4.residual2Pop,
                ds5.residual3Pop
                from finalData3 ds
                left join initialRSE ds2
                on ds.data = ds2.data
                left join doubletRSE ds3
                on ds.data = ds3.data
                left join twoPopRSE ds4
                on ds.data = ds4.data
                left join threePopRSE ds5
                on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
      residualMultiple=ifelse(
        is.na(residual3Pop),
        residual2Pop, residual3Pop
      )
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
      finalRSE=do.call(
        pmin, 
        c(
          subset(., 
                 select=c(
                   residual, 
                   residualDoublet, 
                   residualMultiple)
          ), na.rm=TRUE
        )
      )
    ) %>% dplyr::rename(
      singleRSE=residual,
      doubletRSE=residualDoublet,
      multipleRSE=residualMultiple
    ) %>% dplyr::select(
      -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
      if(
        !is.na(finalData5[i,]$multipleRSE) &
        finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
      ){
        finalData5[i,]$investigate = 1
      }
      
    }
    
    ##Flagging to investigate if the RSE is considered an outlier       
    outlierRSE = as.numeric(
      quantile(finalData5$finalRSE)[4]+
        1.5*stats::IQR(finalData5$finalRSE)
    )
    
    for(i in seq_len(nrow(finalData5))){
      if( finalData5$finalRSE[i] >= outlierRSE ){
        finalData5[i,]$investigate = 1
      }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
      "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
      finalData6=finalData6 %>%
        dplyr::select(
          -c(
            `doublet G1+G2`,
            `doublet G1+G2 count`,
            `doublet G2+G2`,
            `doublet G2+G2 count`
          )
        )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)
    experimentName <- basename(dirname(getwd()))
    experimentName <- sub(" ", "_", experimentName)   # replaces spaces with an underscore
    write.csv(
      finalData6,
      paste0(file.path(dirname(getwd()), subDir), "/", experimentName, "_ploidyPeaksOutput.csv"
      ),
      row.names = FALSE
    )
    
  }else if(
    purrr::is_empty(singleDs) &
    !purrr::is_empty(finishedDs) &
    !purrr::is_empty(investigateDs)
  ){
    ##If the algorithm classified each sub population before prior the 4th
    ##algorithm
    
    ##Formatting the data from the other peak algorithms
    finalPart1<-finishedDs %>% data.frame() %>%
      dplyr::select(
        -c(
          "g3LL",
          "g3UL",
          "g4LL",
          "g4UL",
          "propCellsUsed",
          "cluster",
          "distToNext",
          "LL",
          "UL"
        )
      ) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    finalPart2<-investigateDs %>% data.frame() %>%
      dplyr::select(
        -c(
          "g3LL",
          "g3UL",
          "g4LL",
          "g4UL",
          "cluster",
          "distToNext",
          "LL",
          "UL"
        )
      ) %>%
      dplyr::mutate(
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    ##Merging the two datasets
    finalData<-rbind(
      finalPart1,
      finalPart2
    ) %>% dplyr::select(
      c(
        "data",
        "G1",
        "G2",
        "g1G2Doublet",
        "g2G2Doublet",
        "G1Count",
        "G2Count",
        "g1G2DoubletCount",
        "g2G2DoubletCount",
        "investigate"
      )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    if(!purrr::is_empty(peak1Check)){
      finalData$investigate[finalData$data %in% peak1Check$data] = 1
    }
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
      dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
      dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
      dplyr::mutate(id=rowid(data)) %>%
      tidyr::pivot_wider(
        names_from=id, values_from=c(G1, G1Count, G2, G2Count)
      )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    finalData3=finalData3 %>%
      dplyr::rename(
        `doublet G1+G2`=g1G2Doublet,
        `doublet G1+G2 count`=g1G2DoubletCount,
        `doublet G2+G2`=g2G2Doublet,
        `doublet G2+G2 count`=g2G2DoubletCount
      )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
      doubletRSE<-.popConfidenceDoublet(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      doubletRSE<-finalData3 %>% dplyr::mutate(
        residualDoublet=NA
      )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
      if(TRUE %in% grepl("_3", names(finalData3))){
        ds1<-finalData3 %>% dplyr::filter(is.na(G1_3))    
      }else{
        ds1<-finalData3
      }
      
      if(TRUE %in% !is.na(ds1$G1_2)){
        twoPopRSE<-.popConfidence2Pop(
          flowDir, ds=ds1, xVariable, saveGraph
        )
      }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
          residual2Pop=NA
        )
      }
    }else{
      twoPopRSE<-finalData3 %>% dplyr::mutate(
        residual2Pop=NA
      )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
      threePopRSE<-.popConfidence3Pop(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      threePopRSE<-finalData3 %>% dplyr::mutate(
        residual3Pop=NA
      )
    }
    
    finalData4<-sqldf::sqldf(
      "select ds.*,
                ds2.residual,
                ds3.residualDoublet,
                ds4.residual2Pop,
                ds5.residual3Pop
                from finalData3 ds
                left join initialRSE ds2
                on ds.data = ds2.data
                left join doubletRSE ds3
                on ds.data = ds3.data
                left join twoPopRSE ds4
                on ds.data = ds4.data
                left join threePopRSE ds5
                on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
      residualMultiple=ifelse(
        is.na(residual3Pop), residual2Pop, residual3Pop)
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
      finalRSE=do.call(
        pmin, 
        c(
          subset(
            ., 
            select=c(
              residual, 
              residualDoublet, 
              residualMultiple
            )
          ),
          na.rm=TRUE
        )
      )
    ) %>% dplyr::rename(
      singleRSE=residual,
      doubletRSE=residualDoublet,
      multipleRSE=residualMultiple
    ) %>% dplyr::select(
      -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
      if(
        !is.na(finalData5[i,]$multipleRSE) &
        finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
      ){
        finalData5[i,]$investigate = 1
      }
      
    }
    
    ##Flagging to investigate if the RSE is considered an outlier       
    outlierRSE = as.numeric(
      quantile(finalData5$finalRSE)[4]+
        1.5*stats::IQR(finalData5$finalRSE)
    )
    
    for(i in seq_len(nrow(finalData5))){
      if( finalData5$finalRSE[i] >= outlierRSE ){
        finalData5[i,]$investigate = 1
      }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
      "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
      finalData6=finalData6 %>%
        dplyr::select(
          -c(
            `doublet G1+G2`,
            `doublet G1+G2 count`,
            `doublet G2+G2`,
            `doublet G2+G2 count`
          )
        )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)
    experimentName <- basename(dirname(getwd()))
    experimentName <- sub(" ", "_", experimentName)   # replaces spaces with an underscore
    write.csv(
      finalData6,
      paste0(file.path(dirname(getwd()), subDir), "/", experimentName, "_ploidyPeaksOutput.csv"
      ),
      row.names = FALSE
    )
    
  }else if(
    purrr::is_empty(investigateDs) &
    !purrr::is_empty(finishedDs) &
    purrr::is_empty(singleDs)
  ){
    
    ##Formatting the data from the other peak algorithms
    finalPart1<-finishedDs %>% data.frame() %>%
      dplyr::select(
        -c(
          "g3LL",
          "g3UL",
          "g4LL",
          "g4UL",
          "propCellsUsed",
          "cluster",
          "distToNext",
          "LL",
          "UL"
        )
      ) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    ##Removing columns
    finalData<-finalPart1 %>% dplyr::select(
      c(
        "data",
        "G1",
        "G2",
        "g1G2Doublet",
        "g2G2Doublet",
        "G1Count",
        "G2Count",
        "g1G2DoubletCount",
        "g2G2DoubletCount",
        "investigate"
      )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
      dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
      dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
      dplyr::mutate(id=rowid(data)) %>%
      tidyr::pivot_wider(
        names_from=id,
        values_from=c(G1, G1Count, G2, G2Count)
      )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    finalData3=finalData3 %>%
      dplyr::rename(
        `doublet G1+G2`=g1G2Doublet,
        `doublet G1+G2 count`=g1G2DoubletCount,
        `doublet G2+G2`=g2G2Doublet,
        `doublet G2+G2 count`=g2G2DoubletCount
      )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
      doubletRSE<-.popConfidenceDoublet(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      doubletRSE<-finalData3 %>% dplyr::mutate(
        residualDoublet=NA
      )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
      if(TRUE %in% grepl("_3", names(finalData3))){
        ds1<-finalData3 %>% dplyr::filter(is.na(G1_3))    
      }else{
        ds1<-finalData3
      }
      
      if(TRUE %in% !is.na(ds1$G1_2)){
        twoPopRSE<-.popConfidence2Pop(
          flowDir, ds=ds1, xVariable, saveGraph
        )
      }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
          residual2Pop=NA
        )
      }
    }else{
      twoPopRSE<-finalData3 %>% dplyr::mutate(
        residual2Pop=NA
      )
    }  
    if(TRUE %in% grepl("_3", names(finalData3))){
      threePopRSE<-.popConfidence3Pop(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      threePopRSE<-finalData3 %>% dplyr::mutate(
        residual3Pop=NA
      )
    }
    
    finalData4<-sqldf::sqldf(
      "select ds.*,
                ds2.residual,
                ds3.residualDoublet,
                ds4.residual2Pop,
                ds5.residual3Pop
                from finalData3 ds
                left join initialRSE ds2
                on ds.data = ds2.data
                left join doubletRSE ds3
                on ds.data = ds3.data
                left join twoPopRSE ds4
                on ds.data = ds4.data
                left join threePopRSE ds5
                on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
      residualMultiple=ifelse(
        is.na(residual3Pop),
        residual2Pop,
        residual3Pop
      )
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
      finalRSE=do.call(
        pmin, 
        c(
          subset(
            ., 
            select=c(residual, residualDoublet, residualMultiple)
          ), na.rm=TRUE
        )
      )
    ) %>% dplyr::rename(
      singleRSE=residual,
      doubletRSE=residualDoublet,
      multipleRSE=residualMultiple
    ) %>% dplyr::select(
      -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
      if(
        !is.na(finalData5[i,]$multipleRSE) &
        finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
      ){
        finalData5[i,]$investigate = 1
      }
      
    }
    
    ##Flagging to investigate if the RSE is considered an outlier       
    outlierRSE = as.numeric(
      quantile(finalData5$finalRSE)[4]+
        1.5*stats::IQR(finalData5$finalRSE)
    )
    
    for(i in seq_len(nrow(finalData5))){
      if( finalData5$finalRSE[i] >= outlierRSE ){
        finalData5[i,]$investigate = 1
      }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
      "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
      finalData6=finalData6 %>%
        dplyr::select(
          -c(
            `doublet G1+G2`,
            `doublet G1+G2 count`,
            `doublet G2+G2`,
            `doublet G2+G2 count`
          )
        )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)
    experimentName <- basename(dirname(getwd()))
    experimentName <- sub(" ", "_", experimentName)   # replaces spaces with an underscore
    write.csv(
      finalData6,
      paste0(file.path(dirname(getwd()), subDir), "/", experimentName, "_ploidyPeaksOutput.csv"
      ),
      row.names = FALSE
    )
  }else{
    
    ##Formatting the diploid data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
      dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    ##Formatting the data from peak algorithm 2-4
    finalPart2=finishedDs %>% data.frame() %>%
      dplyr::select(
        -c(
          "g3LL",
          "g3UL",
          "g4LL",
          "g4UL",
          "propCellsUsed",
          "cluster",
          "distToNext",
          "LL",
          "UL"
        )
      ) %>%
      dplyr::mutate(
        investigate=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    ##Formatting the data from the last peak algorthm
    finalPart3=investigateDs %>% data.frame() %>%
      dplyr::select(
        -c(
          "g3LL",
          "g3UL",
          "g4LL",
          "g4UL",
          "cluster",
          "distToNext",
          "LL",
          "UL"
        )
      ) %>%
      dplyr::mutate(
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )
    
    ##Merging the three datasets
    finalData=rbind(
      finalPart1,
      finalPart2,
      finalPart3
    ) %>% dplyr::select(
      c(
        "data",
        "G1",
        "G2",
        "g1G2Doublet",
        "g2G2Doublet",
        "G1Count",
        "G2Count",
        "g1G2DoubletCount",
        "g2G2DoubletCount",
        "investigate"
      )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    if(!purrr::is_empty(peak1Check)){
      finalData$investigate[finalData$data %in% peak1Check$data] = 1
    }
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
      dplyr::select(data, G1, G2, G1Count, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
      dplyr::select(-c(G1, G2,G1Count, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
      dplyr::mutate(id=rowid(data)) %>%
      tidyr::pivot_wider(
        names_from=id,
        values_from=c(G1, G2,G1Count, G2Count)
      )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    
    finalData3=finalData3 %>%
      dplyr::rename(
        `doublet G1+G2`=g1G2Doublet,
        `doublet G1+G2 count`=g1G2DoubletCount,
        `doublet G2+G2`=g2G2Doublet,
        `doublet G2+G2 count`=g2G2DoubletCount
      )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
      doubletRSE<-.popConfidenceDoublet(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      doubletRSE<-finalData3 %>% dplyr::mutate(
        residualDoublet=NA
      )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
      if(TRUE %in% grepl("_3", names(finalData3))){
        ds1<-finalData3 %>% dplyr::filter(is.na(G1_3))    
      }else{
        ds1<-finalData3
      }
      
      if(TRUE %in% !is.na(ds1$G1_2)){
        twoPopRSE<-.popConfidence2Pop(
          flowDir, ds=ds1, xVariable, saveGraph
        )
      }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
          residual2Pop=NA
        )
      }
    }else{
      twoPopRSE<-finalData3 %>% dplyr::mutate(
        residual2Pop=NA
      )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
      threePopRSE<-.popConfidence3Pop(
        flowDir, ds=finalData3, xVariable, saveGraph
      )
    }else{
      threePopRSE<-finalData3 %>% dplyr::mutate(
        residual3Pop=NA
      )
    }
    
    finalData4<-sqldf::sqldf(
      "select ds.*,
                ds2.residual,
                ds3.residualDoublet,
                ds4.residual2Pop,
                ds5.residual3Pop
                from finalData3 ds
                left join initialRSE ds2
                on ds.data = ds2.data
                left join doubletRSE ds3
                on ds.data = ds3.data
                left join twoPopRSE ds4
                on ds.data = ds4.data
                left join threePopRSE ds5
                on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
      residualMultiple=ifelse(
        is.na(residual3Pop),
        residual2Pop,
        residual3Pop
      )
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
      finalRSE=do.call(
        pmin, 
        c(
          subset(
            .,
            select=c(residual, residualDoublet, residualMultiple)
          ), na.rm=TRUE
        )
      )
    ) %>% dplyr::rename(
      singleRSE=residual,
      doubletRSE=residualDoublet,
      multipleRSE=residualMultiple
    ) %>% dplyr::select(
      -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
      if(
        !is.na(finalData5[i,]$multipleRSE) &
        finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
      ){
        finalData5[i,]$investigate = 1
      }
      
    }
    
    ##Flagging to investigate if the RSE is considered an outlier       
    outlierRSE = as.numeric(
      quantile(finalData5$finalRSE)[4]+
        1.5*stats::IQR(finalData5$finalRSE)
    )
    
    for(i in seq_len(nrow(finalData5))){
      if( finalData5$finalRSE[i] >= outlierRSE ){
        finalData5[i,]$investigate = 1
      }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
      "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
      finalData6=finalData6 %>%
        dplyr::select(
          -c(
            `doublet G1+G2`,
            `doublet G1+G2 count`,
            `doublet G2+G2`,
            `doublet G2+G2 count`
          )
        )
    }
    
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings=FALSE)
    experimentName <- basename(dirname(getwd()))
    experimentName <- sub(" ", "_", experimentName)   # replaces spaces with an underscore
    write.csv(
      finalData6,
      paste0(
        file.path(dirname(getwd()), subDir),
        "/", experimentName, "_ploidyPeaksOutput.csv"
      ),
      row.names = FALSE
    )
    
  }
}


## popConfidenceInitial
.popConfidenceInitial = function(flowDir, ds, xVariable, saveGraph = TRUE){
  
  ##Removing NOTE 'no visible binding for global variable'
  G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
  
  modelData <- ds %>% dplyr::select(
    data,
    G1_1,
    G2_1,
    G1Count_1,
    G2Count_1
  )
  
  modelData02 <- modelData %>% dplyr::mutate(
    sdG1Count=G1Count_1*0.6,
    sdG2Count=G2Count_1*0.6
  )
  
  flowNameDs <- unique(modelData02$data)
  
  residualDs <- data.frame(
    matrix(nrow=0, ncol=14)
  )
  
  colnames(residualDs) <- c(
    "data",
    "G1_1",
    "G2_1",
    "G1Count_1",
    "G2Count_1",
    "sdG1Count",
    "sdG2Count",
    "g1Mean",
    "g2Mean",
    "g1SD",
    "g2SD",
    "numG1",
    "numG2",
    "residual"
  )
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    flowData <- .smoothData( flowName, xVariable, 5)
    flowDataMeans <- modelData02%>% dplyr::filter(
      data == flowNameDs[k]
    )
    
    if(0 %in% flowData$x){
      flowData <- flowData %>% dplyr::filter(x != 0)
    }
    
    if(flowDataMeans$G2_1 != 0){
      midPoint <- flowData[
        which(
          abs(
            flowData$x-mean(
              c(flowDataMeans$G1_1, flowDataMeans$G2_1))
          ) == min(
            abs(flowData$x-mean(
              c(flowDataMeans$G1_1, flowDataMeans$G2_1)))
          )
        ),
      ]
      
      if(nrow(midPoint)>1){
        midPoint <- midPoint[nrow(midPoint), ]
      }
      
      #G1 standard deviation 
      g1LeftFlowData <- flowData[
        seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
      g1RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G1_1)+1):
          which(flowData$x == midPoint$x)-1,
      ]
      
      leftPeak1 <- g1LeftFlowData[
        which(
          abs(
            g1LeftFlowData$y-flowDataMeans$sdG1Count
          ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
        ), 
      ]
      
      if(nrow(leftPeak1)>1){
        leftPeak1 <- leftPeak1[1, ]
      }
      
      rightPeak1 <- g1RightFlowData[
        which(
          abs(
            g1RightFlowData$y-flowDataMeans$sdG1Count
          ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      if(nrow(rightPeak1)>1 ){
        rightPeak1 <- rightPeak1[1, ]
      }
      
      #G2 standard deviation
      g2LeftFlowData <- flowData[
        which(
          flowData$x == midPoint$x
        ):which(flowData$x == flowDataMeans$G2_1)-1,
      ]
      g2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G2_1)+1):nrow(flowData),
      ]
      
      leftPeak2 <- g2LeftFlowData[
        which(
          abs(
            g2LeftFlowData$y-flowDataMeans$sdG2Count
          ) == min(abs(g2LeftFlowData$y-flowDataMeans$sdG2Count)
          )
        ),
      ]
      
      if(nrow(leftPeak2)>1 ){
        leftPeak2 <- leftPeak2[1, ]
      }
      
      rightPeak2 <- g2RightFlowData[
        which(
          abs(
            g2RightFlowData$y-flowDataMeans$sdG2Count
          ) == min(abs(g2RightFlowData$y-flowDataMeans$sdG2Count)
          )
        ),
      ]
      
      if(nrow(rightPeak2)>1 ){
        rightPeak2 <- rightPeak2[1, ]
      }
      
      flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
        g1Mean=G1_1,
        g2Mean=G2_1,
        g1SD=mean(
          c(
            rightPeak1$x-flowDataMeans$G1_1,
            flowDataMeans$G1_1-leftPeak1$x
          )
        ),
        g2SD=mean(
          c(
            rightPeak2$x-flowDataMeans$G2_1,
            flowDataMeans$G2_1-leftPeak2$x
          )
        ),
        numG1=sum(
          flowData$y[
            c(
              which(flowData$x == leftPeak1$x):
                which(flowData$x == rightPeak1$x)
            )
          ]
        ),
        numG2=sum(
          flowData$y[
            c(
              which(flowData$x == leftPeak2$x):
                which(flowData$x == rightPeak2$x)
            )
          ]
        )
      )
      g1Mean <- flowDataMeans01$g1Mean
      g2Mean <- flowDataMeans01$g2Mean
      g1SD <- flowDataMeans01$g1SD
      g2SD <- flowDataMeans01$g2SD
      
      xpectr::suppress_mw(
        singlePopNLS <- nls(
          formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                           exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
            (N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
            (A + B*x + C*(x^2))*
            (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                  (2*(g1SD*(x/g1Mean))^2))),
          data=flowData,
          start=c(
            N1=flowDataMeans01$numG1,
            N2=flowDataMeans01$numG2,
            A=0, B=0, C=0
          ),
          nls.control(warnOnly=TRUE)
        )
      )
    }else{
      #G1 standard deviation 
      g1LeftFlowData <- flowData[
        seq_len(which(flowData$x == flowDataMeans$G1_1)-1),
      ]
      g1RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G1_1)+1):nrow(flowData),
      ]
      
      leftPeak1 <- g1LeftFlowData[
        which(
          abs(
            g1LeftFlowData$y-flowDataMeans$sdG1Count
          ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      rightPeak1 <- g1RightFlowData[
        which(
          abs(
            g1RightFlowData$y-flowDataMeans$sdG1Count
          ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
        g1Mean=G1_1,
        g2Mean=G2_1,
        g1SD=mean(
          c(
            rightPeak1$x-flowDataMeans$G1_1,
            flowDataMeans$G1_1-leftPeak1$x
          )
        ),
        g2SD=NA,
        numG1=sum(
          flowData$y[
            c(
              which(flowData$x == leftPeak1$x):
                which(flowData$x == rightPeak1$x)
            )
          ]
        ),
        numG2=0
      )
      
      g1Mean <- flowDataMeans01$g1Mean
      g1SD <- flowDataMeans01$g1SD
      
      xpectr::suppress_mw(
        singlePopNLS <- nls(
          formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                           exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
            (A + B*x + C*(x^2))*
            (1/(sqrt(2 * pi)*g1SD*(x/g1Mean))*
               exp(-((x-g1Mean)^2)/(2*(g1SD*(x/g1Mean))^2))),
          data=flowData,
          start=c(
            N1=flowDataMeans01$numG1,
            A=0, B=0, C=0
          ),
          nls.control(warnOnly=TRUE)
        )
      )
    }
    
    if(saveGraph == TRUE){
      nlsGraph <- ggplot()+
        geom_line(data=flowData, aes(x=x, y=y, col ='Raw Data'), linewidth=1)+
        geom_line(aes(flowData$x,predict(singlePopNLS), col='NLS'))+
        labs(y="Counts", x=xVariable)+
        theme(legend.title=element_blank())
      
      plotDir <- "nlsGraphs"
      dir.create(file.path(dirname(flowDir), plotDir), showWarnings=FALSE)
      plotInitDir <- "nlsSingle_graphs"
      dir.create(
        file.path(dirname(flowDir) ,plotDir, plotInitDir),
        showWarnings=FALSE
      )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
      
      png(
        paste0(plotOutFile, "/", flowNameDs[k], '.png'),
        width=600, height=400
      )
      print(nlsGraph)
      dev.off()
    }
    
    flowDataMeans01$residual <- summary(singlePopNLS)[["sigma"]]
    
    residualDs <- rbind(
      residualDs,
      flowDataMeans01
    )
    
  }
  
  return(residualDs)
}

## popConfidenceDoublet
.popConfidenceDoublet = function(flowDir, ds, xVariable, saveGraph = TRUE){
  ##Removing NOTE 'no visible binding for global variable'
  doublet<-G1_1<-G2_1<-`doublet G1+G2`<-`doublet G2+G2`<-G1Count_1<-NULL
  G2Count_1<-`doublet G1+G2 count`<-`doublet G2+G2 count`<-x<-y<-NULL
  
  ds2 <- ds %>% dplyr::filter(doublet == 1)
  flowNameDs <- unique(ds2$data)
  
  modelData01 <- ds2 %>% dplyr::select(
    data,
    G1_1,
    G2_1,
    `doublet G1+G2`,
    `doublet G2+G2`,
    G1Count_1,
    G2Count_1,
    `doublet G1+G2 count`,
    `doublet G2+G2 count`,
    doublet
  )
  
  modelData02 <- modelData01 %>% dplyr::mutate(
    sdG1Count=G1Count_1*0.6,
    sdG2Count=G2Count_1*0.6,
    sdG1G2Count=`doublet G1+G2 count`*0.6,
    sdG2G2Count=`doublet G2+G2 count`*0.6
  )
  
  residualDoubletDs <- data.frame(
    matrix(nrow=0, ncol=26)
  )
  colnames(residualDoubletDs) <- c(
    "data",
    "G1_1",
    "G2_1",
    "doublet G1+G2",
    "doublet G2+G2",
    "G1Count_1",
    "G2Count_1",
    "doublet G1+G2 count",
    "doublet G2+G2 count",
    "sdG1Count",
    "sdG2Count",
    "sdG1G2Count",
    "sdG2G2Count",
    "g1Mean",
    "g2Mean",
    "G1G2Mean",
    "G2G2Mean",
    "g1SD",
    "g2SD",
    "G1g2SD",
    "G2g2SD",
    "numG1",
    "numG2",
    "numDoublet1",
    "numDoublet2",
    "residualDoublet"
  )
  
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    flowData <- .smoothData( flowName, xVariable, 5)
    flowDataMeans <- modelData02%>% dplyr::filter(
      data == flowNameDs[k]
    )
    
    if(0 %in% flowData$x){
      flowData <- flowData %>% dplyr::filter(x != 0)
    }
    
    if(flowDataMeans$G2_1 != 0){
      midPoint <- flowData[
        which(
          abs(
            flowData$x-mean(c(flowDataMeans$G1_1,
                              flowDataMeans$G2_1))
          ) == min(
            abs(
              flowData$x-mean(c(flowDataMeans$G1_1,
                                flowDataMeans$G2_1))
            )
          )
        ),
      ]
      
      if(nrow(midPoint)>1){
        midPoint <- midPoint[nrow(midPoint), ]
      }
      
      #G1 standard deviation 
      G1LeftFlowData <- flowData[
        seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
      G1RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G1_1)+1):
          which(flowData$x == midPoint$x),
      ]
      
      leftPeak1 <- G1LeftFlowData[
        which(
          abs(G1LeftFlowData$y-flowDataMeans$sdG1Count) == 
            min(abs(G1LeftFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      if(nrow(leftPeak1)>1 ){
        leftPeak1 <- leftPeak1[1, ]
      }
      
      rightPeak1 <- G1RightFlowData[
        which(
          abs(G1RightFlowData$y-flowDataMeans$sdG1Count) == 
            min(abs(G1RightFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      if(nrow(rightPeak1)>1 ){
        rightPeak1 <- rightPeak1[1, ]
      }
      
      #G2 standard deviation
      G2LeftFlowData <- flowData[
        which(flowData$x == midPoint$x):
          (which(flowData$x == flowDataMeans$G2_1)-1),
      ]
      G2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G2_1)+1):nrow(flowData),
      ]
      
      leftPeak2 <- G2LeftFlowData[
        which(
          abs(G2LeftFlowData$y-flowDataMeans$sdG2Count) == 
            min(abs(G2LeftFlowData$y-flowDataMeans$sdG2Count))
        ),
      ]
      
      if(nrow(leftPeak2)>1 ){
        leftPeak2 <- leftPeak2[1, ]
      }
      
      rightPeak2 <- G2RightFlowData[
        which(
          abs(G2RightFlowData$y-flowDataMeans$sdG2Count) == 
            min(abs(G2RightFlowData$y-flowDataMeans$sdG2Count))
        ),
      ]
      
      if(nrow(rightPeak2)>1 ){
        rightPeak2 <- rightPeak2[1, ]
      }
      
      #G1+G2 doublet standard deviation
      midPointDoublet <- flowData[
        which(
          abs(
            flowData$x-mean(c(flowDataMeans$`doublet G1+G2`,
                              flowDataMeans$G2_1))
          ) == min(
            abs(
              flowData$x-mean(c(flowDataMeans$`doublet G1+G2`,
                                flowDataMeans$G2_1))
            )
          )
        ),
      ]
      
      if(nrow(midPointDoublet)>1){
        midPointDoublet <- midPointDoublet[nrow(midPointDoublet), ]
      }
      
      doubletG1G2LeftFlowData <- flowData[
        which(flowData$x == midPointDoublet$x):
          (which(flowData$x == flowDataMeans$`doublet G1+G2`)-1),
      ]
      doubletG1G2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$`doublet G1+G2`)+1):
          nrow(flowData),
      ]
      
      doubletG1G2Left <- doubletG1G2LeftFlowData[
        which(
          abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count) == 
            min(
              abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count)
            )
        ),
      ]
      
      if(nrow(doubletG1G2Left)>1){
        doubletG1G2Left <- doubletG1G2Left[1, ]
      }
      
      doubletG1G2Right <- doubletG1G2RightFlowData[
        which(
          abs(
            doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count
          ) == min(
            abs(
              doubletG1G2RightFlowData$y-
                flowDataMeans$sdG1G2Count
            )
          )
        ),
      ]
      
      if(nrow(doubletG1G2Right)>1){
        doubletG1G2Right <- doubletG1G2Right[1, ]
      }
      
      #G2+G2 doublet standard deviation
      if(!is.na(flowDataMeans$`doublet G2+G2`)){
        midPointDoublet2 <- flowData[
          which(
            abs(
              flowData$x-mean(
                c(
                  flowDataMeans$`doublet G1+G2`,
                  flowDataMeans$`doublet G2+G2`
                )
              )
            ) == min(
              abs(
                flowData$x-mean(
                  c(
                    flowDataMeans$`doublet G1+G2`,
                    flowDataMeans$`doublet G2+G2`
                  )
                )
              )
            )
          ),
        ]
        
        if(nrow(midPointDoublet2)>1){
          midPointDoublet2 <- midPointDoublet2[
            nrow(midPointDoublet2), ]
        }
        
        doubletG1G2LeftFlowData <- flowData[
          which(flowData$x == midPointDoublet2$x):
            (which(flowData$x == flowDataMeans$`doublet G2+G2`)-1),
        ]
        
        doubletG2G2RightFlowData <- flowData[
          (which(flowData$x == flowDataMeans$`doublet G2+G2`)+1):
            nrow(flowData),
        ]
        
        doubletG2G2Left <- doubletG1G2LeftFlowData[
          which(
            abs(
              doubletG1G2LeftFlowData$y-
                flowDataMeans$sdG2G2Count
            ) == min(
              abs(
                doubletG1G2LeftFlowData$y-
                  flowDataMeans$sdG2G2Count
              )
            )
          ),
        ]
        
        if(nrow(doubletG2G2Left)>1){
          doubletG2G2Left <- doubletG2G2Left[1, ]
        }
        
        doubletG2G2Right <- doubletG2G2RightFlowData[
          which(
            abs(
              doubletG2G2RightFlowData$y-flowDataMeans$sdG2G2Count
            ) == min(
              abs(
                doubletG2G2RightFlowData$y-
                  flowDataMeans$sdG2G2Count
              )
            )
          ),
        ]
        
        if(nrow(doubletG2G2Right)>1){
          doubletG2G2Right <- doubletG2G2Right[1, ]
        }
        flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
          g1Mean=G1_1,
          g2Mean=G2_1,
          G1G2Mean=`doublet G1+G2`,
          G2G2Mean=`doublet G2+G2`,
          g1SD=mean(
            c(
              rightPeak1$x-flowDataMeans$G1_1,
              flowDataMeans$G1_1-leftPeak1$x
            )
          ),
          g2SD=mean(
            c(
              rightPeak2$x-flowDataMeans$G2_1,
              flowDataMeans$G2_1-leftPeak2$x
            )
          ),
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2` - doubletG1G2Left$x
            )
          ),
          G2g2SD=mean(
            c(
              doubletG2G2Right$x-flowDataMeans$`doublet G2+G2`,
              flowDataMeans$`doublet G2+G2` - doubletG2G2Left$x
            )
          ),
          numG1=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak1$x):
                  which(flowData$x == rightPeak1$x)
              )
            ]
          ),
          numG2=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak2$x):
                  which(flowData$x == rightPeak2$x)
              )
            ]
          ),
          numDoublet1=sum(
            flowData$y[
              c(
                which(flowData$x == doubletG1G2Left$x):
                  which(flowData$x == doubletG1G2Right$x)
              )
            ]
          ),
          numDoublet2=sum(
            flowData$y[
              c(
                which(flowData$x == doubletG2G2Left$x):
                  which(flowData$x == doubletG2G2Right$x)
              )
            ]
          )
        )
        
        g1Mean <- flowDataMeans01$g1Mean
        g2Mean <- flowDataMeans01$g2Mean
        G1G2Mean <- flowDataMeans01$G1G2Mean
        G2G2Mean <- flowDataMeans01$G2G2Mean
        
        g1SD <- flowDataMeans01$g1SD
        g2SD <- flowDataMeans01$g2SD
        G1g2SD <- flowDataMeans01$G1g2SD
        G2g2SD <- flowDataMeans01$G2g2SD
        
        xpectr::suppress_mw(
          singlePopNLS <- nls(
            formula = y ~ (N1/(sqrt(2*pi) * g1SD) * 
                             exp(-((x-g1Mean)^2)/(2 *g1SD^2))) +
              (A1 + B1*x + C1*(x^2))*
              (1/(sqrt(2 * pi) * g1SD * (x/g1Mean)) * 
                 exp(-((x - g1Mean)^2)/(2 *(g1SD* (x/g1Mean))^2)))+
              (N2/(sqrt(2 * pi) * g2SD) * exp(-((x - g2Mean)^2)/
                                                (2 *g2SD^2)))+
              (A2 + B2*x + C2*(x^2))*
              (1/(sqrt(2 * pi) * g2SD * (x/g2Mean)) * 
                 exp(-((x - g2Mean)^2)/(2 *(g2SD* (x/g2Mean))^2)))+
              (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
                                                     (2*G1g2SD^2)))+
              (A3 + B3*x + C3*(x^2))*
              (1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*
                 exp(-((x-G1G2Mean)^2)/(2*(G1g2SD*(x/G1G2Mean))^2)))+
              (numDoublet2/(sqrt(2*pi)*G2g2SD)*exp(-((x-G2G2Mean)^2)/
                                                     (2 *G2g2SD^2))),
            data=flowData,
            start=c(
              N1=flowDataMeans01$numG1,
              N2=flowDataMeans01$numG2,
              numDoublet1=flowDataMeans01$numDoublet1,
              numDoublet2=flowDataMeans01$numDoublet2,
              A1=0, B1=0, C1=0,
              A2=0, B2=0, C2=0,
              A3=0, B3=0, C3=0
            ),
            nls.control(warnOnly=TRUE)
          )
        )
      }else{
        
        flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
          g1Mean=G1_1,
          g2Mean=G2_1,
          G1G2Mean=`doublet G1+G2`,
          G2G2Mean=`doublet G2+G2`,
          g1SD=mean(
            c(
              rightPeak1$x-flowDataMeans$G1_1,
              flowDataMeans$G1_1-leftPeak1$x
            )
          ),
          g2SD=mean(
            c(
              rightPeak2$x-flowDataMeans$G2_1,
              flowDataMeans$G2_1-leftPeak2$x
            )
          ),
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x
            )
          ),
          G2g2SD=NA,
          numG1=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak1$x):
                  which(flowData$x == rightPeak1$x)
              )
            ]
          ),
          numG2=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak2$x):
                  which(flowData$x == rightPeak2$x)
              )
            ]
          ),
          numDoublet1=sum(
            flowData$y[
              c(
                which(flowData$x == doubletG1G2Left$x):
                  which(flowData$x == doubletG1G2Right$x)
              )
            ]
          ),
          numDoublet2=0
        )
        
        g1Mean <- flowDataMeans01$g1Mean
        g2Mean <- flowDataMeans01$g2Mean
        G1G2Mean <- flowDataMeans01$G1G2Mean
        g1SD <- flowDataMeans01$g1SD
        g2SD <- flowDataMeans01$g2SD
        G1g2SD <- flowDataMeans01$G1g2SD
        
        xpectr::suppress_mw(
          singlePopNLS <- nls(
            formula=y ~ (N1/(sqrt(2*pi)*g1SD)*
                           exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
              (A1 + B1*x + C1*(x^2))*
              (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                    (2*(g1SD*(x/g1Mean))^2)))+
              (N2/(sqrt(2*pi)*g2SD)*exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
              (A2 + B2*x + C2*(x^2))*
              (1/(sqrt(2*pi)*g2SD*(x/g2Mean))*exp(-((x-g2Mean)^2)/
                                                    (2*(g2SD*(x/g2Mean))^2)))+
              (numDoublet1/(sqrt(2*pi)*G1g2SD)*
                 exp(-((x-G1G2Mean)^2)/(2*G1g2SD^2))),
            data=flowData,
            start=c(
              N1=flowDataMeans01$numG1,
              N2=flowDataMeans01$numG2,
              numDoublet1=flowDataMeans01$numDoublet1,
              A1=0, B1=0, C1=0,
              A2=0, B2=0, C2=0
            ),
            nls.control(warnOnly=TRUE)
          )
        )
        
      }
      
    }else{
      
      #G1 standard deviation 
      G1LeftFlowData <- flowData[
        seq_len(which(flowData$x == flowDataMeans$G1_1)), ]
      G1RightFlowData <- flowData[
        which(flowData$x == flowDataMeans$G1_1):nrow(flowData),
      ]
      
      leftPeak1 <- G1LeftFlowData[
        which(
          abs(
            G1LeftFlowData$y-flowDataMeans$sdG1Count
          ) == min(abs(G1LeftFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      if(nrow(leftPeak1)>1 ){
        leftPeak1 <- leftPeak1[1, ]
      }
      
      rightPeak1 <- G1RightFlowData[
        which(
          abs(
            G1RightFlowData$y-flowDataMeans$sdG1Count
          ) == min(abs(G1RightFlowData$y-flowDataMeans$sdG1Count))
        ),
      ]
      
      if(nrow(rightPeak1)>1 ){
        rightPeak1 <- rightPeak1[1, ]
      }
      
      
      #G1+G2 doublet standard deviation
      midPointDoublet <- flowData[
        which(
          abs(
            flowData$x-mean(
              c(flowDataMeans$`doublet G1+G2`, flowDataMeans$G1_1)
            )
          ) == min(
            abs(
              flowData$x-mean(
                c(
                  flowDataMeans$`doublet G1+G2`,
                  flowDataMeans$G1_1)
              )
            )
          )
        ),
      ]
      
      if(nrow(midPointDoublet)>1){
        midPointDoublet <- midPointDoublet[nrow(midPointDoublet), ]
      }
      
      doubletG1G2LeftFlowData <- flowData[
        which(flowData$x == midPointDoublet$x):
          which(flowData$x == flowDataMeans$`doublet G1+G2`),
      ]
      
      doubletG1G2RightFlowData <- flowData[
        which(flowData$x == flowDataMeans$`doublet G1+G2`):
          nrow(flowData),
      ]
      
      doubletG1G2Left <- doubletG1G2LeftFlowData[
        which(
          abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count) == 
            min(
              abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count)
            )
        ),
      ]
      
      if(nrow(doubletG1G2Left)>1){
        doubletG1G2Left <- doubletG1G2Left[1, ]
      }
      
      doubletG1G2Right <- doubletG1G2RightFlowData[
        which(
          abs(
            doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count
          ) == min(
            abs(
              doubletG1G2RightFlowData$y-
                flowDataMeans$sdG1G2Count
            )
          )
        ),
      ]
      
      if(nrow(doubletG1G2Right)>1){
        doubletG1G2Right <- doubletG1G2Right[1, ]
      }
      
      if(!is.na(flowDataMeans$`doublet G2+G2`)){
        #G2+G2 doublet standard deviation
        midPointDoublet2 <- flowData[
          which(
            abs(
              flowData$x-mean(
                c(
                  flowDataMeans$`doublet G1+G2`,
                  flowDataMeans$`doublet G2+G2`)
              )
            ) == 
              min(
                abs(
                  flowData$x-mean(
                    c(
                      flowDataMeans$`doublet G1+G2`,
                      flowDataMeans$`doublet G2+G2`
                    )
                  )
                )
              )
          ),
        ]
        
        if(nrow(midPointDoublet2)>1){
          midPointDoublet2 <- midPointDoublet2[
            nrow(midPointDoublet2), ]
        }
        
        doubletG1G2LeftFlowData <- flowData[
          which(flowData$x == midPointDoublet2$x):
            which(flowData$x == flowDataMeans$`doublet G2+G2`), 
        ]
        doubletG2G2RightFlowData <- flowData[
          which(flowData$x == flowDataMeans$`doublet G2+G2`):
            nrow(flowData),
        ]
        
        doubletG2G2Left <- doubletG1G2LeftFlowData[
          which(
            abs(
              doubletG1G2LeftFlowData$y-
                flowDataMeans$sdG2G2Count
            ) == min(
              abs(
                doubletG1G2LeftFlowData$y-
                  flowDataMeans$sdG2G2Count
              )
            )
          ),
        ]
        
        if(nrow(doubletG2G2Left)>1){
          doubletG2G2Left <- doubletG2G2Left[1, ]
        }
        
        doubletG2G2Right <- doubletG2G2RightFlowData[
          which(
            abs(
              doubletG2G2RightFlowData$y-
                flowDataMeans$sdG2G2Count
            ) == min(
              abs(
                doubletG2G2RightFlowData$y-
                  flowDataMeans$sdG2G2Count
              )
            )
          ),
        ]
        
        if(nrow(doubletG2G2Right)>1){
          doubletG2G2Right <- doubletG2G2Right[1, ]
        }
        flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
          g1Mean=G1_1,
          g2Mean=G2_1,
          G1G2Mean=`doublet G1+G2`,
          G2G2Mean=`doublet G2+G2`,
          g1SD=mean(
            c(
              rightPeak1$x-flowDataMeans$G1_1,
              flowDataMeans$G1_1-leftPeak1$x
            )
          ),
          g2SD=NA,
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x
            )
          ),
          G2g2SD=mean(
            c(
              doubletG2G2Right$x-flowDataMeans$`doublet G2+G2`,
              flowDataMeans$`doublet G2+G2`-doubletG2G2Left$x
            )
          ),
          numG1=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak1$x):
                  which(flowData$x == rightPeak1$x)
              )
            ]
          ),
          numG2=0,
          numDoublet1=sum(
            flowData$y[
              c(
                which(flowData$x == doubletG1G2Left$x):
                  which(flowData$x == doubletG1G2Right$x)
              )
            ]
          ),
          numDoublet2=sum(
            flowData$y[
              c(
                which(flowData$x == doubletG2G2Left$x):
                  which(flowData$x == doubletG2G2Right$x)
              )
            ]
          )
        )
        
        g1Mean <- flowDataMeans01$g1Mean
        G1G2Mean <- flowDataMeans01$G1G2Mean
        G2G2Mean <- flowDataMeans01$G2G2Mean
        g1SD <- flowDataMeans01$g1SD
        G1g2SD <- flowDataMeans01$G1g2SD
        G2g2SD <- flowDataMeans01$G2g2SD
        
        xpectr::suppress_mw(
          singlePopNLS <- nls(
            formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                             exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
              (A1 + B1*x + C1*(x^2))*
              (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x - g1Mean)^2)/
                                                    (2*(g1SD*(x/g1Mean))^2)))+
              (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
                                                     (2*G1g2SD^2)))+
              (A2 + B2*x + C2*(x^2))*
              (1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*
                 exp(-((x-G1G2Mean)^2)/(2*(G1g2SD*(x/G1G2Mean))^2)))+
              (numDoublet2/(sqrt(2*pi)*G2g2SD)*exp(-((x-G2G2Mean)^2)/
                                                     (2*G2g2SD^2))),
            data=flowData,
            start=c(
              N1=flowDataMeans01$numG1,
              numDoublet1=flowDataMeans01$numDoublet1,
              numDoublet2=flowDataMeans01$numDoublet2,
              A1=0, B1=0, C1=0,
              A2=0, B2=0, C2=0
            ),
            nls.control(warnOnly=TRUE)
          )
        )
      }else{
        
        flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
          g1Mean=G1_1,
          g2Mean=G2_1,
          G1G2Mean=`doublet G1+G2`,
          G2G2Mean=`doublet G2+G2`,
          g1SD=mean(
            c(
              rightPeak1$x-flowDataMeans$G1_1,
              flowDataMeans$G1_1-leftPeak1$x
            )
          ),
          g2SD=NA,
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x
            )
          ),
          G2g2SD=NA,
          numG1=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak1$x):
                  which(flowData$x == rightPeak1$x)
              )
            ]
          ),
          numG2=0,
          numDoublet1=sum(
            flowData$y[
              c(
                which(flowData$x == doubletG1G2Left$x):
                  which(flowData$x == doubletG1G2Right$x)
              )
            ]
          ),
          numDoublet2=0
        )
        
        g1Mean <- flowDataMeans01$g1Mean
        G1G2Mean <- flowDataMeans01$G1G2Mean
        g1SD <- flowDataMeans01$g1SD
        G1g2SD <- flowDataMeans01$G1g2SD
        
        xpectr::suppress_mw(
          singlePopNLS <- nls(
            formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                             exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
              (A1 + B1*x + C1*(x^2))*
              (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                    (2*(g1SD*(x/g1Mean))^2)))+
              (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
                                                     (2*G1g2SD^2))),
            data=flowData,
            start=c(
              N1=flowDataMeans01$numG1,
              numDoublet1=flowDataMeans01$numDoublet1,
              A1=0, B1=0, C1=0
            ),
            nls.control(warnOnly=TRUE)
          )
        )
        
      }
      
    }
    
    if(saveGraph == TRUE){
      nlsGraph <- ggplot()+
        geom_line(data=flowData, aes(x=x, y=y, col = 'Raw Data'), linewidth=1)+
        geom_line(aes(flowData$x, predict(singlePopNLS), col ='NLS'))+
        labs(y ="Counts", x =xVariable)+
        theme(legend.title = element_blank())
      
      plotDir <- "nlsGraphs"
      plotInitDir <- "nlsDoublet_graphs"
      dir.create(
        file.path(dirname(flowDir), plotDir, plotInitDir),
        showWarnings=FALSE
      )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
      
      png(
        paste0(plotOutFile, "/", flowNameDs[k], '.png'),
        width=600, height=400
      )
      print(nlsGraph)
      dev.off()
    }
    
    flowDataMeans01$residualDoublet <- summary(singlePopNLS)[["sigma"]]
    
    residualDoubletDs <- rbind(
      residualDoubletDs,
      flowDataMeans01
    )
    
  }
  return(residualDoubletDs)
  
}

## popConfidence2Pop
.popConfidence2Pop = function(flowDir, ds, xVariable, saveGraph = TRUE){
  
  ##Removing NOTE 'no visible binding for global variable'
  G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
  G1_2<-G2_2<-G1_3<-G1Count_2<-G2Count_2<-NULL
  
  modelData <- ds %>% dplyr::select(
    data,
    G1_1,
    G2_1,
    G1_2,
    G2_2,
    G1Count_1,
    G2Count_1,
    G1Count_2,
    G2Count_2
  )
  
  modelData01<-modelData %>% dplyr::filter(
    !is.na(modelData$G1_2)
  )
  modelData02 <- modelData01 %>% dplyr::mutate(
    sdG1Count=G1Count_1*0.6,
    sdG2Count=G2Count_1*0.6,
    sdG1Count2=G1Count_2*0.6,
    sdG2Count2=G2Count_2*0.6
  )
  
  flowNameDs <- unique(modelData02$data)
  
  residualMultiDs <- data.frame(
    matrix(nrow=0, ncol=26)
  )
  
  colnames(residualMultiDs) <- c(
    "data",
    "G1_1",
    "G2_1",
    "G1_2",
    "G2_2", 
    "G1Count_1",
    "G2Count_1",
    "G1Count_2",
    "G2Count_2", 
    "sdG1Count",
    "sdG2Count",
    "sdG1Count2",
    "sdG2Count2",
    "g1Mean",
    "g2Mean",
    "g1SD",
    "g2SD",
    "numG1",
    "numG2",
    "g1Mean2",
    "g2Mean2",
    "g1SD2",
    "g2SD2",
    "pop2NumG1",
    "pop2NumG2",
    "residual"
  )
  
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    flowData <- .smoothData( flowName, xVariable, 5)
    flowDataMeans <- modelData02%>% dplyr::filter(
      data == flowNameDs[k]
    )
    
    if(0 %in% flowData$x){
      flowData <- flowData %>% dplyr::filter(x != 0)
    }
    
    midPoint <- flowData[
      which(
        abs(
          flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2))
        ) == min(
          abs(flowData$x-mean(
            c(flowDataMeans$G1_1,flowDataMeans$G1_2)
          )
          )
        )
      ),
    ]
    
    if(nrow(midPoint)>1){
      midPoint <- midPoint[nrow(midPoint), ]
    }
    
    #G1 standard deviation 
    g1LeftFlowData <- flowData[
      seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
    g1RightFlowData <- flowData[
      (which(flowData$x == flowDataMeans$G1_1)+1):
        which(flowData$x == midPoint$x)-1,
    ]
    
    leftPeak1 <- g1LeftFlowData[
      which(
        abs(
          g1LeftFlowData$y-flowDataMeans$sdG1Count
        ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
      ),
    ]
    
    if(nrow(leftPeak1)>1){
      leftPeak1 <- leftPeak1[1, ]
    }
    
    rightPeak1 <- g1RightFlowData[
      which(
        abs(
          g1RightFlowData$y-flowDataMeans$sdG1Count
        ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
      ),
    ]
    
    if(nrow(rightPeak1)>1 ){
      rightPeak1 <- rightPeak1[1, ]
    }
    midPointLeft <- flowData[
      which(
        abs(
          flowData$x-mean(c(flowDataMeans$G2_1, flowDataMeans$G1_2))
        ) == min(
          abs(
            flowData$x-mean(c(flowDataMeans$G2_1,flowDataMeans$G1_2))
          )
        )
      ),
    ]
    
    if(nrow(midPointLeft)>1){
      midPointLeft <- midPointLeft[nrow(midPointLeft), ]
    }
    
    midPointRight <- flowData[
      which(
        abs(
          flowData$x-mean(c(flowDataMeans$G2_1, flowDataMeans$G2_2))
        ) == min(
          abs(
            flowData$x-
              mean(c(flowDataMeans$G2_1,flowDataMeans$G2_2))
          )
        )
      ),
    ]
    
    if(nrow(midPointRight)>1){
      midPointRight <- midPointRight[nrow(midPointRight), ]
    }
    
    #G2 standard deviation
    if(
      length(
        unique(
          c(
            flowDataMeans$G1_1,
            flowDataMeans$G1_2,
            flowDataMeans$G2_2,
            flowDataMeans$G2_1 
          )
        )
      ) == 3
    ){
      g2LeftFlowData <- flowData[
        which(flowData$x == midPoint$x):
          which(flowData$x == flowDataMeans$G2_1)-1,
      ]
      g2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G2_1)+1):
          nrow(flowData),
      ]
      
      leftPeak2 <- g2LeftFlowData[
        which(
          abs(g2LeftFlowData$y-flowDataMeans$sdG2Count) == 
            min(abs(g2LeftFlowData$y-flowDataMeans$sdG2Count))
        ),
      ]
      
      if(nrow(leftPeak2)>1 ){
        leftPeak2 <- leftPeak2[1, ]
      }
    }else{
      g2LeftFlowData <- flowData[
        which(flowData$x == midPointLeft$x):
          which(flowData$x == flowDataMeans$G2_1)-1,
      ]
      g2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G2_1)+1):
          which(flowData$x == midPointRight$x),
      ]
      
      leftPeak2 <- g2LeftFlowData[
        which(
          abs(g2LeftFlowData$y-flowDataMeans$sdG2Count) == 
            min(abs(g2LeftFlowData$y-flowDataMeans$sdG2Count))
        ),
      ]
      
      if(nrow(leftPeak2)>1 ){
        leftPeak2 <- leftPeak2[1, ]
      }
    }
    
    rightPeak2 <- g2RightFlowData[
      which(
        abs(g2RightFlowData$y-flowDataMeans$sdG2Count) == 
          min(abs(g2RightFlowData$y-flowDataMeans$sdG2Count))
      ),
    ]
    
    if(nrow(rightPeak2)>1 ){
      rightPeak2 <- rightPeak2[1, ]
    }
    
    #pop2
    
    #G1 standard deviation
    g1LeftFlowData2 <- flowData[
      (which(flowData$x == midPoint$x)+1):
        which(flowData$x == flowDataMeans$G1_2)-1,
    ]
    
    g1RightFlowData2 <- flowData[
      (which(flowData$x == flowDataMeans$G1_2)+1):
        which(flowData$x == midPointLeft$x)-1,
    ]
    
    pop2leftPeak1 <- g1LeftFlowData2[
      which(
        abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2) == 
          min(abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2))
      ),
    ]
    
    if(nrow(pop2leftPeak1)>1){
      pop2leftPeak1 <- pop2leftPeak1[1, ]
    }
    
    pop2rightPeak1 <- g1RightFlowData2[
      which(
        abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2) == 
          min(abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2))
      ),
    ]
    
    if(nrow(pop2rightPeak1)>1 ){
      pop2rightPeak1 <- pop2rightPeak1[1, ]
    }
    
    #G2 standard deviation
    g2LeftFlowData2 <- flowData[
      which(flowData$x == midPointRight$x):
        which(flowData$x == flowDataMeans$G2_2)-1,
    ]
    g2RightFlowData2 <- flowData[
      (which(flowData$x == flowDataMeans$G2_2)+1):nrow(flowData),
    ]
    
    pop2leftPeak2 <- g2LeftFlowData2[
      which(
        abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2) == 
          min(abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2))
      ),
    ]
    
    if(nrow(pop2leftPeak2)>1 ){
      pop2leftPeak2 <- pop2leftPeak2[1, ]
    }
    
    pop2rightPeak2 <- g2RightFlowData2[
      which(
        abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2) == 
          min(abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2))
      ),
    ]
    
    if(nrow(pop2rightPeak2)>1 ){
      pop2rightPeak2 <- pop2rightPeak2[1, ]
    }
    
    flowDataMeans01 = flowDataMeans %>% dplyr::mutate(
      g1Mean=G1_1,
      g2Mean=G2_1,
      g1SD=mean(
        c(
          rightPeak1$x-flowDataMeans$G1_1,
          flowDataMeans$G1_1-leftPeak1$x
        )
      ),
      g2SD=mean(
        c(
          rightPeak2$x-flowDataMeans$G2_1,
          flowDataMeans$G2_1-leftPeak2$x
        )
      ),
      pop1NumG1=sum(
        flowData$y[
          c(
            which(flowData$x == leftPeak1$x):
              which(flowData$x == rightPeak1$x)
          )
        ]
      ),
      pop1NumG2=sum(
        flowData$y[
          c(
            which(flowData$x == leftPeak2$x):
              which(flowData$x == rightPeak2$x)
          )
        ]
      ),
      g1Mean2=G1_2,
      g2Mean2=G2_2,
      g1SD2=mean(
        c(
          abs(pop2rightPeak1$x-flowDataMeans$G1_2),
          flowDataMeans$G1_2-pop2leftPeak1$x
        )
      ),
      g2SD2=mean(
        c(
          pop2rightPeak2$x-flowDataMeans$G2_2,
          flowDataMeans$G2_2-pop2leftPeak2$x
        )
      ),
      pop2NumG1=sum(
        flowData$y[
          c(
            which(flowData$x == pop2leftPeak1$x):
              which(flowData$x == pop2rightPeak1$x)
          )
        ]
      ),
      pop2NumG2=sum(
        flowData$y[
          c(
            which(flowData$x == pop2leftPeak2$x):
              which(flowData$x == pop2rightPeak2$x)
          )
        ]
      )
    )
    
    g1Mean <- flowDataMeans01$g1Mean
    g2Mean <- flowDataMeans01$g2Mean
    g1SD <- flowDataMeans01$g1SD
    g2SD <- flowDataMeans01$g2SD
    g1Mean2 <- flowDataMeans01$g1Mean2
    g2Mean2 <- flowDataMeans01$g2Mean2
    g1SD2 <- flowDataMeans01$g1SD2
    g2SD2 <- flowDataMeans01$g2SD2
    
    if(
      length(
        unique(
          c(
            flowDataMeans$G1_1,
            flowDataMeans$G1_2,
            flowDataMeans$G2_2,
            flowDataMeans$G2_1 
          )
        )
      ) == 3
    ){
      xpectr::suppress_mw(
        multiPopNLS <- nls(
          formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*
                           exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
            (A + B*x + C*(x^2))*
            (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                  (2*(g1SD*(x/g1Mean))^2)))+
            (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
            (A2 + B2*x + C2*(x^2))*
            (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                                                    (2*(g2SD2*(x/g2Mean2))^2)))+
            (pop2N2/(sqrt(2*pi)*g2SD2)*
               exp(-((x-g2Mean2)^2)/(2*g2SD2^2))),
          data=flowData,
          start=c(
            pop1N1=flowDataMeans01$pop1NumG1,
            pop1N2=flowDataMeans01$pop1NumG2,
            pop2N2=flowDataMeans01$pop2NumG2,
            A=0, B=0, C=0,
            A2=0, B2=0, C2=0
          ),
          nls.control(warnOnly=TRUE)
        )
      )
    }else{
      xpectr::suppress_mw(
        multiPopNLS <- nls(
          formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*
                           exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
            (A + B*x + C*(x^2))*
            (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                  (2*(g1SD*(x/g1Mean))^2)))+
            (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
            (A2 + B2*x + C2*(x^2))*
            (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                                                    (2*(g2SD2*(x/g2Mean2))^2)))+(pop2N1/(sqrt(2*pi)*g1SD2)*
                                                                                   exp(-((x-g1Mean2)^2)/(2*g1SD2^2)))+
            (pop2N2/(sqrt(2*pi)*g2SD2)*
               exp(-((x-g2Mean2)^2)/(2*g2SD2^2))),
          data=flowData,
          start=c(
            pop1N1=flowDataMeans01$pop1NumG1,
            pop1N2=flowDataMeans01$pop1NumG2,
            pop2N1=flowDataMeans01$pop2NumG1,
            pop2N2=flowDataMeans01$pop2NumG2,
            A=0, B=0, C=0,
            A2=0, B2=0, C2=0
          ),
          nls.control(warnOnly=TRUE)
        )
      )
    }
    
    if(saveGraph == TRUE){
      nlsGraph <- ggplot()+
        geom_line(data=flowData, aes(x=x, y=y, col ='Raw Data'), linewidth=1)+
        geom_line(aes(flowData$x,predict(multiPopNLS), col='NLS'))+
        labs(y="Counts", x=xVariable)+
        theme(legend.title=element_blank())
      
      plotDir <- "nlsGraphs"
      dir.create(file.path(dirname(flowDir), plotDir), showWarnings=FALSE)
      plotInitDir <- "nlsMultiple_graphs"
      dir.create(
        file.path(dirname(flowDir) ,plotDir, plotInitDir),
        showWarnings=FALSE
      )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
      
      png(
        paste0(plotOutFile, "/", flowNameDs[k], '.png'),
        width=600, height=400
      )
      print(nlsGraph)
      dev.off()
    }
    
    flowDataMeans01$residual2Pop <- summary(multiPopNLS)[["sigma"]]
    
    residualMultiDs <- rbind(
      residualMultiDs,
      flowDataMeans01
    )
    
  }
  
  return(residualMultiDs)
}

## popConfidence3Pop
.popConfidence3Pop = function(flowDir, ds, xVariable, saveGraph = TRUE){
  
  ##Removing NOTE 'no visible binding for global variable'
  G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
  G1_2<-G2_2<-G1Count_2<-G2Count_2<-NULL
  G1_3<-G2_3<-G1Count_3<-G2Count_3<-NULL
  
  ds1<-ds %>% dplyr::filter(
    !is.na(ds$G1_3)
  )
  
  modelData <- ds1 %>% dplyr::select(
    data,
    G1_1,
    G2_1,
    G1_2,
    G2_2,
    G1_3,
    G2_3,
    G1Count_1,
    G2Count_1,
    G1Count_2,
    G2Count_2,
    G1Count_3,
    G2Count_3
  )
  
  modelData02 <- modelData %>% dplyr::mutate(
    sdG1Count=G1Count_1*0.6,
    sdG2Count=G2Count_1*0.6,
    sdG1Count2=G1Count_2*0.6,
    sdG2Count2=G2Count_2*0.6,
    sdG1Count3=G1Count_3*0.6,
    sdG2Count3=G2Count_3*0.6
  )
  
  flowNameDs <- unique(modelData02$data)
  
  residualMultiDs <- data.frame(
    matrix(nrow=0, ncol=38)
  )
  
  colnames(residualMultiDs) <- c(
    "data",
    "G1_1",
    "G2_1",
    "G1_2",
    "G2_2", 
    "G1_3", 
    "G2_3", 
    "G1Count_1",
    "G2Count_1",
    "G1Count_2",
    "G2Count_2", 
    "G1Count_3",
    "G2Count_3", 
    "sdG1Count",
    "sdG2Count",
    "sdG1Count2",
    "sdG2Count2",
    "sdG1Count3",
    "sdG2Count3",
    "g1Mean",
    "g2Mean",
    "g1SD",
    "g2SD",
    "numG1",
    "numG2",
    "g1Mean2",
    "g2Mean2",
    "g1SD2",
    "g2SD2",
    "pop2NumG1",
    "pop2NumG2",
    "g1Mean3",
    "g2Mean3",
    "g1SD3",
    "g2SD3",
    "pop3NumG1",
    "pop3NumG2",
    "residual"
  )
  
  
  for(k in seq_len(length(flowNameDs))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    flowData <- .smoothData( flowName, xVariable, 5)
    flowDataMeans <- modelData02%>% dplyr::filter(
      data == flowNameDs[k]
    )
    
    if(0 %in% flowData$x){
      flowData <- flowData %>% dplyr::filter(x != 0)
    }
    
    midpoint1<- flowData[
      which(
        abs(
          flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2))
        ) == min(
          abs(flowData$x-mean(
            c(flowDataMeans$G1_1, flowDataMeans$G1_2)
          )
          )
        )
      ),
    ]
    
    if(nrow(midpoint1)>1){
      midpoint1 <- midpoint1[nrow(midpoint1), ]
    }
    
    midpoint2<- flowData[
      which(
        abs(
          flowData$x-mean(c(midpoint1$x, flowDataMeans$G1_3))
        ) == min(
          abs(flowData$x-mean(c(midpoint1$x, flowDataMeans$G1_3)))
        )
      ),
    ]
    
    if(nrow(midpoint2)>1){
      midpoint2 <- midpoint2[nrow(midpoint2), ]
    }
    
    midpoint3<- flowData[
      which(
        abs(
          flowData$x-mean(c(midpoint2$x, flowDataMeans$G2_2))
        ) == min(
          abs(flowData$x-mean(c(midpoint2$x, flowDataMeans$G2_2)))
        )
      ),
    ]
    
    if(nrow(midpoint3)>1){
      midpoint3 <- midpoint3[nrow(midpoint3), ]
    }
    
    midpoint4<-flowData[
      which(
        abs(
          flowData$x-mean(
            c(
              max(c(midpoint3$x,flowDataMeans$G2_2)),
              flowDataMeans$G2_3
            )
          )
        ) == min(
          abs(
            flowData$x-mean(
              c(
                max(c(midpoint3$x,flowDataMeans$G2_2)),
                flowDataMeans$G2_3
              )
            )
          )
        )
      ),
    ]
    
    if(nrow(midpoint4)>1){
      midpoint4 <- midpoint4[nrow(midpoint4), ]
    }
    
    #pop 1 G1 standard deviation 
    g1LeftFlowData <- flowData[
      seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
    g1RightFlowData <- flowData[
      (which(flowData$x == flowDataMeans$G1_1)+1): 
        (which(flowData$x == midpoint1$x)-1),
    ]
    
    g1LeftPeak1 <- g1LeftFlowData[
      which(
        abs(
          g1LeftFlowData$y-flowDataMeans$sdG1Count
        ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
      ), 
    ]
    
    if(nrow(g1LeftPeak1)>1){
      g1LeftPeak1 <- g1LeftPeak1[1, ]
    }
    
    g1RightPeak1 <- g1RightFlowData[
      which(
        abs(
          g1RightFlowData$y-flowDataMeans$sdG1Count
        ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
      ),
    ]
    
    if(nrow(g1RightPeak1)>1 ){
      g1RightPeak1 <- g1RightPeak1[1, ]
    }
    
    #pop 2 G1 standard deviation
    g1LeftFlowData2 <- flowData[
      which(flowData$x == midpoint1$x):
        (which(flowData$x == flowDataMeans$G1_2)-1),
    ]
    
    g1RightFlowData2 <- flowData[
      (which(flowData$x == flowDataMeans$G1_2)+1):
        (which(flowData$x == midpoint2$x)-1),
    ]
    
    pop2g1LeftPeak1 <- g1LeftFlowData2[
      which(
        abs(
          g1LeftFlowData2$y-flowDataMeans$sdG1Count2
        ) == min(abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2))
      ), 
    ]
    
    if(nrow(pop2g1LeftPeak1)>1){
      pop2g1LeftPeak1 <- pop2g1LeftPeak1[1, ]
    }
    
    pop2g1RightPeak1 <- g1RightFlowData2[
      which(
        abs(
          g1RightFlowData2$y-flowDataMeans$sdG1Count2
        ) == min(abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2))
      ),
    ]
    
    if(nrow(pop2g1RightPeak1)>1 ){
      pop2g1RightPeak1 <- pop2g1RightPeak1[1, ]
    }
    
    #pop 3 G1 standard deviation
    g1LeftFlowData3 <- flowData[
      which(flowData$x == midpoint2$x):
        (which(flowData$x == flowDataMeans$G1_3)-1),
    ]
    
    g1RightFlowData3 <- flowData[
      (which(flowData$x == flowDataMeans$G1_3)+1):
        (which(flowData$x == midpoint3$x)-1),
    ]
    
    #check for overlapping pops, only keep points right of peak
    g1RightFlowData3 <- g1RightFlowData3[
      (which(g1RightFlowData3$x > flowDataMeans$G1_3)),
    ]
    
    pop3g1LeftPeak1 <- g1LeftFlowData3[
      which(
        abs(
          g1LeftFlowData3$y-flowDataMeans$sdG1Count3
        ) == min(abs(g1LeftFlowData3$y-flowDataMeans$sdG1Count3))
      ), 
    ]
    
    if(nrow(pop3g1LeftPeak1)>1){
      pop3g1LeftPeak1 <- pop3g1LeftPeak1[1, ]
    }
    
    pop3g1RightPeak1 <- g1RightFlowData3[
      which(
        abs(
          g1RightFlowData3$y-flowDataMeans$sdG1Count3
        ) == min(abs(g1RightFlowData3$y-flowDataMeans$sdG1Count3))
      ),
    ]
    
    if(nrow(pop3g1RightPeak1)>1 ){
      pop3g1RightPeak1 <- pop3g1RightPeak1[1, ]
    }
    
    #pop 2 G2 standard deviation
    g2LeftFlowData2 <- flowData[
      which(flowData$x == midpoint3$x):
        (which(flowData$x == flowDataMeans$G2_2)-1),
    ]
    
    g2RightFlowData2 <- flowData[
      (which(flowData$x == flowDataMeans$G2_2)+1):
        (which(flowData$x == midpoint4$x)-1),
    ]
    
    pop2g2LeftPeak1 <- g2LeftFlowData2[
      which(
        abs(
          g2LeftFlowData2$y-flowDataMeans$sdG2Count2
        ) == min(abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2))
      ), 
    ]
    
    if(nrow(pop2g2LeftPeak1)>1){
      pop2g2LeftPeak1 <- pop2g2LeftPeak1[1, ]
    }
    
    pop2g2RightPeak1 <- g2RightFlowData2[
      which(
        abs(
          g2RightFlowData2$y-flowDataMeans$sdG2Count2
        ) == min(abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2))
      ),
    ]
    
    if(nrow(pop2g2RightPeak1)>1 ){
      pop2g2RightPeak1 <- pop2g2RightPeak1[1, ]
    }
    
    #pop 2 G3 standard deviation
    g2LeftFlowData3 <- flowData[
      which(flowData$x == midpoint4$x):
        (which(flowData$x == flowDataMeans$G2_3)-1),
    ]
    
    g2RightFlowData3 <- flowData[
      (which(flowData$x == flowDataMeans$G2_3)+1):
        nrow(flowData),
    ]
    
    pop3g2LeftPeak1 <- g2LeftFlowData3[
      which(
        abs(
          g2LeftFlowData3$y-flowDataMeans$sdG2Count3
        ) == min(abs(g2LeftFlowData3$y-flowDataMeans$sdG2Count3))
      ), 
    ]
    
    if(nrow(pop3g2LeftPeak1)>1){
      pop3g2LeftPeak1 <- pop3g2LeftPeak1[1, ]
    }
    
    pop3g2RightPeak1 <- g2RightFlowData3[
      which(
        abs(
          g2RightFlowData3$y-flowDataMeans$sdG2Count3
        ) == min(abs(g2RightFlowData3$y-flowDataMeans$sdG2Count3))
      ),
    ]
    
    if(nrow(pop3g2RightPeak1)>1 ){
      pop3g2RightPeak1 <- pop3g2RightPeak1[1, ]
    }
    
    
    flowDataMeans01 = flowDataMeans %>% dplyr::mutate(
      g1Mean=G1_1,
      g2Mean=G2_1,
      g1SD=mean(
        c(
          g1RightPeak1$x-flowDataMeans$G1_1,
          flowDataMeans$G1_1-g1LeftPeak1$x
        )
      ),
      g2SD=mean(
        c(
          pop3g1RightPeak1$x-flowDataMeans$G2_1, #added abs
          flowDataMeans$G2_1-pop3g1LeftPeak1$x
        )
      ),
      pop1NumG1=sum(
        flowData$y[
          c(
            which(flowData$x == g1LeftPeak1$x):
              which(flowData$x == pop3g1RightPeak1$x)
          )
        ]
      ),
      pop1NumG2=sum(
        flowData$y[
          c(
            which(flowData$x == pop3g1LeftPeak1$x):
              which(flowData$x == pop3g1RightPeak1$x)
          )
        ]
      ),
      g1Mean2=G1_2,
      g2Mean2=G2_2,
      g1SD2=mean(
        c(
          pop2g1RightPeak1$x-flowDataMeans$G1_2,
          flowDataMeans$G1_2-pop2g1LeftPeak1$x
        )
      ),
      g2SD2=mean(
        c(
          pop2g2RightPeak1$x-flowDataMeans$G2_2,
          flowDataMeans$G2_2-pop2g2LeftPeak1$x
        )
      ),
      pop2NumG1=sum(
        flowData$y[
          c(
            which(flowData$x == pop2g1LeftPeak1$x):
              which(flowData$x == pop2g1RightPeak1$x)
          )
        ]
      ),
      pop2NumG2=sum(
        flowData$y[
          c(
            which(flowData$x == pop2g2LeftPeak1$x):
              which(flowData$x == pop2g2RightPeak1$x)
          )
        ]
      ),
      
      g1Mean3=G1_3,
      g2Mean3=G2_3,
      g1SD3=mean(
        c(
          abs(pop3g1RightPeak1$x-flowDataMeans$G1_3), #added abs()
          flowDataMeans$G1_3-pop3g1LeftPeak1$x
        )
      ),
      g2SD3=mean(
        c(
          pop3g2RightPeak1$x-flowDataMeans$G2_3,
          flowDataMeans$G2_3-pop3g2LeftPeak1$x
        )
      ),
      pop3NumG1=sum(
        flowData$y[
          c(
            which(flowData$x == pop3g1LeftPeak1$x):
              which(flowData$x == pop3g2RightPeak1$x)
          )
        ]
      ),
      pop3NumG2=sum(
        flowData$y[
          c(
            which(flowData$x == pop3g2LeftPeak1$x):
              which(flowData$x == pop3g2RightPeak1$x)
          )
        ]
      )
    )
    
    g1Mean <- flowDataMeans01$g1Mean
    g2Mean <- flowDataMeans01$g2Mean
    g1SD <- flowDataMeans01$g1SD
    g2SD <- flowDataMeans01$g2SD
    g1Mean2 <- flowDataMeans01$g1Mean2
    g2Mean2 <- flowDataMeans01$g2Mean2
    g1SD2 <- flowDataMeans01$g1SD2
    g2SD2 <- flowDataMeans01$g2SD2
    g1Mean3 <- flowDataMeans01$g1Mean3
    g2Mean3 <- flowDataMeans01$g2Mean3
    g1SD3 <- flowDataMeans01$g1SD3
    g2SD3 <- flowDataMeans01$g2SD3
    
    xpectr::suppress_mw(
      multiPopNLS <- nls(
        formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*
                         exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
          (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
          (pop2N1/(sqrt(2*pi)*g1SD2)*exp(-((x-g1Mean2)^2)/(2*g1SD2^2)))+
          (pop2N2/(sqrt(2*pi)*g2SD2)* exp(-((x-g2Mean2)^2)/(2*g2SD2^2)))+
          (pop3N2/(sqrt(2*pi)*g2SD3)* exp(-((x-g2Mean3)^2)/(2*g2SD3^2)))+
          (A + B*x + C*(x^2))*
          (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                (2*(g1SD*(x/g1Mean))^2)))+
          (A2 + B2*x + C2*(x^2))*
          (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                                                  (2*(g2SD2*(x/g2Mean2))^2)))+
          (A3 + B3*x + C3*(x^2))*
          (1/(sqrt(2*pi)*g2SD3*(x/g2Mean3))*exp(-((x-g2Mean3)^2)/
                                                  (2*(g2SD3*(x/g2Mean3))^2)))+
          (A4 + B4*x + C4*(x^2))*
          (1/(sqrt(2*pi)*g2SD*(x/g2Mean))*exp(-((x-g2Mean)^2)/
                                                (2*(g2SD*(x/g2Mean))^2))),
        data=flowData,
        start=c(
          pop1N1=flowDataMeans01$pop1NumG1,
          pop1N2=flowDataMeans01$pop1NumG2,
          pop2N1=flowDataMeans01$pop2NumG1,
          pop2N2=flowDataMeans01$pop2NumG2,
          pop3N2=flowDataMeans01$pop3NumG2,
          A=0, B=0, C=0,
          A2=0, B2=0, C2=0,
          A3=0, B3=0, C3=0,
          A4=0, B4=0, C4=0
        ),
        nls.control(warnOnly=TRUE)
      )
    )
    
    if(saveGraph == TRUE){
      nlsGraph <- ggplot()+
        geom_line(data=flowData, aes(x=x, y=y, col ='Raw Data'), linewidth=1)+
        geom_line(aes(flowData$x,predict(multiPopNLS), col='NLS'))+
        labs(y="Counts", x=xVariable)+
        theme(legend.title=element_blank())
      
      plotDir <- "nlsGraphs"
      dir.create(
        file.path(dirname(flowDir), plotDir),
        showWarnings=FALSE
      )
      plotInitDir <- "nlsMultiple_graphs"
      dir.create(
        file.path(
          dirname(flowDir) ,plotDir, plotInitDir),
        showWarnings=FALSE
      )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
      
      png(
        paste0(plotOutFile, "/", flowNameDs[k], '.png'),
        width=600, height=400
      )
      print(nlsGraph)
      dev.off()
    }
    
    flowDataMeans01$residual3Pop <- summary(multiPopNLS)[["sigma"]]
    
    residualMultiDs <- rbind(
      residualMultiDs,
      flowDataMeans01
    )
  }
  
  return(residualMultiDs)
}
