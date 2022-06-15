#' peakAlgorithm1
#'
#' peakAlgorithm1 is the first branching point to identify samples with single
#' populations (which contain only two peaks and can be easily analysed) and to
#' flag everything else. The algorithm will find a single G1 and G2 pairing.
#' With that population, a boundary is set to the right of the G2 peak and the
#' ratio of cells in that area are calculated. If this ratio is within a certain
#' threshold, the population is marked as a single population. If the ratio exceeds
#' the threshold, the population is flagged to pass to the next algorithm, PeakAlgorithm2.
#'
#' @param flowDir The directory of the gated .fcs data
#' @param xVariable The fluorescence channel on the x axis
#' @param flowSet List of names of the flow frames in the flow set
#' @param singleThreshold threshold for classifying single populations
#'
#' @export
#'
#' @examples
#' peakAlgorithm1(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  flowSet = flow_names_ds ,
#'  xVariable = "FITC-A"
#' )
#'

peakAlgorithm1 = function(flowDir, flowSet, xVariable, singleThreshold = 8){

  singleData <- c()
  flaggedData <- c()

  #Adding flow frame to log
  logFlow <- data.frame(
    matrix(nrow = 0, ncol = 3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  algorithmNum <- 1

  #Looping through each flow frame
  for(k in 1:length(flowSet)){
    #Reading in and smoothing data
    flowName <- flowCore::read.FCS( paste0(flowDir,"/",flowSet[k]), transformation=FALSE)
    flowData <- smoothData( flowName, xVariable, 13)

    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)

    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )

    #Finding local peaks
    localPeaks <- detect_localmaxima(flowData$y,5)
    possiblePeaks <- flowData[localPeaks,]

    #Removing the peaks that are identified at the base of the histogram
    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > quantile(flowData$y)[3]+5),
      ]
    xVarMax <- max(flowData$x)
    #Removing the peaks that are identified at the extreme left side of the
    #Histogram that could be caused by debris/improper gating
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/9.5), ]

    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > quantile(flowData$y)[3]+5),
        ]
    }

    possiblePeaks3 <- findClusters(possiblePeaks2, 40, xVarMax)

    #Ordering the peaks
    possiblePeaks4 <- possiblePeaks3[order(possiblePeaks3$x, decreasing=FALSE),]
    #Removing the peaks that are identified at the extreme right side of the
    #Histogram
    possiblePeaks4 <- possiblePeaks4[
      which(possiblePeaks4$x < quantile(flowData$x)[4]),
      ]
    peaksFix <- possiblePeaks4
    #Identifying the first peak in the dataset
    firstPeak <- possiblePeaks4[1,1:2]
    #If there is only one peak identified, try to identify the pair for this peak
    if(nrow(possiblePeaks4) == 1){

      initialPeak <- firstPeak
      #The peak is G1 and G2 is missing
      g2PeakRadiusUL <- initialPeak$x*2.2 + 1
      g2PeakRadiusLL <- initialPeak$x*1.75 - 1
      g2ToTestDs2 <- which(
        g2PeakRadiusLL < possiblePeaks$x & possiblePeaks$x < g2PeakRadiusUL &
        possiblePeaks$x < quantile(flowData$x)[4]
        )
      g2ToTestDs3 <- possiblePeaks[g2ToTestDs2,]
      if(nrow(g2ToTestDs3) == 0){
        g2ToTestDs3 <- data.frame(
          x = 0,
          y = 0
        )
      }
      #The peak is G2 and G1 is missing
      g1PeakRadiusLL <- initialPeak$x/2.2 + 1
      g1PeakRadiusUL <- initialPeak$x/1.75 - 1
      g1ToTestDs2 <- which(
        g1PeakRadiusLL < possiblePeaks$x & possiblePeaks$x < g1PeakRadiusUL
        )
      g1ToTestDs3 <- possiblePeaks[g1ToTestDs2,]
      if(nrow(g1ToTestDs3) == 0){
        g1ToTestDs3 <- data.frame(
          x = 0,
          y = 0
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

    #checking if we have the proper G1 peak
    peaksFix <- peaksFix[order(peaksFix$x, decreasing = FALSE),]
    firstPeak <- peaksFix[1,1:2]
    if(firstPeak$x > xVarMax/3.1 | firstPeak$x < xVarMax/6 ){

      initialPeak <- firstPeak
      otherPeak <- peaksFix[-1,1:2]

      #the peak is G2 missing G1
      g1PeakRadiusLL <- initialPeak$x/2.3 - 1
      g1PeakRadiusUL <- initialPeak$x/1.7 + 1
      g1ToTestDs2 <- which(
        g1PeakRadiusLL <= possiblePeaks$x & possiblePeaks$x <= g1PeakRadiusUL
        )
      g1ToTestDs3 <- possiblePeaks[g1ToTestDs2,]

      if(nrow(g1ToTestDs3) == 0){
        g1ToTestDs2 <- which(
          g1PeakRadiusLL <= flowData$x & flowData$x <= g1PeakRadiusUL
          )
        g1ToTestDs3 <- flowData[g1ToTestDs2,]
      }
      g1ToTestDs3 <- g1ToTestDs3[order(g1ToTestDs3$y, decreasing = TRUE),]
      g1ToTestDs3 <- g1ToTestDs3[which(g1ToTestDs3$y > quantile(flowData$y)[2]),]
      peak2 <- g1ToTestDs3[1,]
      peak1 <- initialPeak
      peaksFix <- rbind(
        peak1,
        peak2,
        otherPeak
      )
    }

    #Finding the are to the right of the G2 peak
    peaksFix <- peaksFix[order(peaksFix$x, decreasing = FALSE),]
    peaksFix <- peaksFix %>% dplyr::distinct()

    #Finding the distance between the G1 and G2 peak
    #
    epsilon <- abs((peaksFix$x[2] - peaksFix$x[1])/2)
    epsilonLeft <- peaksFix$x[2] - epsilon
    epsilonRight <- peaksFix$x[2] + epsilon

    xEpsilon <- which(
      abs(flowData$x-epsilonRight)==min(abs(flowData$x-epsilonRight))
      )
    if(length(xEpsilon)>1){
      xEpsilon <- xEpsilon[1]
    }
    postEpsilon <- flowData$y[c(xEpsilon:length(flowData$x))]
    preEpsilon <- flowData$y[c(1:xEpsilon)]
    cellCountPreEpsilon <- sum(preEpsilon)
    cellCountPostEpsilon <- sum(postEpsilon)
    totalCellCount <- sum(flowData$y)
    #Proportion of cells that are to the right of the G2 peak
    detect <- round((cellCountPostEpsilon/totalCellCount)*100,4)

    if(detect <= singleThreshold){

      if(peaksFix$x[1] == 0){
        peaksFix[1,] <- peaksFix[2,]
        peaksFix$x[2] <- 0
        peaksFix$y[2] <- 0
      }

      singleDs <- data.frame(
        data = flowName@description[["GUID"]],
        x = peaksFix$x[1],
        y = peaksFix$y[1],
        possiblePairX = peaksFix$x[2],
        possiblePairY = peaksFix$y[2]
      )

      rangeLength <- nchar(format(xVarMax, scientific=F))
      multiplier <- 10^(rangeLength-3)
      singleDs2 <- doubletCheck(
        singleDs,
        possiblePeaks,
        10*multiplier,
        15*multiplier
        )

      singleData <- rbind(
        singleData,
        singleDs2
      )
      singleDataUpdated <- updatedMeans(singleData, flowDir, xVariable)
    }else{
      flaggedData <- rbind(
        flaggedData,
        flowName@description[["GUID"]]
      )
    }

    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs),]$Success <- 1

  }

  returnedList <- list(flaggedData,singleDataUpdated)

  return(returnedList)
}

