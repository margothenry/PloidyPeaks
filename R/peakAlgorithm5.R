#' peakAlgorithm5
#'
#' The last branching point is in peakAlgorithm5, which seeks to identify all
#' possible sub populations that got flagged from PeakAlgorithm4. The algorithm
#' looks at the proportion of cells used in the subpopulations (CellProp).
#' In other words, once the populations have been identified for a given flow frame,
#' peakAlgorithm5 will find the cells associated with the G1/G2 pairs. The difference
#' from PeakAlgorithm4 is that peakAlgorithm5 works with more granular data, i.e, smooth level = 3, as well try and find a G2 pair in a (1.5*G1, 2.5*G1)
#'
#' @param flowDir The directory of the gated .fcs data
#' @param xVariable The fluorescence channel on the x axis
#' @param flaggedData_ List of names of the flow frames that got flagged by PeakAlgorithm3
#' @param appendData The data set of flow frames that are finished being analyzed
#' from the previous peak algorithms.
#'
#' @export
#'
#' @examples
#' peakAlgorithm5(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  flaggedData_ = flagged_data ,
#'  xVariable = "FITC-A",
#'  appendData = analyzed_ds
#'  )
#'

peakAlgorithm5 = function(flowDir, flaggedData_, xVariable, appendData){
  naDs <- flaggedData_ %>% dplyr::filter(is.na(x))
  flowNameDs <- unique(naDs$data)

  logFlow <- data.frame(
    matrix(nrow = 0, ncol = 3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")

  algorithmNum <- 5

  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir,"/",flowNameDs[k]), transformation=FALSE
      )
    flowData <- smoothData( flowName, xVariable, 5)

    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)

    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )

    localPeaks <- detect_localmaxima(flowData$y,3)
    possiblePeaks <- flowData[localPeaks,]

    possiblePeaks2 = possiblePeaks[
      which(possiblePeaks$y > quantile(flowData$y)[3]),
      ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/9.5), ]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > quantile(flowData$y)[3]),
        ]
    }
    possiblePeaks3 <- findClusters(possiblePeaks2, 30, xVarMax)

    possiblePeaks4 <- possiblePeaks3[
      order(possiblePeaks3$x, decreasing = FALSE),
      ]
    possiblePeaks4 <- possiblePeaks4[
      which(possiblePeaks4$x < quantile(flowData$x)[4]),
      ]
    peaksFix <- possiblePeaks4
    firstPeak <- possiblePeaks4[1,1:2]

    if(nrow(possiblePeaks4) == 1){

      initialPeak <- firstPeak
      #the peak is G1 missing G2
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
      #the peak is G2 missing G1
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
    possiblePeaks5 <- findPairs(peaksFix, 1.4, 2.5)

    possiblePeaks5 <- possiblePeaks5 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        y > quantile(flowData$y)[3] & possiblePairY > quantile(flowData$y)[3]
        )
    possiblePeaks6 <- possiblePeaks5[!duplicated(possiblePeaks5$possiblePairY),]

    if(nrow(possiblePeaks6) > 1){
      possiblePeaks7 <- findClusters(possiblePeaks6, 30, xVarMax)
    }else{
      possiblePeaks7 <- possiblePeaks6 %>% dplyr::mutate(
        cluster = 1,
        distToNext = 0
      )
    }

    if(nrow(possiblePeaks7) == 0){
      maxPeaksFix <- peaksFix[which(max(peaksFix$y) == peaksFix$y),]
      possiblePeaks7 <- data.frame(
        x = maxPeaksFix$x,
        y = maxPeaksFix$y,
        cluster = 1,
        distToNext = 0,
        LL = NA,
        UL = NA,
        possiblePairX = 0,
        possiblePairY = 0
      )

    }

    if(nrow(possiblePeaks7) == 1){
      rangeLength <- nchar(format(xVarMax, scientific=F))
      multiplier <- 10^(rangeLength-3)
      possiblePeaks8 <- doubletCheck(
        possiblePeaks7,
        possiblePeaks,
        15*multiplier,
        15*multiplier
        )
    }else{
      possiblePeaks8 <- possiblePeaks7 %>%
        dplyr::mutate(
          g3LL = NA,
          g3UL = NA,
          g4LL = NA,
          g4UL = NA,
          g1G2Doublet = NA,
          g1G2DoubletCount = NA,
          g2G2Doublet = NA,
          g2G2DoubletCount = NA
        )
    }

    possiblePeaks9 <- possiblePeaks8 %>% dplyr::mutate(
      data = flowName@description[["GUID"]],
      messy = 1
    )

    appendData <- rbind(
      appendData,
      possiblePeaks9
    )

    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs),]$Success <- 1
  }

  returnedList <- list(appendData)

  return(returnedList)
}

