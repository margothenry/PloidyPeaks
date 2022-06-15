#' peakAlgorithm2
#'
#' The second branching point is in peakAlgorithm2, which seeks to identify all
#' possible sub populations. The algorithm looks at the proportion of cells used
#' in the subpopulations (CellProp). In other words, once the populations have been
#' identified for a given flow frame, peakAlgorithm2 will find the cells associated with the G1/G2 pairs.
#'
#' @param flowDir The directory of the gated .fcs data
#' @param xVariable The fluorescence channel on the x axis
#' @param flaggedData_ List of names of the flow frames that got flagged by PeakAlgorithm1
#' @param usedCellsThreshold Threshold for classifying multiple populations
#'
#' @export
#'
#' @examples
#' peakAlgorithm2(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  flaggedData_ = flaggedData ,
#'  xVariable = "FITC-A",
#'  usedCellsThreshold = 86
#'  )
#'


peakAlgorithm2 = function(flowDir, flaggedData_, xVariable, usedCellsThreshold = 86){
  flowNameDs <- flaggedData_$data
  finishedData <- c()
  flaggedData <- c()

  logFlow <- data.frame(
    matrix(nrow = 0, ncol = 3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")
  algorithmNum <- 2

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

    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > quantile(flowData$y)[3]+5),
      ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/10), ]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > quantile(flowData$y)[3]+5),
        ]
    }
    possiblePeaks3 <- findClusters(possiblePeaks2, 40, xVarMax)

    possiblePeaks4 <- findPairs(possiblePeaks3, 1.7, 2.3)

    possiblePeaks4 <- possiblePeaks4 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        y > quantile(flowData$y)[4] & possiblePairY > quantile(flowData$y)[3]+7
        )
    possiblePeaks5 <- possiblePeaks4[!duplicated(possiblePeaks4$possiblePairY),]

    if(nrow(possiblePeaks5) == 1){
      rangeLength <- nchar(format(xVarMax, scientific=F))
      multiplier <- 10^(rangeLength-3)
      possiblePeaks6 <- doubletCheck(
        possiblePeaks5,
        possiblePeaks,
        10*multiplier,
        15*multiplier
        )
    }else{
      possiblePeaks6 <- possiblePeaks5 %>%
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


    if(nrow(possiblePeaks6) != 0){
      #function start
      if(!is.na(possiblePeaks6$g1G2Doublet[1])){
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
          peaks = c(
            possiblePeaks6$x,
            possiblePeaks6$possiblePairX
            )
          )
      }

      peakRow <- peakRow %>%
        dplyr::distinct() %>% tidyr::drop_na()
      cellsUsed <- c()
      for(j in 1:nrow(peakRow)){

        xPeak <- which(flowData$x == peakRow$peaks[j])
        if(length(xPeak)>1){
          xPeak <- xPeak[1]
        }

        flowData$keep <- FALSE
        flowData$keep[xPeak] <- TRUE

        #Forwards
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

        #Backwards
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

        xEpsilonRight <- which(flowData$x==intPeak$x[3])
        xEpsilonLeft <- which(flowData$x==intPeak$x[1])
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
      propCellsUsed <- round((cellsUsed01/totalCellCount)*100,2)
      possiblePeaks6$propCellsUsed <- propCellsUsed

    }else{
      possiblePeaks6 <- data.frame(
        x = NA,
        y = NA,
        cluster = NA,
        LL = NA,
        UL = NA,
        possiblePairX= NA,
        possiblePairY= NA,
        data = flowName@description[["GUID"]],
        propCellsUsed = NA
      )
    }

    if(
      possiblePeaks6$propCellsUsed[1] >= usedCellsThreshold &
      !is.na(possiblePeaks6$propCellsUsed[1])
      ){

      possiblePeaks7 <- possiblePeaks6 %>% dplyr::mutate(
        data = flowName@description[["GUID"]]
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
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs),]$Success <- 1
  }

  returnedList <- list(flaggedData,finishedData)

  return(returnedList)

}




