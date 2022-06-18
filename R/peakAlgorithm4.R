#' peakAlgorithm4
#'
#' The forth branching point is in peakAlgorithm4, which seeks to identify all
#' possible sub populations that got flagged from PeakAlgorithm3. The algorithm
#' looks at the proportion of cells used in the subpopulations (CellProp).
#' In other words, once the populations have been identified for a given flow frame,
#' peakAlgorithm4 will find the cells associated with the G1/G2 pairs. The difference
#' from PeakAlgorithm3 is that peakAlgorithm4 works with more granular data, i.e, smooth level = 4.
#' Additionally, those that get flagged from peakAlgorithm4 will get a 'messy' indicator.
#' Letting the user know that they may want to manually look those who got flagged as messy.
#'
#' @param flowDir The directory of the gated .fcs data
#' @param x_variable The fluorescence channel on the x axis
#' @param flagged_data_ds List of names of the flow frames that got flagged by PeakAlgorithm3
#' @param append_data The data set of flow frames that are finished being analyzed
#' from the previous peak algorithms. The flow frames that are done after peakAlgorithm4
#' will append to this data set.
#' @param usedCellsThreshold Threshold for classifying multiple populations
#'
#' @export
#'
#' @examples
#' peakAlgorithm4(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  flagged_data_ds = flagged_data ,
#'  x_variable = "FITC-A",
#'  append_data = analyzed_ds,
#'  usedCellsThreshold = 86
#'  )
#'


peakAlgorithm4 = function(flowDir, flaggedData_, xVariable, appendData, usedCellsThreshold = 86){
  flowNameDs <- flaggedData_$data
  flaggedData <- c()

  logFlow <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(logFlow) <- c("Algorithm", "Data", "Success")

  algorithmNum <- 4

  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
      )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    
    flowData <- smoothData( flowName, xVariable, 4)

    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)

    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )

    localPeaks <- detect_localmaxima(flowData$y, 3)
    possiblePeaks <- flowData[localPeaks, ]

    possiblePeaks2 <- possiblePeaks[
      which(possiblePeaks$y > quantile(flowData$y)[3]+5),
      ]
    xVarMax <- max(flowData$x)
    possiblePeaks2 <- possiblePeaks2[which(possiblePeaks2$x > xVarMax/9.5), ]
    if(nrow(possiblePeaks2) == 0){
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > quantile(flowData$y)[3]+5),
        ]
    }
    possiblePeaks3 <- findClusters(possiblePeaks2, 40, xVarMax)

    possiblePeaks4 <- findPairs(possiblePeaks3, 1.7, 2.3)

    possiblePeaks5 = possiblePeaks4 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        y > quantile(flowData$y)[4]-15 & possiblePairY > quantile(flowData$y)[3]
        )
    possiblePeaks6 <- possiblePeaks5[!duplicated(possiblePeaks5$possiblePairY), ]


    if(nrow(possiblePeaks6) > 1){
      possiblePeaks7 <- findClusters(possiblePeaks6, 40, xVarMax)
    }else{
      possiblePeaks7 <- possiblePeaks6
    }

    if(nrow(possiblePeaks7) == 1){
      rangeLength <- nchar(format(xVarMax, scientific=FALSE))
      multiplier <- 10^(rangeLength-3)
      possiblePeaks8 <- doubletCheck(
        possiblePeaks7,
        possiblePeaks,
        10*multiplier,
        15*multiplier
        )
    }else{
      possiblePeaks8 <- possiblePeaks7 %>%
        dplyr::mutate(
          g3LL=NA,
          g3UL=NA,
          g4LL=NA,
          g4UL=NA,
          g1G2Doublet=NA,
          g1G2DoubletCount=NA,
          g2G2Doublet=NA,
          g2G2DoubletCount=NA
        )
    }

    if(nrow(possiblePeaks8) != 0){
      if(!is.na(possiblePeaks8$g1G2Doublet[1])){
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
      smoothedData <- smoothData(flowName, xVariable, 15)

      for(j in 1:nrow(peakRow)){
        xPeak <- which(smoothedData$x == peakRow$peaks[j])
        if(length(xPeak)>1){
          xPeak <-xPeak[1]
        }

        if(xPeak+10 > nrow(smoothedData) & xPeak-10 < 1){
          peakRowSmoothedRange <- smoothedData[seq(1, nrow(smoothedData), 1), ]
        }else if(xPeak+10 > nrow(smoothedData) & xPeak-10 >= 1){
          peakRowSmoothedRange <- smoothedData[
            seq(xPeak-10, nrow(smoothedData), 1),
            ]
        }else if(xPeak+10 <= nrow(smoothedData) & xPeak-10 < 1){
          peakRowSmoothedRange <- smoothedData[seq(1, xPeak+10, 1), ]
        }else{
          peakRowSmoothedRange <- smoothedData[seq(xPeak-10, xPeak+10, 1), ]
        }

        xPeakSmoothed <- which(max(peakRowSmoothedRange$y) == smoothedData$y)
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

      possiblePeaks9 <- possiblePeaks8 %>% dplyr::mutate(
        data=flowName@description[["GUID"]],
        messy=1
      )

      flaggedData <- rbind(
        flaggedData,
        possiblePeaks9
      )


    }
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
  }

  returnedList <- list(flaggedData, appendData)

  return(returnedList)

}



