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
#' @param maxDoubletHeight The maximum height a doublet can be. If left as NA
#'  the algorithm will find a value based on the other peaks
#'
#' @export
#'
#' @examples
#' \dontrun{
#' peakAlgorithm2(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  flaggedData_ = flaggedData ,
#'  xVariable = "FITC-A",
#'  usedCellsThreshold = 86,
#'  maxDoubletHeight = 50
#'  )
#'  }
#'


peakAlgorithm2 = function(
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

  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    if(!xVariable %in% flowName@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    flowData <- smoothData( flowName, xVariable, 5)
    
    logFlow[1, ] <- c(algorithmNum, flowName@description[["GUID"]], 0)
    
    .GlobalEnv$logDs <- rbind(
      .GlobalEnv$logDs,
      logFlow
    )
    
    localPeaks <- detect_localmaxima(flowData$y, 3)
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
    possiblePeaks3 <- findClusters(possiblePeaks2, 40, xVarMax)
    
    possiblePeaks4 <- findPairs(possiblePeaks3, possiblePeaks3, 1.75, 2.2)
    
    possiblePeaks4 <- possiblePeaks4 %>%
      tidyr::drop_na() %>%
      dplyr::filter(
        possiblePairY > quantile(flowData$y)[4]
      )
    possiblePeaks5 <- possiblePeaks4[!duplicated(possiblePeaks4$possiblePairY), ]
    
    if(nrow(possiblePeaks5) > 1){
      tempDs<-possiblePeaks5
      tempDs2<-c()
      for(i in 1:nrow(tempDs)){
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
    possiblePeaks6 <- doubletCheck(
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
      for(j in 1:nrow(peakRow)){
        
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