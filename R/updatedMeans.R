#' updatedMeans
#'
#' A function that will update the peaks identified for the single populations
#' to find their true peaks. The single populations are smoothed out by a factor
#' of 13 and so updatedMeans will find the peaks associated with the data at
#' a smoothing factor of 5.
#'
#' @param ds The dataset that has a smoothing factor of 13.
#' @param flowDir The directory of the gated .fcs data
#' @param xVariable The fluorescence channel on the x axis
#' @export
#'
#' @examples
#'updatedMeans(
#' ds = singleData,
#' flowDir = "FlowData/T10_FLC/gated_data",
#' xVariable = "FITC-A"
#' )


updatedMeans = function(ds, flowDir, xVariable){

  flowNameDs <- ds$data
  singlePopUpdated <- c()

  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir,"/",flowNameDs[k]), transformation=FALSE
    )

    flowData <- smoothData( flowName, xVariable, 5)
    singlePop01 <- ds %>% dplyr::filter(
      data == flowNameDs[k]
    ) %>% rename(
      xOld = x,
      yOld = y,
      possiblePairXOld = possiblePairX,
      possiblePairYOld = possiblePairY
    )

    xPeak1 <- which(flowData$x == singlePop01$xOld)
    if(length(xPeak1)>1){
      xPeak1 <- xPeak1[1]
    }

    if(xPeak1+10 > nrow(flowData) & xPeak1-10 < 1){
      peakRowSmoothedRange <- flowData[seq(1, nrow(flowData),1), ]
    }else if(xPeak1+10 > nrow(flowData) & xPeak1-10 >= 1){
      peakRowSmoothedRange <- flowData[
        seq(xPeak1-10, nrow(flowData),1),
      ]
    }else if(xPeak1+10 <= nrow(flowData) & xPeak1-10 < 1){
      peakRowSmoothedRange <- flowData[seq(1, xPeak1+10,1), ]
    }else{
      peakRowSmoothedRange <- flowData[seq(xPeak1-10, xPeak1+10,1), ]
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

      if(xPeak2+10 > nrow(flowData) & xPeak2-10 < 1){
        peakRowSmoothedRange <- flowData[seq(xPeak1, nrow(flowData),1), ]
      }else if(xPeak2+10 > nrow(flowData) & xPeak2-10 >= 1){
        peakRowSmoothedRange <- flowData[
          seq(xPeak2-10, nrow(flowData),1),
        ]
      }else if(xPeak2+10 <= nrow(flowData) & xPeak2-10 < 1){
        peakRowSmoothedRange <- flowData[seq(xPeak1, xPeak2+10,1), ]
      }else{
        peakRowSmoothedRange <- flowData[seq(xPeak2-10, xPeak2+10,1), ]
      }

      xPeak2Smoothed <- which(max(peakRowSmoothedRange$y) == flowData$y)
      if(length(xPeak2Smoothed) > 1){
        xPeak2Smoothed =  xPeak2Smoothed[which(xPeak1Smoothed<xPeak2Smoothed)]
        xPeak2Smoothed <- xPeak2Smoothed[1]
      }

      singlePop02 <- singlePop01 %>% dplyr::mutate(
        x = flowData[xPeak1Smoothed,]$x,
        y = flowData[xPeak1Smoothed,]$y,
        possiblePairX = flowData[xPeak2Smoothed,]$x,
        possiblePairY = flowData[xPeak2Smoothed,]$y
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
        x = flowData[xPeak1Smoothed,]$x,
        y = flowData[xPeak1Smoothed,]$y,
        possiblePairX = possiblePairXOld,
        possiblePairY = possiblePairYOld
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




