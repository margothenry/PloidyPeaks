#' DJFMeans
#'
#' A function that will update the peaks identified for the throughout the
#' analysis for the Dean-Jett-Fox model. The model is sensitive to noisy data
#' and so we smooth out the data to run the model. DJFMeans will find the peaks
#' associated with the data at the smoothing level indicated.
#'
#' @param ds The dataset 
#' @param flowDir The directory of the gated .fcs data
#' @param xVariable The fluorescence channel on the x axis
#' @param smoothLevel The level of smoothing applied to the data
#' @export
#'
#' @examples
#'DJFMeans(
#' ds = flowData,
#' flowDir = "FlowData/T10_FLC/gated_data",
#' xVariable = "FITC-A",
#' smoothLevel = 10
#' )


DJFMeans = function(ds, flowDir, xVariable, smoothLevel){
  ##Removing NOTE 'no visible binding for global variable'
  G1_1<-G1Count_1<-G2_1<-G2Count_1<-NULL
  G1Old<-G1CountOld<-G2Old<-G2CountOld<-NULL

  
  flowNameDs <- ds$data
  allPopUpdated <- c()
  
  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    
    flowData <- smoothData( flowName, xVariable, smoothLevel)
    allPop01 <- ds %>% dplyr::filter(
      data == flowNameDs[k]
    ) %>% dplyr::rename(
      G1Old=G1_1,
      G1CountOld=G1Count_1,
      G2Old=G2_1,
      G2CountOld=G2Count_1
    )
    
    midPoint <- flowData[
      which(
        abs(
          flowData$x-mean(c(allPop01$G1Old, allPop01$G2Old))
        ) == min(
          abs(flowData$x-mean(c(allPop01$G1Old, allPop01$G2Old)))
        )
      ),
    ]
    
    if(nrow(midPoint)>1){
      midPoint <- midPoint[nrow(midPoint), ]
    }
    
    xPeak1 <- which(flowData$x == allPop01$G1Old)
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
    
    if(allPop01$G2Old != 0){
      xPeak2 <- which(flowData$x == allPop01$G2Old)
      if(length(xPeak2)>1){
        xPeak2 <- xPeak2[1]
      }
      midPointX = which(midPoint$x == flowData$x)
      if(xPeak2+10 > nrow(flowData) & xPeak2-10 < midPointX){
        peakRowSmoothedRange <- flowData[seq(midPointX, nrow(flowData), 1), ]
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
            xPeak2Smoothed> midPointX
          )
        ]
        xPeak2Smoothed <- xPeak2Smoothed[1]
      }
      
      allPop02 <- allPop01 %>% dplyr::mutate(
        G1_1=flowData[xPeak1Smoothed, ]$x,
        G1Count_1=flowData[xPeak1Smoothed, ]$y,
        G2_1=flowData[xPeak2Smoothed, ]$x,
        G2Count_1=flowData[xPeak2Smoothed, ]$y
      ) %>% dplyr::select(
        data,
        G1_1,
        G2_1,
        G1Count_1,
        G2Count_1
      )
    }else{
      
      allPop02 <- allPop01 %>% dplyr::mutate(
        G1_1=flowData[xPeak1Smoothed, ]$x,
        G1Count_1=flowData[xPeak1Smoothed, ]$y,
        G2_1=G2Old,
        G2Count_1=G2CountOld
      ) %>% dplyr::select(
        data,
        G1_1,
        G2_1,
        G1Count_1,
        G2Count_1
      )
    }
    
    
    allPopUpdated <- rbind(
      allPopUpdated,
      allPop02
    )
    
  }
  returnedList <- data.frame(allPopUpdated)
}




