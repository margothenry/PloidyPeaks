#' peakCorrection
#'
#' peakCorrection is a function created for the user to rerun the analysis on 
#' a sample with a certain number of subpopulations.
#'
#' @param xVariable The fluorescence channel on the x axis
#' @param flowDir The directory of the gated .fcs data
#' @param sampleName the name of the sample you want to analyze
#' @param numSubPop the number of subpopulations you want the algorithm to identify
#' @export
#'
#' @examples
#' \dontrun{
#' peakCorrection(
#'  xVariable = "FITC-A",
#'  flowDir = NA,
#'  sampleName = "A01-A01",
#'  numSubPop = 2
#'  )
#'  }
#'
peakCorrection = function(
  xVariable,
  flowDir = NA,
  sampleName,
  numSubPop = 2
){
  
  .GlobalEnv$logDs <- data.frame(
    matrix(nrow=0, ncol=3)
  )
  colnames(.GlobalEnv$logDs) <- c("Algorithm", "Data", "Success")
  
  algorithmNum <- 1
  
  .GlobalEnv$logDs[1, ] <- c(algorithmNum, sampleName, 0)
  
  ##Directory: If user do not hard codes their directory in the function
  ##a window will open and ask the user to pick their directory.
  ##This should be where the data is located
  if(is.na(flowDir)){
    getwd()
    flowDir <- tclvalue(tkchooseDirectory())
  }
  
  flowName<-flowCore::read.FCS(
    paste0(flowDir, "/", sampleName), transformation=FALSE
  )
  
  # Warning Check for X variable
  if(!xVariable %in% flowName@parameters@data$name){
    stop("Your X variable is not in the dataset")
  }
  
  # Warning Check for Sample Name
  if(!sampleName %in% list.files(flowDir)){
    stop("Your Sample is not in the list of flow frames in the directory provided")
  }
  
  # Smoothing Data
  flowData<-smoothData( flowName, xVariable, 5)
  
  # Applying peak algorithm
  localPeaks<-detect_localmaxima(flowData$y, 3)
  possiblePeaks<-flowData[localPeaks, ]
  
  possiblePeaks2<-possiblePeaks[
    which(possiblePeaks$y > max(possiblePeaks$y)/20),
  ]
  
  xVarMax <- max(flowData$x)
  possiblePeaks3<-findClusters(possiblePeaks2, 40, xVarMax)
  
  possiblePeaks4<-findPairs(possiblePeaks3, possiblePeaks3, 1.75, 2.2)
  possiblePeaks4<-possiblePeaks4 %>%
    tidyr::drop_na() 
  possiblePeaks5<-possiblePeaks4[!duplicated(possiblePeaks4$possiblePairY), ]
  
  # If sub populations < numSubPop, find missing population
  if(nrow(possiblePeaks5) < numSubPop){
    range_increase<-0
    range_increase_val = 0.2
    while(!nrow(possiblePeaks5) == numSubPop){
      range_increase = range_increase + range_increase_val
      
      flowData <- smoothData( flowName, xVariable, 4)
      
      
      localPeaks <- detect_localmaxima(flowData$y, 3)
      possiblePeaks <- flowData[localPeaks, ]
      
      possiblePeaks2 <- possiblePeaks[
        which(possiblePeaks$y > max(possiblePeaks$y)/20),
      ]
      
      plot(flowData, type = "l", main = "Local Peaks")
      points(possiblePeaks2$x, possiblePeaks2$y, col = "red", pch = 19)
      
      possiblePeaks3 <- findClusters(possiblePeaks2, 40, xVarMax)
      
      plot(flowData, type = "l", main = "Local Peaks")
      points(possiblePeaks3$x, possiblePeaks3$y, col = "red", pch = 19)
      
      possiblePeaks4 <- findPairs(possiblePeaks3, possiblePeaks3, 1.75-range_increase, 2.2+range_increase)
      possiblePeaks4 <- possiblePeaks4 %>%
        tidyr::drop_na() 
      possiblePeaks5 <- possiblePeaks4[!duplicated(possiblePeaks4$possiblePairY), ]
    }
  }else if(nrow(possiblePeaks5) > numSubPop){
    possiblePeaks5<-possiblePeaks5[1:numSubPop,] 
  }
  
  rangeLength <- nchar(format(xVarMax, scientific=FALSE))
  multiplier <- 10^(rangeLength-3)
  possiblePeaks6 <- doubletCheck(
    possiblePeaks5,
    possiblePeaks,
    10*multiplier,
    15*multiplier
  )
  
  tempDs<-possiblePeaks6
  if(nrow(tempDs) > 1 & !is.na(tempDs$g1G2DoubletCount[1])){
    pop1<-tempDs[1,]
    for(i in 2:nrow(tempDs)){
      popInQuestion<-tempDs[i,]
      if(!is.na(pop1$g2G2Doublet)){
        if(pop1$g2G2Doublet == popInQuestion$possiblePairX ){
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
  
  possiblePeaks6<-possiblePeaks6 %>% dplyr::mutate(
    data = sampleName,
    "propCellsUsed" = NA
  )
  
  if(!purrr::is_empty(possiblePeaks6)){
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
  }else{
    .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 0
    stop(paste0("Could not find ", numSubPop," subpopulations in this sample"))
  }
  

  if(nrow(possiblePeaks6) == 1){
    singleData<-possiblePeaks6
    finishedData<- NULL
    investigateData<- NULL
  }else{
    finishedData<-possiblePeaks6
    singleData<-NULL
    investigateData<- NULL
  }
  
  outputData(
    flowDir,
    singleData,
    finishedData,
    investigateData,
    xVariable,
    doublet = FALSE,
    saveGraph = TRUE
  )
  print("Done! - Check 'analysis' folder for results")
}

