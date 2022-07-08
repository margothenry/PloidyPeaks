#' flowLineGraph
#'
#' A function that permits the users to plot as many flow frames as they would
#' like on the same graph panel. Each of the colourful lines represent a
#' different flow frame, and the black line represents the control sample used
#' for comparison (control sample is indicated by the user).
#'
#' @param flowDir The directory where the data is located. If NA then a window
#' will prompt the user to select a folder.
#' @param flowControl The control sample
#' @param flowSamples The vector of samples
#' @param xVariable The fluorescence channel you are interested in visualizing
#' @import tcltk
#' @import ggplot2
#' @export
#'
#' @examples
#'  flowLineGraph(
#'   flowControl = "FH-B11",
#'   flowSamples = c("T1-D08","FH-B08", "T1-D10", "A13-E10"),
#'   xVariable = "FITC-A",
#'   flowDir = "FlowData/gated_data"
#'  )


flowLineGraph = function(
  flowDir = NA,
  flowControl = NA,
  flowSamples = NA,
  xVariable
){
  ##Removing NOTE 'no visible binding for global variable'
  x<-y<-Data<-NULL
  if(is.na(flowDir)){
    getwd()
    flowDir <- tclvalue(tkchooseDirectory())
  }
  
  ##Checking if the folder they selected is empty
  if(purrr::is_empty(flowDir)){
    stop("Your directory is empty")
  }
  
  ##samples dataset
  if(TRUE %in% !is.na(flowSamples)){
    if(length(flowSamples) > 1){
      
      sampleDs <- c()
      for(k in 1:length(flowSamples)){
        
        if(!flowSamples[k] %in% list.files(flowDir)){
          stop(paste0("The flow frame ",flowSamples[k]," is not in the folder, check on the spelling of flowName and/or make sure you selected the proper folder"))
        }
        xpectr::suppress_mw(
          flowName <- flowCore::read.FCS(
            paste0(flowDir, "/", flowSamples[k]), transformation=FALSE
          )
        )
        
        ##Checking to see if the user input the correct X variable
        if(!xVariable %in% flowName@parameters@data$name){
          stop("Your X variable is not in the dataset")
        }
        
        ds <- smoothData( flowName, xVariable, 5)
        ds$Data <- flowSamples[k]
        sampleDs <- rbind(
          sampleDs,
          ds
        )
      }
      
    }else{
      xpectr::suppress_mw(
        flowName <- flowCore::read.FCS(
          paste0(flowDir, "/", flowSamples), transformation=FALSE
        ) 
      )
      
      ##Checking to see if the user input the correct X variable
      if(!xVariable %in% flowName@parameters@data$name){
        stop("Your X variable is not in the dataset")
      }
      
      sampleDs <- smoothData(flowName, xVariable, 5)
      sampleDs$Data <- flowSamples
    }
    
    if(!is.na(flowControl)){
      if(!flowControl %in% list.files(flowDir)){
        stop(paste0("The flow frame ",flowControl," is not in the folder, check on the spelling of flowName and/or make sure you selected the proper folder"))
      }
      ##control dataset
      xpectr::suppress_mw(
        flowNameControl <- flowCore::read.FCS(
          paste0(flowDir, "/", flowControl), transformation=FALSE
        ) 
      )
      
      ##Checking to see if the user input the correct X variable
      if(!xVariable %in% flowNameControl@parameters@data$name){
        stop("Your X variable is not in the dataset")
      }
      
      controlDs <- smoothData(flowNameControl, xVariable, 5)
      controlDs$Data <- flowControl
      
      ##plotting
      flowPlot <- ggplot() +
        geom_line(data=sampleDs, aes(x=x, y=y, group=Data, color = Data))+
        geom_line(data=controlDs, aes(x=x, y=y, group=2), size = 1, color='black')+
        ylab("Counts")+
        xlab(xVariable)+
        theme_bw()
      
    }else{
      if(length(flowSamples) > 1){
        ##plotting
        flowPlot <- ggplot() +
          geom_line(data=sampleDs, aes(x=x, y=y, group=Data, color = Data))+
          ylab("Counts")+
          xlab(xVariable)+
          theme_bw()
      }else if(length(flowSamples) == 1){
        ##plotting
        flowPlot <- ggplot() +
          geom_line(data=sampleDs, aes(x=x, y=y), size = 1, color='black')+
          ylab("Counts")+
          xlab(xVariable)+
          theme_bw()
      }
      
    }
  }else{
    if(!flowControl %in% list.files(flowDir)){
      stop(paste0("The flow frame ",flowControl," is not in the folder, check on the spelling of flowName and/or make sure you selected the proper folder"))
    }
    ##control dataset
    xpectr::suppress_mw(
      flowNameControl <- flowCore::read.FCS(
        paste0(flowDir, "/", flowControl), transformation=FALSE
      ) 
    )
    
    ##Checking to see if the user input the correct X variable
    if(!xVariable %in% flowNameControl@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    controlDs <- smoothData(flowNameControl, xVariable, 5)
    controlDs$Data <- flowControl
    
    ##plotting
    flowPlot <- ggplot() +
      geom_line(data=controlDs, aes(x=x, y=y, group=2), size = 1, color='black')+
      ylab("Counts")+
      xlab(xVariable)+
      theme_bw()
  }
  
  
  return(flowPlot)
}


