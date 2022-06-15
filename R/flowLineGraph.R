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
#'   flowControl = "gated_FH-B11",
#'   flowSamples = c("gated_T1-D08","gated_FH-B08", "gated_T1-D10", "gated_A13-E10"),
#'   xVariable = "FITC-A",
#'   flowDir = "FlowData/gated_data"
#'  )


flowLineGraph = function(flowDir = NA, flowControl, flowSamples, xVariable){

  if(is.na(flowDir)){
    getwd()
    flowDir <- tclvalue(tkchooseDirectory())
  }

  #samples dataset
  sampleDs <- c()
  for(k in 1:length(flowSamples)){
    flowName <- read.FCS(
      paste0(flowDir,"/",flowSamples[k]), transformation=FALSE
      )
    ds <- smoothData( flowName, xVariable, 5)
    ds$Data <- flowName@description[["GUID"]]
    sampleDs <- rbind(
      sampleDs,
      ds
    )
  }

  #control dataset
  flowNameControl <- read.FCS(
    paste0(flowDir,"/",flowControl), transformation=FALSE
    )
  controlDs <- smoothData(flowNameControl, xVariable, 5)
  controlDs$Data <- flowNameControl@description[["GUID"]]

  #plotting
  flowPlot <- ggplot() +
    geom_line(data=sampleDs, aes(x=x, y=y, group=Data, color = Data))+
    geom_line(data=controlDs, aes(x=x, y=y, group=2), size = 1, color='black')+
    ylab("Counts")+
    xlab(xVariable)+
    theme_bw()

  return(flowPlot)
}




