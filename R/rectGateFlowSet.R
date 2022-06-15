#' rectGateFlowSet
#'
#' This function will gate your raw flow set and save your gated data, a plot
#' of your gated data, and a plot of the original data
#'
#' @param rawDir The directory of the raw .fcs data
#' @param xVariable The fluorescence channel on the x axis - by default "FL1-A"
#' @param yVariable The fluorescence channel on the y axis - by default "SSC-A"
#' @param xMinValue The lower bound x value for the gate - by default 10000
#' @param xMaxValue The upper bound x value for the gate - by default 900000
#' @param yMinValue The lower bound y value for the gate - by default 10000
#' @param yMaxValue The upper bound y value for the gate - by default 900000
#' @param savePlot A side by side graph comparison the raw data and the gated
#'  data - by default TRUE
#'
#' @return A .fcs of the gated data
#' @import ggcyto
#' @import patchwork
#' @import tcltk
#' @export
#'
#' @examples
#' rectGateFlowSet(
#'  rawDir = NA,
#'  xVariable = "FL1-A",
#'  yVariable = "SSC-A",
#'  xMinValue = 10000,
#'  xMaxValue = 900000,
#'  yMinValue = 10000,
#'  yMaxValue = 900000,
#'  savePlot = TRUE
#')
#

rectGateFlowSet = function(
  rawDir = NA,
  xVariable = "FL1-A",
  yVariable = "SSC-A",
  xMinValue = 10000,
  xMaxValue = 900000,
  yMinValue = 10000,
  yMaxValue = 900000,
  savePlot = TRUE
){

  if(is.na(rawDir)){
    getwd()
    rawDir <- tcltk::tclvalue(tkchooseDirectory())
  }

  flowSet <- read.flowSet(
    path = rawDir,
    transformation=FALSE,
    truncate_max_range = TRUE
    )
  #Progress bar iterations
  total <- length(flowSet)
  #Create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  setwd(rawDir)
  subDir <- "gated_data"
  dir.create(file.path(dirname(rawDir), subDir), showWarnings = FALSE)

  if(savePlot == TRUE){
    gatePlotDir <- "plotted_data"
    dir.create(file.path(dirname(rawDir), gatePlotDir), showWarnings = FALSE)
    plotOutFile <- file.path(dirname(rawDir),gatePlotDir)
  }

  gatedCellsOut <- data.frame(
    matrix(nrow = 0, ncol = 2)
  )
  colnames(gatedCellsOut) <- c("Data", "% of cells gated out")

  for(i in 1:length(flowSet)){

    flowData <- flowSet[[i]]

    autoGate <- paste0(
      'rectGate <- rectangleGate(
          filterId=\"Fluorescence Region\",\"',
      xVariable,'\" = c(',xMinValue,',', xMaxValue,'),\"',
      yVariable, '\" = c(',yMinValue,',', yMaxValue,')
        )'
    )

    eval(parse(text = autoGate))

    gatedFlowData <- Subset(flowData, rectGate)
    frameName <- gatedFlowData@description[["GUID"]]
    numCellsGatedOut <- round(
      100 - (
        length(gatedFlowData@exprs[,xVariable])/
        length(flowData@exprs[,xVariable])
        )*100,1
      )

    gatedCellsFlow <- data.frame(
      matrix(nrow = 0, ncol = 2)
    )
    colnames(gatedCellsFlow) <- c("Data", "% of cells gated out")
    gatedCellsFlow[1, ] <- c(frameName, numCellsGatedOut)
    gatedCellsOut <- rbind(
      gatedCellsOut,
      gatedCellsFlow
    )

    outFile <- file.path(dirname(rawDir),subDir, frameName)
    write.FCS(gatedFlowData, outFile)

    if(savePlot == TRUE){

      flowData@description[["GUID"]] <- "Raw data"
      rawDataPlot <- autoplot(
        flowData, xVariable, yVariable, bins = 64
        ) + geom_gate(rectGate) + ggcyto_par_set(limits = "data")
      gatedFlowData@description[["GUID"]] <- "Gated data"
      gatedDataPlot <- autoplot(
        gatedFlowData, xVariable, yVariable,  bins = 64
        ) + ggcyto_par_set(limits = "instrument")
      combinedPlot <- as.ggplot(rawDataPlot) + as.ggplot(gatedDataPlot)
      gatedFlowData@description[["GUID"]] <- frameName
      png(paste0(plotOutFile,"/",frameName,'.png'), width = 600, height = 400)
      print(combinedPlot)
      dev.off()

    }
    setTxtProgressBar(pb, i)
  }
  write.csv(
    gatedCellsOut, paste0(dirname(getwd()),"/ProportionOfCellsGatedOut.csv")
    )
  close(pb)
}
