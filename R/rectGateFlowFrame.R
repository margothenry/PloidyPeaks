#' rectGateFlowFrame
#'
#' This function will gate your raw .fcs file and save your gated data, a plot
#' of your gated data, and a plot of the original data
#'
#' @param rawDir The directory of the raw .fcs data
#' @param flowName The name of the .fcs file
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
#' 
#' @import patchwork
#' @import tcltk
#' 
#' @export
#'
#' @examples
#' rectGateFlowFrame(
#'  rawDir = NA,
#'  flowName = "A10-t0.fcs",
#'  xVariable = "FL1-A",
#'  yVariable = "SSC-A",
#'  xMinValue = 10000,
#'  xMaxValue = 900000,
#'  yMinValue = 10000,
#'  yMaxValue = 900000,
#'  savePlot = TRUE
#')
#
rectGateFlowFrame = function(
  rawDir = NA,
  flowName,
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
    rawDir <- tclvalue(tkchooseDirectory())
  }

  flowData <- flowCore::read.FCS(
    paste0(rawDir,"/",flowName),
    transformation=FALSE,
    truncate_max_range = TRUE
    )

  autoGate <- paste0(
    'rectGate <- flowCore::rectangleGate(
          filterId=\"Fluorescence Region\",\"',
    xVariable,'\" = c(',xMinValue,',', xMaxValue,'),\"',
    yVariable, '\" = c(',yMinValue,',', yMaxValue,')
        )'
  )

  eval(parse(text = autoGate))

  gatedFlowData <- flowCore::Subset(flowData, rectGate)
  numCellsGatedOut <-  round(
    100 - (
      length(gatedFlowData@exprs[,xVariable]) /
      length(flowData@exprs[,xVariable]))*100,1
    )
  print( paste0(numCellsGatedOut,"% of the cells were gated out") )
  setwd(rawDir)
  subDir <- "gated_data"
  dir.create(file.path(dirname(rawDir), subDir), showWarnings = FALSE)

  outFile <- file.path(dirname(rawDir),subDir, flowName)
  flowCore::write.FCS(gatedFlowData, outFile)

  if(savePlot == TRUE){

    gatePlotDir <- "plotted_data"
    dir.create(file.path(dirname(rawDir), gatePlotDir), showWarnings = FALSE)
    plotOutFile <- file.path(dirname(rawDir),gatePlotDir)

    flowData@description[["GUID"]] <- "Raw data"
    rawDataPlot <- ggcyto::autoplot(
      flowData, xVariable, yVariable, bins = 64
      ) + ggcyto::geom_gate(rectGate) + ggcyto::ggcyto_par_set(limits = "data")
    gatedFlowData@description[["GUID"]] <- "Gated data"
    gatedDataPlot <- ggcyto::autoplot(
      gatedFlowData, xVariable, yVariable,  bins = 64
      ) + ggcyto::ggcyto_par_set(limits = "instrument")
    combinedPlot <- ggcyto::as.ggplot(rawDataPlot) + ggcyto::as.ggplot(gatedDataPlot)

    gatedFlowData@description[["GUID"]] <- flowName
    png(paste0(plotOutFile,"/",flowName,'.png'), width = 600, height = 400)
    print(combinedPlot)
    dev.off()

  }

}
