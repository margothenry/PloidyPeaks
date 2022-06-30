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

  ##Removing NOTE 'no visible binding for global variable'
  rectGate<-NULL
  
  if(is.na(rawDir)){
    getwd()
    rawDir <- tclvalue(tkchooseDirectory())
  }
  
  if(purrr::is_empty(rawDir)){
    stop("Your directory is empty")
  }
  
  xpectr::suppress_mw(
    flowSet <- flowCore::read.flowSet(
      path=rawDir,
      transformation=FALSE,
      truncate_max_range=TRUE
    )
  )
  
  #Progress bar iterations
  total <- length(flowSet)
  #Create progress bar
  pb <- txtProgressBar(min=0, max=total, style=3)

  setwd(rawDir)
  subDir <- "gated_data"
  dir.create(file.path(dirname(rawDir), subDir), showWarnings=FALSE)

  if(savePlot == TRUE){
    gatePlotDir <- "plotted_data"
    dir.create(file.path(dirname(rawDir), gatePlotDir), showWarnings=FALSE)
    plotOutFile <- file.path(dirname(rawDir), gatePlotDir)
  }

  gatedCellsOut <- data.frame(
    matrix(nrow=0, ncol=2)
  )
  colnames(gatedCellsOut) <- c("Data", "% of cells gated out")

  for(i in 1:length(flowSet)){
    
    # if(
    #   tools::file_ext(flowSet[[i]]) %in% c(
    #     "csv", "xls", "xlsx","html", "ppt", "pptx" ,"pdf", "doc", "docx"
    #   )
    # ){
    #   stop(paste0("The flow frame ",flowSet[i]," is does not seem to be in the valid flow format. Consider changing the format or removing from the folder."))
    # }
    
    flowData <- flowSet[[i]]
    
    ##Checking to see if the user input the correct X and Y variable
    if(!xVariable %in% flowData@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    if(!yVariable %in% flowData@parameters@data$name){
      stop("Your Y variable is not in the dataset")
    }
    
    ##Checking the gating parameters are in the dataset
    if(xMaxValue > max(flowData@exprs[,xVariable])){
      stop("Your xMaxValue exceeds the range of the flow frame, consider a new value")
    }
    
    if(xMinValue < min(flowData@exprs[,xVariable])){
      stop("Your xMinValue exceeds the range of the flow frame, consider a new value")
    }
    
    if(yMaxValue > max(flowData@exprs[,yVariable])){
      stop("Your yMaxValue exceeds the range of the flow frame, consider a new value")
    }
    
    if(yMinValue < min(flowData@exprs[,yVariable])){
      stop("Your yMinValue exceeds the range of the flow frame, consider a new value")
    }
    
    ##Creating the gate
    autoGate <- paste0(
      'rectGate <- flowCore::rectangleGate(
          filterId=\"Fluorescence Region\",\"',
      xVariable,'\" = c(',xMinValue,',', xMaxValue,'),\"',
      yVariable, '\" = c(',yMinValue,',', yMaxValue,')
        )'
    )

    eval(parse(text=autoGate))
    ##Subsetting the data that is in the gate
    gatedFlowData <- flowCore::Subset(flowData, rectGate)
    ##Finding the % of cells gated out
    frameName <- gatedFlowData@description[["GUID"]]
    numCellsGatedOut <- round(
      100 - (
        length(gatedFlowData@exprs[, xVariable])/
        length(flowData@exprs[, xVariable])
        )*100,1
      )

    gatedCellsFlow <- data.frame(
      matrix(nrow=0, ncol=2)
    )
    colnames(gatedCellsFlow) <- c("Data", "% of cells gated out")
    gatedCellsFlow[1, ] <- c(frameName, numCellsGatedOut)
    gatedCellsOut <- rbind(
      gatedCellsOut,
      gatedCellsFlow
    )
    ##Saving the gated data intob a folder
    outFile <- file.path(dirname(rawDir), subDir, frameName)
    flowCore::write.FCS(gatedFlowData, outFile)
    ##If TRUE plots will be created for the user and saved in a folder
    if(savePlot == TRUE){

      flowData@description[["GUID"]] <- "Raw data"
      rawDataPlot <- ggcyto::autoplot(
        flowData, xVariable, yVariable, bins=64
        ) + ggcyto::geom_gate(rectGate) + ggcyto::ggcyto_par_set(limits="data")
      gatedFlowData@description[["GUID"]] <- "Gated data"
      gatedDataPlot <- xpectr::suppress_mw(
        ggcyto::autoplot(
          gatedFlowData, xVariable, yVariable,  bins=64
        ) + ggcyto::ggcyto_par_set(limits="instrument")
      )
      combinedPlot <- ggcyto::as.ggplot(rawDataPlot) + ggcyto::as.ggplot(gatedDataPlot)
      gatedFlowData@description[["GUID"]] <- frameName
      png(paste0(plotOutFile, "/", frameName, '.png'), width=600, height=400)
      print(combinedPlot)
      dev.off()

    }
    setTxtProgressBar(pb, i)
  }
  write.csv(
    gatedCellsOut, paste0(dirname(getwd()), "/ProportionOfCellsGatedOut.csv")
    )
  close(pb)
}
