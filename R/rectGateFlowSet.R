#' rectGateFlowSet
#'
#' This function will gate your raw flow set and save your gated data, a plot
#' of your gated data, and a plot of the original data
#'
#' @param raw_dir The directory of the raw .fcs data
#' @param x_variable The fluorescence channel on the x axis - by default "FL1-A"
#' @param y_variable The fluorescence channel on the y axis - by default "SSC-A"
#' @param x_min_value The lower bound x value for the gate - by default 10000
#' @param x_max_value The upper bound x value for the gate - by default 900000
#' @param y_min_value The lower bound y value for the gate - by default 10000
#' @param y_max_value The upper bound y value for the gate - by default 900000
#' @param save_plot A side by side graph comparison the raw data and the gated
#'  data - by default TRUE
#'
#' @return A .fcs of the gated data
#' @export
#'
#' @examples
#' rectGateFlowSet(
#'  raw_dir = NA,
#'  x_variable = "FL1-A",
#'  y_variable = "SSC-A",
#'  x_min_value = 10000,
#'  x_max_value = 900000,
#'  y_min_value = 10000,
#'  y_max_value = 900000,
#'  save_plot = TRUE
#')
#

rectGateFlowSet = function(
  raw_dir = NA,
  x_variable = "FL1-A",
  y_variable = "SSC-A",
  x_min_value = 10000,
  x_max_value = 900000,
  y_min_value = 10000,
  y_max_value = 900000,
  save_plot = TRUE
){

  if(is.na(raw_dir)){
    getwd()
    raw_dir <- tclvalue(tkchooseDirectory())
  }

  flow_set <- read.flowSet(
    path = raw_dir,
    transformation=FALSE,
    truncate_max_range = TRUE
    )
  #Progress bar iterations
  total <- length(flow_set)
  #Create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  setwd(raw_dir)
  subDir <- "gated_data"
  dir.create(file.path(dirname(raw_dir), subDir), showWarnings = FALSE)

  if(save_plot == TRUE){
    gate_plotDir <- "plotted_data"
    dir.create(file.path(dirname(raw_dir), gate_plotDir), showWarnings = FALSE)
    plotOutFile <- file.path(dirname(raw_dir),gate_plotDir)
  }

  GatedCellsOut <- data.frame(
    matrix(nrow = 0, ncol = 2)
  )
  colnames(GatedCellsOut) <- c("Data", "% of cells gated out")

  for(i in 1:length(flow_set)){

    flow_data <- flow_set[[i]]

    auto_gate <- paste0(
      'rectGate <- rectangleGate(
          filterId=\"Fluorescence Region\",\"',
      x_variable,'\" = c(',x_min_value,',', x_max_value,'),\"',
      y_variable, '\" = c(',y_min_value,',', y_max_value,')
        )'
    )

    eval(parse(text = auto_gate))

    gated_flow_data <- Subset(flow_data, rectGate)
    frame_name <- gated_flow_data@description[["GUID"]]
    num_cells_gated_out <- round(
      100 - (
        length(gated_flow_data@exprs[,x_variable])/
        length(flow_data@exprs[,x_variable])
        )*100,1
      )

    GatedCellsFlow <- data.frame(
      matrix(nrow = 0, ncol = 2)
    )
    colnames(GatedCellsFlow) <- c("Data", "% of cells gated out")
    GatedCellsFlow[1, ] <- c(frame_name, num_cells_gated_out)
    GatedCellsOut <- rbind(
      GatedCellsOut,
      GatedCellsFlow
    )

    outFile <- file.path(dirname(raw_dir),subDir, frame_name)
    write.FCS(gated_flow_data, outFile)

    if(save_plot == TRUE){

      flow_data@description[["GUID"]] <- "Raw data"
      raw_data_plot <- autoplot(
        flow_data, x_variable, y_variable, bins = 64
        ) + geom_gate(rectGate) + ggcyto_par_set(limits = "data")
      gated_flow_data@description[["GUID"]] <- "Gated data"
      gated_data_plot <- autoplot(
        gated_flow_data, x_variable, y_variable,  bins = 64
        ) + ggcyto_par_set(limits = "instrument")
      combined_plot <- as.ggplot(raw_data_plot) + as.ggplot(gated_data_plot)
      gated_flow_data@description[["GUID"]] <- frame_name
      png(paste0(plotOutFile,"/",frame_name,'.png'), width = 600, height = 400)
      print(combined_plot)
      dev.off()

    }
    setTxtProgressBar(pb, i)
  }
  write.csv(
    GatedCellsOut, paste0(dirname(getwd()),"/ProportionOfCellsGatedOut.csv")
    )
  close(pb)
}
