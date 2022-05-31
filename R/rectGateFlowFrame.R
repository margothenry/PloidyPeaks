#' rectGateFlowFrame
#'
#' This function will gate your raw .fcs file and save your gated data, a plot
#' of your gated data, and a plot of the original data
#'
#' @param raw_dir The directory of the raw .fcs data
#' @param flow_name The name of the .fcs file
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
#' rectGateFlowFrame(
#'  raw_dir = NA,
#'  flow_name = "A10-t0.fcs",
#'  x_variable = "FL1-A",
#'  y_variable = "SSC-A",
#'  x_min_value = 10000,
#'  x_max_value = 900000,
#'  y_min_value = 10000,
#'  y_max_value = 900000,
#'  save_plot = TRUE
#')
#
rectGateFlowFrame = function(
  raw_dir = NA,
  flow_name,
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

  flow_data <- read.FCS( paste0(raw_dir,"/",flow_name), transformation=FALSE)

  auto_gate <- paste0(
    'rectGate <- rectangleGate(
          filterId=\"Fluorescence Region\",\"',
    x_variable,'\" = c(',x_min_value,',', x_max_value,'),\"',
    y_variable, '\" = c(',y_min_value,',', y_max_value,')
        )'
  )

  eval(parse(text = auto_gate))

  gated_flow_data <- Subset(flow_data, rectGate)
  num_cells_gated_out <-  round(
    100 - (
      length(gated_flow_data@exprs[,x_variable]) /
      length(flow_data@exprs[,x_variable]))*100,1
    )
  print( paste0(num_cells_gated_out,"% of the cells were gated out") )
  setwd(raw_dir)
  subDir <- "gated_data"
  dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)

  outFile <- file.path(dirname(getwd()),subDir, flow_name)
  write.FCS(gated_flow_data, outFile)

  if(save_plot == TRUE){

    gate_plotDir <- "plotted_data"
    dir.create(file.path(dirname(getwd()), gate_plotDir), showWarnings = FALSE)
    plotOutFile <- file.path(dirname(getwd()),gate_plotDir)

    flow_data@description[["GUID"]] <- "Raw data"
    raw_data_plot <- autoplot(
      flow_data, x_variable, y_variable, bins = 64
      ) + geom_gate(rectGate) + ggcyto_par_set(limits = "data")
    gated_flow_data@description[["GUID"]] <- "Gated data"
    gated_data_plot <- autoplot(
      gated_flow_data, x_variable, y_variable,  bins = 64
      ) + ggcyto_par_set(limits = "instrument")
    combined_plot <- as.ggplot(raw_data_plot) + as.ggplot(gated_data_plot)

    gated_flow_data@description[["GUID"]] <- flow_name
    png(paste0(plotOutFile,"/",flow_name,'.png'), width = 600, height = 400)
    print(combined_plot)
    dev.off()

  }

}
