#' massFlowGraph
#'
#' A function that permits the user to plot flow frames using the output from
#' `flowPeakDetection()`. The user is prompted to choose which set of frames to
#' plot (all, flagged, or doublet flagged). The chosen flow frames are plotted
#' using the grid functionality in `flowLineGraph()`. The user must ensure the
#' flow directory and file path correspond to the same set of samples for plotting.
#'
#' @param xVariable The fluorescence channel on the x-axis
#' @param flowDir The directory of the gated .fcs data
#' @param filePath The directory of the output .csv file from
#' `flowPeakDetection()` that corresponds to the gated data
#' @param fileName To label the PDF file produced internally by `flowLineGraph()`
#' 
#' @import ggplot2
#' @import tidyverse
#' @import gridExtra
#' @import dplyr
#' 
#' @return description
#' @export
#' 
#' @examples
#' massFlowGraph(
#'  xVariable = "FITC-A",
#'  flowDir = paste0(system.file(package = "PloidyPeaks"), "/gated_data/"),
#'  filePath = paste0(system.file(package = "PloidyPeaks"), "/analysis/"),
#'  fileName = "gatedDataRSE"
#'  )

massFlowGraph <- function(
    xVariable = "FITC-A",
    flowDir = NA,
    filePath = NA,
    fileName = NA
){
  # check for valid file path
  if(!file.exists(filePath)){
    stop("Invalid filepath")
  }
  
  # read in the indicated .csv file 
  data <- read.csv(filePath, header = TRUE)

  # round the final RSE values to 1 decimal
  data$finalRSE <- round(data$finalRSE, 1)
  
  # prompt the user to select what samples they want to plot
  input <- as.numeric(menu(choices = c("Plot All Samples", "Plot Flagged Samples", "Plot Doublet Flagged Samples"), 
                           title = "What would you like to plot?"))
  
  # check if user specifies file name 
  if(!is.na(fileName)){
    output_fileName <- fileName
  }
  
  if(input == 1){
    # check if user specifies file name, give generic name if they don't
    if(is.na(fileName)){
      output_fileName <- "allSamples"
      
    }
  }else if(input == 2){
    # filter data to flagged samples
    data <- data[data$investigate == 1, ]
    
    # check if user specifies file name, give generic name if they don't
    if(is.na(fileName)){
      output_fileName <- "flaggedSamples"
    }
  }else if(input == 3){
    # filter data to doublet flagged samples
    data <- data[data$doublet == 1, ]
    
    # check if user specifies file name, give generic name if they don't
    if(is.na(fileName)){
      output_fileName <- "doubletSamples"
    }
  } else if(input == 0){
    message <- "Exiting the mass plot function"
  }
  
  # get final RSE values for selected samples
  annote <- append(data$finalRSE, "RSE", after = 0)
  
  # if G1_2 peaks are present
  if("G1_2" %in% colnames(data) & (input == 1 | input == 2 | input == 3)){
    flowLineGraph(
      flowDir = flowDir,
      flowSamples = data$Sample,
      xVariable = xVariable,
      grid = TRUE,
      vertLine1 = data$G1_1,
      vertLine2 = data$G1_2,
      annotations = annote,
      fileName = output_fileName
    )
    message <- "Plotting of samples have been completed."
  } else if(input == 1 | input == 2 | input == 3){
    flowLineGraph(
      flowDir = flowDir,
      flowSamples = data$Sample,
      xVariable = xVariable,
      grid = TRUE,
      vertLine1 = data$G1_1,
      annotations = annote,
      fileName = output_fileName
    )
    message <- "Plotting of samples have been completed."
  }
  return(message)
}
