#' RSEOutlierDetection
#'
#' A function that permits the user to analyze the RSE values provided by the
#' output from `flowPeakDetection()`. A new column is added to the .csv file
#' that indicates whether each sample is an outlier based on its RSE. The alpha
#' parameter allows the user to specify the quantile percentage they'd like to
#' use as a cutoff. Each sample that is an outlier is plotted in a grid
#' formation using `flowLineGraph()` where the dots on each plot are the G1
#' (blue) and G2 (green) peaks with the respective RSE in the top right corner.
#'
#' @param xVariable The fluorescence channel on the x axis
#' @param flowDir The directory of the gated .fcs data
#' @param filePath The csv analysis file that pairs with the gated .fcs data
#' @param fileName To label the PDF file of plotted outlier samples
#' @param alpha The top percentage used to identify outliers (i.e 0.05)
#' 
#' @import ggplot2
#' @import tidyverse
#' @import gridExtra
#' @import cowplot
#' @import moments
#' @import dplyr
#' 
#' @return a .csv with an added column for outlier RSE about each sample, a
#' histogram to show the distribution of RSE values, and a PDF file of plotted
#' samples that were identified as outliers. Note: only the histogram is saved
#' if the skewness is less than zero.
#' @export
#'
#' @examples
#' RSEOutlierDetection(
#'  xVariable = "FITC-A",
#'  flowDir = paste0(system.file(package = "PloidyPeaks"), "/gated_data/"),
#'  filePath = paste0(system.file(package = "PloidyPeaks"), "/analysis/"),
#'  fileName = "gatedDataRSE",
#'  alpha = 0.05
#' )

RSEOutlierDetection = function(
    xVariable = "FL1-A",
    flowDir = NA,
    filePath = NA,
    fileName = "RSEOutlier",
    alpha = 0.05
){
  # check if the file path exists
  if(!file.exists(filePath)){
    stop("Invalid filepath")
  }
  # ensure alpha is between 0 and 1
  if(!is.numeric(alpha)){
    stop("Your alpha is not numeric")
  }else if(!between(alpha, 0, 1)){
    stop("Your alpha is not a decimal between 0 and 1")
  }
  
  # read in data
  data <- read.csv(filePath)
  
  # compute skew
  skewed <- skewness(data$finalRSE)
  
  # plot RSE values as a histogram w/ skew value
  RSEHist <- ggplot(data, aes(x = finalRSE)) +
    stat_function(fun = dnorm, args = list(mean = mean(data$finalRSE),
                                           sd = sd(data$finalRSE))) +
    geom_vline(aes(xintercept = mean(data$finalRSE), linetype = 'Mean'),
               col = "red") +
    geom_vline(aes(xintercept = median(data$finalRSE), linetype = 'Median'),
               col = "blue") +
    scale_linetype_manual(name = 'Lines', values = c('Mean' = 1, 'Median' = 1),
                          guide = guide_legend(override.aes = list(colour =
                                                          c("red", "blue")))) +
    labs(title = "RSE Values") +
    theme_bw() + annotate("text", x = (max(data$finalRSE) * 0.85), y = Inf,
                          label = paste("Skewness = ", round(skewed, 4)),
                          colour = "red", size = 3, vjust = 2)
  
  # if skew, add "outlier RSE" column to flag top alpha % quantile
  if(skewed > 0){
    quantileValue <- quantile(data$finalRSE, probs = (1 - alpha))
    topPercent <- as.numeric(data$finalRSE > quantileValue)
    # add to data frame
    data$outlierRSE <- topPercent
    
    # plot grid of samples
    outliers <- data[data$outlierRSE == 1, ]
    outliers$finalRSE <- round(outliers$finalRSE, 2)
    if("G1_3" %in% colnames(outliers)){   # has 3 peaks
      peakVals <- outliers[, 2:7]
    }else{     # only has 2 peaks
      peakVals <- outliers[, 2:5]
    }
    row.names(peakVals) <- NULL
    annote <- append(outliers$finalRSE, "RSE", after = 0)
    flowLineGraph(
      flowDir = flowDir,
      flowSamples = outliers$Sample,
      xVariable = "FITC-A",
      grid = TRUE,
      annotations = annote,
      samplePeaks = peakVals,
      fileName = fileName
    )
    
    # save csv file
    write.csv(data, filePath, row.names = FALSE)
  }else{
    print("The skewness is less than zero. Only the histogram will be saved.")
  }
  
  fileName <- sub(" ", "_", fileName)  # clean up (no white space)
  ggsave(plot = RSEHist, filename = paste0(fileName, "_plot.jpeg"))
  return("Outlier Analysis is Complete")
}
