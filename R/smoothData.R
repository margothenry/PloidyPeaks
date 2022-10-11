#' smoothData
#'
#' A function that reads in the gated data and applies a smoothing factor
#' to the data.
#'
#' @param flowDs The directory where the data is located. If NA then a window
#' will prompt the user to select a folder.
#' @param smoothLevel The level of smoothing applied to the data
#' @param xVariable The fluorescence channel you are interested in analyzing
#' @export
#'
#' @examples
#' \dontrun{
#'smoothData(flowDs, "FITC-A", 5)
#'}

smoothData = function(flowDs, xVariable, smoothLevel){

  ##Get counts and breaks from histogram
  histData <- hist(flowDs@exprs[, xVariable], breaks=256, plot=FALSE)

  ##Data that will be smoothed
  data <- histData$counts
  ##Apply smoothing to the counts with 'smoothLevel'
  smoothedDs <- zoo::rollmean(data, k=smoothLevel, fill=0)

  ##Create data frame with smoothed data
  ds <- data.frame(
    x=histData$breaks,
    y=c(0, smoothedDs)
  )
  return(ds)

}

