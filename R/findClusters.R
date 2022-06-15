#' findClusters
#'
#' findClusters identifies G1 and G2 sub population pairing.
#'
#' @param ds The data set to analyze
#' @param clusterDist The smallest distance between two peaks/clusters
#' @param maxXValue the maximum of the x-axis
#' @export
#'
#' @examples
#' findClusters(
#'  ds = possible_peaks2,
#'  clusterDist = 40,
#'  maxXValue = maxXValue
#'  )
#'
findClusters = function(ds, clusterDist, maxXValue){

  #Finding initial distance between all identified peaks
  tempDs <- ds
  tempDs$cluster <- 1
  tempDs$distToNext <- c(diff(tempDs$x), 0)

  #Creating an distance requirement between two peaks
  maxDist <- maxXValue/clusterDist
  clusterNum <- 1
  #Going through each peak and checking if their distance to the next peak is
  #greater than the distance requirement 'maxDist'. If not, the peaks will be
  #grouped together in the same cluster.
  for (i in 1:nrow(tempDs)) {
    if(tempDs$distToNext[i] > maxDist){
      tempDs$cluster[i] <- clusterNum
      clusterNum <- clusterNum + 1
    } else {
      tempDs$cluster[i] <- clusterNum
    }
  }
  #If clusters have more than one peak identified, the tallest peak will be
  #selected as the peak in that cluster
  tempDs <- tempDs %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice(which.max(y))

  return(tempDs)

}
