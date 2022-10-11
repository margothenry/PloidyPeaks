#' doubletCheck
#'
#' For single populations, doubletCheck is used to identify doublets.
#' Doublets are an excess of cells around the G1+G2 range.
#'
#' @param doubletCheckDs The data set to analyze for doublets
#' @param peaks The data peaks that have been identified
#' @param g1G2Range A range around the G1+G2 location on the x-axis,
#' i.e identifying the G1+G2 doublet in the range g1G2Range - G1+G2, g1G2Range + G1+G2
#' @param g2G2Range A range around the G2+G2 location on the x-axis,
#' i.e identifying the G2+G2 doublet in the range g2G2Range - G2+G2, g2G2Range + G2+G2
#' @export
#'
#' @examples
#' \dontrun{
#' doubletCheck(
#'  doubletCheckDs = single_ds,
#'  peaks = peaks_ds,
#'  g1G2Range = 10,
#'  g2G2Range = 15
#' )
#' }

doubletCheck = function(doubletCheckDs, peaks, g1G2Range, g2G2Range){
  ##Removing NOTE 'no visible binding for global variable'
  x<-possiblePairX<-NULL
  ##creating lower bounds and upper bounds for peaks that could be classified as doublets
  doubletCheckDs2 <- doubletCheckDs %>% dplyr::mutate(
    g3LL=x + possiblePairX - g1G2Range,
    g3UL=x + possiblePairX + g1G2Range,
    g4LL=possiblePairX + possiblePairX - g2G2Range,
    g4UL=possiblePairX + possiblePairX + g2G2Range
  )

  ##finding the peaks that are in the G1+G2 range
  doubletCheckDs2$g1G2Doublet <- NA
  doubletCheckDs2$g1G2DoubletCount <- NA
  for( i in 1:nrow(doubletCheckDs2)){
    possible_ <- which(
      peaks$x >= doubletCheckDs2$g3LL[i] &
        peaks$x <= doubletCheckDs2$g3UL[i]
    )

    if(length(possible_) > 1){
      possible_ <- which(
        peaks$y == max(peaks[possible_,]$y) &
        peaks$x >= doubletCheckDs2$g3LL[i]   
        )
      possible_ <-possible_[1]
    }

    if(!purrr::is_empty(possible_)){
      if(possible_ != i) {
        ##If there are more than one peak identified, we pick the tallest of those peaks
        maxPossible_ <- peaks[possible_, ]
        maxPossibleRows <- maxPossible_[
          order(maxPossible_$y, decreasing = TRUE),
          ][1, ]

        doubletCheckDs2$g1G2Doublet[i] <- maxPossibleRows$x
        doubletCheckDs2$g1G2DoubletCount[i] <- maxPossibleRows$y
      }else{
        possible_ <- NA
      }
    }
  }

  ##similarly for G2+G2
  ##finding the peaks that are in the G2+G2 range
  doubletCheckDs2$g2G2Doublet <- NA
  doubletCheckDs2$g2G2DoubletCount <- NA
  for( i in 1:nrow(doubletCheckDs2)){
    possible_ <- which(
      peaks$x >= doubletCheckDs2$g4LL[i] &
        peaks$x <= doubletCheckDs2$g4UL[i]
    )
    if(length(possible_) > 1){
      possible_ <- which(
        peaks$y == max(peaks[possible_,]$y) &
        peaks$x >= doubletCheckDs2$g4LL[i]  
        )
      possible_ <-possible_[1]
    }

    
    if(!purrr::is_empty(possible_)){
      if(possible_ != i) {
        maxPossible_ <- peaks[possible_,]
        maxPossibleRows <- maxPossible_[
          order(maxPossible_$y, decreasing = TRUE),
          ][1, ]

        doubletCheckDs2$g2G2Doublet[i] <- maxPossibleRows$x
        doubletCheckDs2$g2G2DoubletCount[i] <- maxPossibleRows$y
      }else{
        possible_ <- NA
      }
    }
  }

  return(doubletCheckDs2)
}



