#' findPairs
#'
#' findPairs identifies G1 and G2 sub population pairing.
#'
#' @param ds The data set with the possible peaks
#' @param LL The lower limit multiplier for the range of the G2 peak, i.e, LL*G1
#' @param UL The upper limit multiplier for the range of the G2 peak, i.e, UL*G1
#' @export
#'
#' @examples
#'findPairs(
#' ds = peaks_ds,
#' LL = 1.8,
#' UL = 2.1
#')

findPairs = function(ds, LL, UL){
  #creating upper and lowed bounds that are read into the function
  findingPairsDs <- ds %>% dplyr::mutate(
    LL = x*LL,
    UL = x*UL
  )
  #Finding the peaks that are in the in (LL, UL)
  #These are the possible peaks that will be considered for their G1/G2 pairing.
  #If there is more than one peak identified, we pick the tallest of those peaks
  findingPairsDs$possiblePairX <- NA
  findingPairsDs$possiblePairY <- NA
  for( i in 1:nrow(findingPairsDs)){
    possible_ <- which(
      findingPairsDs$x >= findingPairsDs$LL[i] &
        findingPairsDs$x <= findingPairsDs$UL[i]
    )

    if(length(possible_) > 1){
      possible_ <- possible_[1]
    }

    if(!is_empty(possible_)){
      if(possible_ != i) {
        maxPossible_ <- findingPairsDs[possible_,]
        maxPossibleRows <- maxPossible_[
          order(maxPossible_$y, decreasing = TRUE),
          ][1,]

        findingPairsDs$possiblePairX[i] <- maxPossibleRows$x
        findingPairsDs$possiblePairY[i] <- maxPossibleRows$y
      }else{
        possible_ <- NA
      }
    }
  }

  return(findingPairsDs)
}

