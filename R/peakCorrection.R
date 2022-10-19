#' peakCorrection
#'
#' peakCorrection is a function created for the user to rerun the analysis on 
#' a sample with a certain number of subpopulations.
#'
#' @param xVariable The fluorescence channel on the x axis
#' @param flowDir The directory of the gated .fcs data
#' @param sampleName the name of the sample the user wants to analyze
#' @param numSubPop the number of subpopulations the user wants
#'  the algorithm to identify
#' @return a .csv with information abour each sample and nls graphs
#' @export
#'
#' @examples
#' \dontrun{
#' peakCorrection(
#'  xVariable = "FITC-A",
#'  flowDir = NA,
#'  sampleName = "A01-A01",
#'  numSubPop = 2
#'  )
#'  }
#'
peakCorrection = function(
    xVariable,
    flowDir = NA,
    sampleName,
    numSubPop = 2
){
    
    points<-NULL
    
    .GlobalEnv$logDs <- data.frame(
        matrix(nrow=0, ncol=3)
    )
    colnames(.GlobalEnv$logDs) <- c("Algorithm", "Data", "Success")
    
    algorithmNum <- 1
    
    .GlobalEnv$logDs[1, ] <- c(algorithmNum, sampleName, 0)
    
    ##Directory: If user do not hard codes their directory in the function
    ##a window will open and ask the user to pick their directory.
    ##This should be where the data is located
    if(is.na(flowDir)){
        getwd()
        flowDir <- tclvalue(tkchooseDirectory())
    }
    
    flowName<-flowCore::read.FCS(
        paste0(flowDir, "/", sampleName), transformation=FALSE
    )
    
    # Warning Check for X variable
    if(!xVariable %in% flowName@parameters@data$name){
        stop("Your X variable is not in the dataset")
    }
    
    # Warning Check for Sample Name
    if(!sampleName %in% list.files(flowDir)){
        stop(
            "Your Sample is not in the list of flow frames in the directory
            provided"
        )
    }
    
    # Smoothing Data
    flowData<-.smoothData( flowName, xVariable, 5)
    
    # Applying peak algorithm
    localPeaks<-detect_localmaxima(flowData$y, 3)
    possiblePeaks<-flowData[localPeaks, ]
    
    possiblePeaks2<-possiblePeaks[
        which(possiblePeaks$y > max(possiblePeaks$y)/20),
    ]
    
    xVarMax <- max(flowData$x)
    possiblePeaks3<-.findTruePeaks(possiblePeaks2, 40, xVarMax)
    
    possiblePeaks4<-.findPairs(possiblePeaks3, possiblePeaks3, 1.75, 2.2)
    possiblePeaks4<-possiblePeaks4 %>%
        tidyr::drop_na() 
    possiblePeaks5<-possiblePeaks4[!duplicated(possiblePeaks4$possiblePairY), ]
    
    # If sub populations < numSubPop, find missing population
    if(nrow(possiblePeaks5) < numSubPop){
        range_increase<-0
        range_increase_val = 0.2
        while(!nrow(possiblePeaks5) == numSubPop){
            range_increase = range_increase + range_increase_val
            
            flowData <- .smoothData( flowName, xVariable, 4)
            
            
            localPeaks <- detect_localmaxima(flowData$y, 3)
            possiblePeaks <- flowData[localPeaks, ]
            
            possiblePeaks2 <- possiblePeaks[
                which(possiblePeaks$y > max(possiblePeaks$y)/20),
            ]
            
            plot(flowData, type = "l", main = "Local Peaks")
            points(possiblePeaks2$x, possiblePeaks2$y, col = "red", pch = 19)
            
            possiblePeaks3 <- .findTruePeaks(possiblePeaks2, 40, xVarMax)
            
            plot(flowData, type = "l", main = "Local Peaks")
            points(possiblePeaks3$x, possiblePeaks3$y, col = "red", pch = 19)
            
            possiblePeaks4 <- .findPairs(
                possiblePeaks3, 
                possiblePeaks3, 
                1.75-range_increase, 
                2.2+range_increase
            )
            possiblePeaks4 <- possiblePeaks4 %>%
                tidyr::drop_na() 
            possiblePeaks5 <- possiblePeaks4[
                !duplicated(possiblePeaks4$possiblePairY),
            ]
        }
    }else if(nrow(possiblePeaks5) > numSubPop){
        possiblePeaks5<-possiblePeaks5[seq_len(numSubPop), ] 
    }
    
    rangeLength <- nchar(format(xVarMax, scientific=FALSE))
    multiplier <- 10^(rangeLength-3)
    possiblePeaks6 <- .doubletCheck(
        possiblePeaks5,
        possiblePeaks,
        10*multiplier,
        15*multiplier
    )
    
    tempDs<-possiblePeaks6
    if(nrow(tempDs) > 1 & !is.na(tempDs$g1G2DoubletCount[1])){
        pop1<-tempDs[1,]
        for(i in 2:nrow(tempDs)){
            popInQuestion<-tempDs[i,]
            if(!is.na(pop1$g2G2Doublet)){
                if(pop1$g2G2Doublet == popInQuestion$possiblePairX ){
                    tempDs <- tempDs[-i,]
                }else{
                    tempDs$g1G2Doublet<-NA
                    tempDs$g2G2Doublet<-NA
                    tempDs$g1G2DoubletCount<-NA
                    tempDs$g2G2DoubletCount<-NA
                } 
            }
        }
    }
    possiblePeaks6<-tempDs
    rm(tempDs) 
    
    possiblePeaks6<-possiblePeaks6 %>% dplyr::mutate(
        data = sampleName,
        "propCellsUsed" = NA
    )
    
    if(!purrr::is_empty(possiblePeaks6)){
        .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 1
    }else{
        .GlobalEnv$logDs[nrow(.GlobalEnv$logDs), ]$Success <- 0
        errorMsg<-paste0(
            "Could not find ", numSubPop," subpopulations in this sample")
        stop(errorMsg)
        rm(errorMsg)
    }
    
    
    if(nrow(possiblePeaks6) == 1){
        singleData<-possiblePeaks6
        finishedData<- NULL
        investigateData<- NULL
    }else{
        finishedData<-possiblePeaks6
        singleData<-NULL
        investigateData<- NULL
    }
    
    .outputData(
        flowDir,
        singleData,
        finishedData,
        investigateData,
        xVariable,
        doubletFlag = FALSE,
        saveGraph = TRUE
    )
    print("Done! - Check 'analysis' folder for results")
}


## Helper functions within the wrapper function

## smoothData
.smoothData = function(flowDs, xVariable, smoothLevel){
    
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

## findTruePeaks
.findTruePeaks = function(ds, clusterDist, maxXValue){
    ##Removing NOTE 'no visible binding for global variable'
    cluster<-y<-NULL
    ##Finding initial distance between all identified peaks
    tempDs <- ds
    tempDs$cluster <- 1
    tempDs$distToNext <- c(diff(tempDs$x), 0)
    
    ##Creating an distance requirement between two peaks
    maxDist <- maxXValue/clusterDist
    clusterNum <- 1
    ##Going through each peak and checking if their distance to the next peak is
    ##greater than the distance requirement 'maxDist'. If not, the peaks will be
    ##grouped together in the same cluster.
    for (i in seq_len(nrow(tempDs))) {
        if(tempDs$distToNext[i] > maxDist){
            tempDs$cluster[i] <- clusterNum
            clusterNum <- clusterNum + 1
        } else {
            tempDs$cluster[i] <- clusterNum
        }
    }
    ##If clusters have more than one peak identified, the tallest peak will be
    ##selected as the peak in that cluster
    tempDs <- tempDs %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice(which.max(y))
    
    return(tempDs)
    
}

## findPairs
.findPairs = function(ds, peaks, LL, UL){
    ##Removing NOTE 'no visible binding for global variable'
    x<-NULL
    ##creating upper and lowed bounds that are read into the function
    findingPairsDs <- ds %>% dplyr::mutate(
        LL = x*LL,
        UL = x*UL
    )
    ##Finding the peaks that are in the in (LL, UL)
    ##These are the possible peaks that will be considered for their G1/G2 pairing
    ##If there is more than one peak identified, we pick the tallest
    ##of those peaks
    findingPairsDs$possiblePairX <- NA
    findingPairsDs$possiblePairY <- NA
    for( i in seq_len(nrow(findingPairsDs))){
        possible_ <- which(
            peaks$x >= findingPairsDs$LL[i] &
            peaks$x <= findingPairsDs$UL[i]
        )
        
        if(length(possible_) > 1){
            possible_ <- which(
                peaks$y == max(peaks[possible_,]$y) &
                peaks$x >= findingPairsDs$LL[i]
            )
            possible_<-possible_[1]
        }
        
        if(!purrr::is_empty(possible_)){
            if(possible_ != i) {
                maxPossible_ <- peaks[possible_, ]
                maxPossibleRows <- maxPossible_[
                    order(maxPossible_$y, decreasing = TRUE), ][1,]
                
                findingPairsDs$possiblePairX[i] <- maxPossibleRows$x
                findingPairsDs$possiblePairY[i] <- maxPossibleRows$y
            }else{
                possible_ <- NA
            }
        }
    }
    
    return(findingPairsDs)
}

## doubletCheck
.doubletCheck = function(doubletCheckDs, peaks, g1G2Range, g2G2Range){
    ##Removing NOTE 'no visible binding for global variable'
    x<-possiblePairX<-NULL
    ##creating lower bounds and upper bounds for peaks that
    ##could be classified as doublets
    doubletCheckDs2 <- doubletCheckDs %>% dplyr::mutate(
        g3LL=x + possiblePairX - g1G2Range,
        g3UL=x + possiblePairX + g1G2Range,
        g4LL=possiblePairX + possiblePairX - g2G2Range,
        g4UL=possiblePairX + possiblePairX + g2G2Range
    )
    
    ##finding the peaks that are in the G1+G2 range
    doubletCheckDs2$g1G2Doublet <- NA
    doubletCheckDs2$g1G2DoubletCount <- NA
    for( i in seq_len(nrow(doubletCheckDs2))){
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
                ##If there are more than one peak identified, we pick the 
                ##tallest of those peaks
                maxPossible_ <- peaks[possible_, ]
                maxPossibleRows <- maxPossible_[
                    order(maxPossible_$y, decreasing = TRUE), ][1, ]
                
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
    for( i in seq_len(nrow(doubletCheckDs2))){
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
                    order(maxPossible_$y, decreasing = TRUE), ][1, ]
                
                doubletCheckDs2$g2G2Doublet[i] <- maxPossibleRows$x
                doubletCheckDs2$g2G2DoubletCount[i] <- maxPossibleRows$y
            }else{
                possible_ <- NA
            }
        }
    }
    
    return(doubletCheckDs2)
}


## outputData
.outputData = function(
    flowDir,
    singleDs,
    finishedDs,
    investigateDs,
    xVariable,
    doubletFlag,
    saveGraph
){
    ##Removing NOTE 'no visible binding for global variable'
    x<-y<-.<-possiblePairX<-possiblePairY<-G1<-G1Count<-G2<-G2Count<-id<-NULL
    g1G2Doublet<-g1G2DoubletCount<-g2G2Doublet<-g2G2DoubletCount<-NULL
    residual<-residualDoublet<-Success<-Algorithm<-Data<-`doublet G1+G2`<-NULL
    `doublet G1+G2 count`<-`doublet G2+G2`<- `doublet G2+G2 count`<-NULL
    residual3Pop<-residual2Pop<-residualMultiple<-NULL
    ##If the algorithm classified each sub population in the first algorithm
    if(
        purrr::is_empty(finishedDs) &
        purrr::is_empty(investigateDs) &
        !purrr::is_empty(singleDs)
    ){
    ##Formatting the single pop data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
        dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    finalData=finalPart1 %>% dplyr::select(
        c(
            "data",
            "G1",
            "G2",
            "g1G2Doublet",
            "g2G2Doublet",
            "G1Count",
            "G2Count",
            "g1G2DoubletCount",
            "g2G2DoubletCount",
            "investigate"
        )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
        dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
        dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
        dplyr::mutate(id=rowid(data)) %>%
        tidyr::pivot_wider(
            names_from=id,
            values_from=c(G1, G1Count, G2, G2Count)
        )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    
    finalData3=finalData3 %>%
        dplyr::rename(
            `doublet G1+G2`=g1G2Doublet,
            `doublet G1+G2 count`=g1G2DoubletCount,
            `doublet G2+G2`=g2G2Doublet,
            `doublet G2+G2 count`=g2G2DoubletCount
        )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
        flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
        doubletRSE<-.popConfidenceDoublet(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        doubletRSE<-finalData3 %>% dplyr::mutate(
            residualDoublet=NA
        )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
        twoPopRSE<-.popConfidence2Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
            residual2Pop=NA
        )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
        threePopRSE<-.popConfidence3Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        threePopRSE<-finalData3 %>% dplyr::mutate(
            residual3Pop=NA
        )
    }
    
    finalData4<-sqldf::sqldf(
        "select ds.*,
            ds2.residual,
            ds3.residualDoublet,
            ds4.residual2Pop,
            ds5.residual3Pop
            from finalData3 ds
            left join initialRSE ds2
            on ds.data = ds2.data
            left join doubletRSE ds3
            on ds.data = ds3.data
            left join twoPopRSE ds4
            on ds.data = ds4.data
            left join threePopRSE ds5
            on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
        residualMultiple=ifelse(
            is.na(residual3Pop), residual2Pop, residual3Pop)
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
        finalRSE=do.call(
            pmin, 
            c(
                subset(., select=c(
                    residual,
                    residualDoublet,
                    residualMultiple)
                ), na.rm=TRUE
            )
        )
        ) %>% dplyr::rename(
        singleRSE=residual,
        doubletRSE=residualDoublet,
        multipleRSE=residualMultiple
        ) %>% dplyr::select(
        -c(residual2Pop, residual3Pop)
        )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
        if(
            !is.na(finalData5[i,]$multipleRSE) &
            finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
        ){
            finalData5[i,]$investigate = 1
        }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
        dplyr::select(-Success) %>% 
        dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
        dplyr::mutate(id=rowid(Data)) %>%
        tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
        Algorithm=do.call(
            pmax, 
            c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
        dplyr::select(Data, Algorithm) %>%
        dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
        "Sample" = "data"
    )
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
        finalData6=finalData6 %>%
            dplyr::select(
                -c(
                    `doublet G1+G2`,
                    `doublet G1+G2 count`,
                    `doublet G2+G2`,
                    `doublet G2+G2 count`
                )
            )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings=FALSE)
    write.csv(
        finalData6,
        paste0(file.path(dirname(getwd()), subDir), "/ploidyPeaksOutput.csv")
    )
    
    }else if(
        purrr::is_empty(investigateDs) &
        !purrr::is_empty(finishedDs) &
        !purrr::is_empty(singleDs)
    ){
    ##If the algorithm classified each sub population before prior the 4th
    ##algorithm
    
    ##Formatting the diploid data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
        dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    ##Formatting the data from the other peak algorithms
    finalPart2=finishedDs %>% data.frame() %>%
        dplyr::select(
            -c(
                "g3LL",
                "g3UL",
                "g4LL",
                "g4UL",
                "propCellsUsed",
                "cluster",
                "distToNext",
                "LL",
                "UL"
            )
        ) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    ##Merging the two datasets
    finalData=rbind(
        finalPart1,
        finalPart2
    ) %>% dplyr::select(
        c(
            "data",
            "G1",
            "G2",
            "g1G2Doublet",
            "g2G2Doublet",
            "G1Count",
            "G2Count",
            "g1G2DoubletCount",
            "g2G2DoubletCount",
            "investigate"
        )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
        dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
        dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
        dplyr::mutate(id=rowid(data)) %>%
        tidyr::pivot_wider(
            names_from=id,
            values_from=c(G1, G1Count, G2, G2Count))
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
        finalData3=finalData3 %>%
            dplyr::rename(
                `doublet G1+G2`=g1G2Doublet,
                `doublet G1+G2 count`=g1G2DoubletCount,
                `doublet G2+G2`=g2G2Doublet,
                `doublet G2+G2 count`=g2G2DoubletCount
            )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
        flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
        doubletRSE<-.popConfidenceDoublet(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        doubletRSE<-finalData3 %>% dplyr::mutate(
            residualDoublet=NA
        )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
        twoPopRSE<-.popConfidence2Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
            residual2Pop=NA
        )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
        threePopRSE<-.popConfidence3Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        threePopRSE<-finalData3 %>% dplyr::mutate(
            residual3Pop=NA
        )
    }
    
    finalData4<-sqldf::sqldf(
        "select ds.*,
            ds2.residual,
            ds3.residualDoublet,
            ds4.residual2Pop,
            ds5.residual3Pop
            from finalData3 ds
            left join initialRSE ds2
            on ds.data = ds2.data
            left join doubletRSE ds3
            on ds.data = ds3.data
            left join twoPopRSE ds4
            on ds.data = ds4.data
            left join threePopRSE ds5
            on ds.data = ds5.data"
        )
    
    finalData4<-finalData4 %>% dplyr::mutate(
        residualMultiple=ifelse(
            is.na(residual3Pop),
            residual2Pop, residual3Pop
        )
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
        finalRSE=do.callSS(
            pmin, 
            c(
                subset(.,
                    select=c(
                        residual, 
                        residualDoublet, 
                        residualMultiple)
                ), na.rm=TRUE
            )
        )
        ) %>% dplyr::rename(
            singleRSE=residual,
            doubletRSE=residualDoublet,
            multipleRSE=residualMultiple
        ) %>% dplyr::select(
            -c(residual2Pop, residual3Pop)
        )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
        if(
            !is.na(finalData5[i,]$multipleRSE) &
            finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
        ){
            finalData5[i,]$investigate = 1
        }
        
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
        dplyr::select(-Success) %>% 
        dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
        dplyr::mutate(id=rowid(Data)) %>%
        tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
        Algorithm=do.call(
            pmax, 
            c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
        dplyr::select(Data, Algorithm) %>%
        dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
        "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
        finalData6=finalData6 %>%
            dplyr::select(
                -c(
                    `doublet G1+G2`,
                    `doublet G1+G2 count`,
                    `doublet G2+G2`,
                    `doublet G2+G2 count`
                )
            )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)
    write.csv(
        finalData6,
        paste0(file.path(dirname(getwd()), subDir), "/ploidyPeaksOutput.csv")
    )
    
    }else if(
        purrr::is_empty(singleDs) &
        !purrr::is_empty(finishedDs) &
        !purrr::is_empty(investigateDs)
    ){
    ##If the algorithm classified each sub population before prior the 4th
    ##algorithm
    
    ##Formatting the data from the other peak algorithms
    finalData1<-finishedDs %>% data.frame() %>%
        dplyr::select(
            -c(
                "g3LL",
                "g3UL",
                "g4LL",
                "g4UL",
                "propCellsUsed",
                "cluster",
                "distToNext",
                "LL",
                "UL"
            )
        ) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    finalPart2<-investigateDs %>% data.frame() %>%
        dplyr::select(
            -c(
                "g3LL",
                "g3UL",
                "g4LL",
                "g4UL",
                "cluster",
                "distToNext",
                "LL",
                "UL"
            )
        ) %>%
        dplyr::mutate(
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    ##Merging the two datasets
    finalData<-rbind(
        finalPart1,
        finalPart2
    ) %>% dplyr::select(
        c(
            "data",
            "G1",
            "G2",
            "g1G2Doublet",
            "g2G2Doublet",
            "G1Count",
            "G2Count",
            "g1G2DoubletCount",
            "g2G2DoubletCount",
            "investigate"
        )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
        dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
        dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
        dplyr::mutate(id=rowid(data)) %>%
        tidyr::pivot_wider(
            names_from=id,
            values_from=c(G1, G1Count, G2, G2Count)
        )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    finalData3=finalData3 %>%
        dplyr::rename(
            `doublet G1+G2`=g1G2Doublet,
            `doublet G1+G2 count`=g1G2DoubletCount,
            `doublet G2+G2`=g2G2Doublet,
            `doublet G2+G2 count`=g2G2DoubletCount
        )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
        flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
        doubletRSE<-.popConfidenceDoublet(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        doubletRSE<-finalData3 %>% dplyr::mutate(
            residualDoublet=NA
        )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
        twoPopRSE<-.popConfidence2Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
            residual2Pop=NA
        )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
        threePopRSE<-.popConfidence3Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        threePopRSE<-finalData3 %>% dplyr::mutate(
            residual3Pop=NA
        )
    }
    
    finalData4<-sqldf::sqldf(
        "select ds.*,
            ds2.residual,
            ds3.residualDoublet,
            ds4.residual2Pop,
            ds5.residual3Pop
            from finalData3 ds
            left join initialRSE ds2
            on ds.data = ds2.data
            left join doubletRSE ds3
            on ds.data = ds3.data
            left join twoPopRSE ds4
            on ds.data = ds4.data
            left join threePopRSE ds5
            on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
        residualMultiple=ifelse(
            is.na(residual3Pop), residual2Pop, residual3Pop)
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
        finalRSE=do.call(
            pmin, 
            c(
                subset(
                    ., 
                    select=c(
                    residual, 
                    residualDoublet, 
                    residualMultiple
                    )
                ),
                na.rm=TRUE
            )
        )
    ) %>% dplyr::rename(
        singleRSE=residual,
        doubletRSE=residualDoublet,
        multipleRSE=residualMultiple
    ) %>% dplyr::select(
        -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
        if(
            !is.na(finalData5[i,]$multipleRSE) &
            finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
        ){
            finalData5[i,]$investigate = 1
        }
        
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
        dplyr::select(-Success) %>% 
        dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
        dplyr::mutate(id=rowid(Data)) %>%
        tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
        Algorithm=do.call(
            pmax, 
            c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
        dplyr::select(Data, Algorithm) %>%
        dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
        "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
        finalData6=finalData6 %>%
            dplyr::select(
                -c(
                    `doublet G1+G2`,
                    `doublet G1+G2 count`,
                    `doublet G2+G2`,
                    `doublet G2+G2 count`
                )
            )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)
    write.csv(
        finalData6,
        paste0(file.path(dirname(getwd()), subDir), "/ploidyPeaksOutput.csv")
    )
    
    }else if(
        purrr::is_empty(investigateDs) &
        !purrr::is_empty(finishedDs) &
        purrr::is_empty(singleDs)
    ){
    
    ##Formatting the data from the other peak algorithms
    finalPart1<-finishedDs %>% data.frame() %>%
        dplyr::select(
            -c(
                "g3LL",
                "g3UL",
                "g4LL",
                "g4UL",
                "propCellsUsed",
                "cluster",
                "distToNext",
                "LL",
                "UL"
            )
        ) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    ##Removing columns
    finalData<-finalPart1 %>% dplyr::select(
        c(
            "data",
            "G1",
            "G2",
            "g1G2Doublet",
            "g2G2Doublet",
            "G1Count",
            "G2Count",
            "g1G2DoubletCount",
            "g2G2DoubletCount",
            "investigate"
        )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
        dplyr::select(data, G1, G1Count, G2, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
        dplyr::select(-c(G1, G1Count, G2, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
        dplyr::mutate(id=rowid(data)) %>%
        tidyr::pivot_wider(
            names_from=id,
            values_from=c(G1, G1Count, G2, G2Count)
        )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    finalData3=finalData3 %>%
        dplyr::rename(
            `doublet G1+G2`=g1G2Doublet,
            `doublet G1+G2 count`=g1G2DoubletCount,
            `doublet G2+G2`=g2G2Doublet,
            `doublet G2+G2 count`=g2G2DoubletCount
        )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
        flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
        doubletRSE<-.popConfidenceDoublet(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        doubletRSE<-finalData3 %>% dplyr::mutate(
            residualDoublet=NA
        )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
        twoPopRSE<-.popConfidence2Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
            residual2Pop=NA
        )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
        threePopRSE<-.popConfidence3Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        threePopRSE<-finalData3 %>% dplyr::mutate(
            residual3Pop=NA
        )
    }
    
    finalData4<-sqldf::sqldf(
        "select ds.*,
            ds2.residual,
            ds3.residualDoublet,
            ds4.residual2Pop,
            ds5.residual3Pop
            from finalData3 ds
            left join initialRSE ds2
            on ds.data = ds2.data
            left join doubletRSE ds3
            on ds.data = ds3.data
            left join twoPopRSE ds4
            on ds.data = ds4.data
            left join threePopRSE ds5
            on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
        residualMultiple=ifelse(
            is.na(residual3Pop),
            residual2Pop,
            residual3Pop
        )
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
        finalRSE=do.call(
            pmin, 
            c(
                subset(
                ., 
                select=c(residual, residualDoublet, residualMultiple)
                ), na.rm=TRUE
            )
        )
        ) %>% dplyr::rename(
            singleRSE=residual,
            doubletRSE=residualDoublet,
            multipleRSE=residualMultiple
        ) %>% dplyr::select(
            -c(residual2Pop, residual3Pop)
        )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
        if(
            !is.na(finalData5[i,]$multipleRSE) &
            finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
        ){
            finalData5[i,]$investigate = 1
        }
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
        dplyr::select(-Success) %>% 
        dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
        dplyr::mutate(id=rowid(Data)) %>%
        tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
        Algorithm=do.call(
            pmax, 
            c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
        dplyr::select(Data, Algorithm) %>%
        dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
        "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
        finalData6=finalData6 %>%
            dplyr::select(
                -c(
                    `doublet G1+G2`,
                    `doublet G1+G2 count`,
                    `doublet G2+G2`,
                    `doublet G2+G2 count`
                )
            )
    }
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings = FALSE)
    write.csv(
        finalData6,
        paste0(file.path(dirname(getwd()), subDir), "/ploidyPeaksOutput.csv")
    )
    }else{
    ##Formatting the diploid data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
        dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    ##Formatting the data from peak algorithm 2-4
    finalPart2=finishedDs %>% data.frame() %>%
        dplyr::select(
            -c(
                "g3LL",
                "g3UL",
                "g4LL",
                "g4UL",
                "propCellsUsed",
                "cluster",
                "distToNext",
                "LL",
                "UL"
            )
        ) %>%
        dplyr::mutate(
            investigate=0,
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    ##Formatting the data from the last peak algorthm
    finalPart3=investigateDs %>% data.frame() %>%
        dplyr::select(
            -c(
                "g3LL",
                "g3UL",
                "g4LL",
                "g4UL",
                "cluster",
                "distToNext",
                "LL",
                "UL"
            )
        ) %>%
        dplyr::mutate(
            G1=x,
            G1Count=y,
            G2=possiblePairX,
            G2Count=possiblePairY
        )
    
    ##Merging the three datasets
    finalData=rbind(
        finalPart1,
        finalPart2,
        finalPart3
    ) %>% dplyr::select(
        c(
            "data",
            "G1",
            "G2",
            "g1G2Doublet",
            "g2G2Doublet",
            "G1Count",
            "G2Count",
            "g1G2DoubletCount",
            "g2G2DoubletCount",
            "investigate"
        )
    )
    
    finalData$data=factor(finalData$data)
    finalData=finalData[order(finalData$data), ]
    
    ##Creating the doublet indicator
    finalData$doublet=0
    finalData$doublet[!is.na(finalData$g1G2Doublet)]=1
    finalData$g2G2Doublet[finalData$doublet == 0]=NA
    finalData$g2G2DoubletCount[finalData$doublet == 0]=NA
    
    ##Turning the long dataset into a wide dataset
    ##Selection the variables that need to be transformed from long to wide
    finalDataG1G2=finalData %>%
        dplyr::select(data, G1, G2, G1Count, G2Count)
    ##Selection the variables that do not need to change
    finalDataFlags=finalData %>%
        dplyr::select(-c(G1, G2,G1Count, G2Count))
    
    finalDataFlags=finalDataFlags[!duplicated(finalDataFlags$data), ]
    
    ##Making wide dataset
    finalData2=finalDataG1G2 %>%
        dplyr::mutate(id=rowid(data)) %>%
        tidyr::pivot_wider(
            names_from=id,
            values_from=c(G1, G2,G1Count, G2Count)
        )
    
    ##Merging data together
    finalData3=merge(finalData2, finalDataFlags, by="data")
    
    finalData3=finalData3 %>%
        dplyr::rename(
            `doublet G1+G2`=g1G2Doublet,
            `doublet G1+G2 count`=g1G2DoubletCount,
            `doublet G2+G2`=g2G2Doublet,
            `doublet G2+G2 count`=g2G2DoubletCount
        )
    
    ##Getting RSE value
    ##Testing hypothesis of a single population or multiple population
    initialRSE=.popConfidenceInitial(
        flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    if( TRUE %in% finalData3$doublet ){
        doubletRSE<-.popConfidenceDoublet(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        doubletRSE<-finalData3 %>% dplyr::mutate(
            residualDoublet=NA
        )
    }
    
    if(TRUE %in% grepl("_2", names(finalData3))){
        twoPopRSE<-.popConfidence2Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        twoPopRSE<-finalData3 %>% dplyr::mutate(
            residual2Pop=NA
        )
    }
    
    if(TRUE %in% grepl("_3", names(finalData3))){
        threePopRSE<-.popConfidence3Pop(
            flowDir, ds=finalData3, xVariable, saveGraph
        )
    }else{
        threePopRSE<-finalData3 %>% dplyr::mutate(
            residual3Pop=NA
        )
    }
    
    finalData4<-sqldf::sqldf(
        "select ds.*,
            ds2.residual,
            ds3.residualDoublet,
            ds4.residual2Pop,
            ds5.residual3Pop
            from finalData3 ds
            left join initialRSE ds2
            on ds.data = ds2.data
            left join doubletRSE ds3
            on ds.data = ds3.data
            left join twoPopRSE ds4
            on ds.data = ds4.data
            left join threePopRSE ds5
            on ds.data = ds5.data"
    )
    
    finalData4<-finalData4 %>% dplyr::mutate(
        residualMultiple=ifelse(
            is.na(residual3Pop),
            residual2Pop,
            residual3Pop
        )
    )
    
    finalData5<-finalData4 %>% dplyr::mutate(
        finalRSE=do.call(
            pmin, 
            c(
                subset(
                    .,
                    select=c(residual, residualDoublet, residualMultiple)
                ), na.rm=TRUE
            )
        )
    ) %>% dplyr::rename(
        singleRSE=residual,
        doubletRSE=residualDoublet,
        multipleRSE=residualMultiple
    ) %>% dplyr::select(
        -c(residual2Pop, residual3Pop)
    )
    
    ##Flagging to investigate if the algorithm found multiple subpopulations
    ##But the RSE for single population is less than RSE for multiple subpop
    for(i in seq_len(nrow(finalData5))){
        if(
            !is.na(finalData5[i,]$multipleRSE) &
            finalData5[i,]$multipleRSE > finalData5[i,]$singleRSE
        ){
            finalData5[i,]$investigate = 1
        } 
    }
    
    ##adding which algorithm analyzed the flow frame
    logDs2<-.GlobalEnv$logDs %>% 
        dplyr::select(-Success) %>% 
        dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3<-logDs2 %>%
        dplyr::mutate(id=rowid(Data)) %>%
        tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4<-logDs3 %>% dplyr::mutate(
        Algorithm=do.call(
            pmax, 
            c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5<-logDs4 %>%
        dplyr::select(Data, Algorithm) %>%
        dplyr::rename(data=Data)
    
    ##merging the algorithm used with the final dataset
    finalData6<-merge(finalData5, logDs5, by="data")
    finalData6<-finalData6 %>% dplyr::rename(
        "Sample" = "data"
    )
    
    ##If doublet = FALSE, remove doublet information
    if(doubletFlag == FALSE){
        finalData6=finalData6 %>%
            dplyr::select(
                -c(
                    `doublet G1+G2`,
                    `doublet G1+G2 count`,
                    `doublet G2+G2`,
                    `doublet G2+G2 count`
                )
            )
    }
    
    
    ##Creating a folder called analysis where the dataset will be saved
    setwd(flowDir)
    subDir <- "analysis"
    dir.create(file.path(dirname(getwd()), subDir), showWarnings=FALSE)
    write.csv(
        finalData6,
        paste0(file.path(dirname(getwd()), subDir), "/ploidyPeaksOutput.csv")
    )
    
    }
    
}

## popConfidenceInitial
.popConfidenceInitial = function(flowDir, ds, xVariable, saveGraph = TRUE){
    
    ##Removing NOTE 'no visible binding for global variable'
    G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
    
    modelData <- ds %>% dplyr::select(
        data,
        G1_1,
        G2_1,
        G1Count_1,
        G2Count_1
    )
    
    modelData02 <- modelData %>% dplyr::mutate(
        sdG1Count=G1Count_1*0.6,
        sdG2Count=G2Count_1*0.6
    )
    
    flowNameDs <- unique(modelData02$data)
    
    residualDs <- data.frame(
        matrix(nrow=0, ncol=14)
    )
    
    colnames(residualDs) <- c(
        "data",
        "G1_1",
        "G2_1",
        "G1Count_1",
        "G2Count_1",
        "sdG1Count",
        "sdG2Count",
        "g1Mean",
        "g2Mean",
        "g1SD",
        "g2SD",
        "numG1",
        "numG2",
        "residual"
    )
    
    for(k in seq_len(length(flowNameDs))){
        flowName <- flowCore::read.FCS(
            paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
        )
        flowData <- .smoothData( flowName, xVariable, 5)
        flowDataMeans <- modelData02%>% dplyr::filter(
            data == flowNameDs[k]
        )
        
        if(0 %in% flowData$x){
            flowData <- flowData %>% dplyr::filter(x != 0)
        }
        
        if(flowDataMeans$G2_1 != 0){
            midPoint <- flowData[
                which(
                    abs(
                        flowData$x-mean(
                        c(flowDataMeans$G1_1, flowDataMeans$G2_1))
                    ) == min(
                        abs(flowData$x-mean(
                        c(flowDataMeans$G1_1, flowDataMeans$G2_1)))
                    )
                ),
            ]
            
            if(nrow(midPoint)>1){
                midPoint <- midPoint[nrow(midPoint), ]
            }
            
            #G1 standard deviation 
            g1LeftFlowData <- flowData[
                seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
            g1RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G1_1)+1):
                which(flowData$x == midPoint$x)-1,
            ]
            
            leftPeak1 <- g1LeftFlowData[
                which(
                    abs(
                        g1LeftFlowData$y-flowDataMeans$sdG1Count
                    ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
                ), 
            ]
            
            if(nrow(leftPeak1)>1){
                leftPeak1 <- leftPeak1[1, ]
            }
            
            rightPeak1 <- g1RightFlowData[
                which(
                    abs(
                        g1RightFlowData$y-flowDataMeans$sdG1Count
                    ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            if(nrow(rightPeak1)>1 ){
                rightPeak1 <- rightPeak1[1, ]
            }
            
            #G2 standard deviation
            g2LeftFlowData <- flowData[
                which(flowData$x == midPoint$x):
                which(flowData$x == flowDataMeans$G2_1)-1,
            ]
            g2RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G2_1)+1):
                nrow(flowData),
            ]
            
            leftPeak2 <- g2LeftFlowData[
                which(
                    abs(
                        g2LeftFlowData$y-flowDataMeans$sdG2Count
                    ) == min(abs(g2LeftFlowData$y-flowDataMeans$sdG2Count)
                    )
                ),
            ]
            
            if(nrow(leftPeak2)>1 ){
                leftPeak2 <- leftPeak2[1, ]
            }
            
            rightPeak2 <- g2RightFlowData[
                which(
                    abs(
                        g2RightFlowData$y-flowDataMeans$sdG2Count
                    ) == min(abs(g2RightFlowData$y-flowDataMeans$sdG2Count)
                    )
                ),
            ]
            
            if(nrow(rightPeak2)>1 ){
                rightPeak2 <- rightPeak2[1, ]
            }
            
            flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
                g1Mean=G1_1,
                g2Mean=G2_1,
                g1SD=mean(
                    c(
                        rightPeak1$x-flowDataMeans$G1_1,
                        flowDataMeans$G1_1-leftPeak1$x
                    )
                ),
                g2SD=mean(
                    c(
                        rightPeak2$x-flowDataMeans$G2_1,
                        flowDataMeans$G2_1-leftPeak2$x
                    )
                ),
                numG1=sum(
                    flowData$y[
                        c(
                            which(flowData$x == leftPeak1$x):
                            which(flowData$x == rightPeak1$x)
                        )
                    ]
                ),
                numG2=sum(
                    flowData$y[
                        c(
                            which(flowData$x == leftPeak2$x):
                            which(flowData$x == rightPeak2$x)
                        )
                    ]
                )
            )
            g1Mean <- flowDataMeans01$g1Mean
            g2Mean <- flowDataMeans01$g2Mean
            g1SD <- flowDataMeans01$g1SD
            g2SD <- flowDataMeans01$g2SD
            
            xpectr::suppress_mw(
                singlePopNLS <- nls(
                    formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                    exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
                    (N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
                    (A + B*x + C*(x^2))*
                    (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                    (2*(g1SD*(x/g1Mean))^2))),
                    data=flowData,
                    start=c(
                        N1=flowDataMeans01$numG1,
                        N2=flowDataMeans01$numG2,
                        A=0, B=0, C=0
                    ),
                    nls.control(warnOnly=TRUE)
                )
            )
        }else{
            #G1 standard deviation 
            g1LeftFlowData <- flowData[
                seq_len(which(flowData$x == flowDataMeans$G1_1)-1),
            ]
            g1RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G1_1)+1):nrow(flowData),
            ]
            
            leftPeak1 <- g1LeftFlowData[
                which(
                    abs(
                        g1LeftFlowData$y-flowDataMeans$sdG1Count
                    ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            rightPeak1 <- g1RightFlowData[
                which(
                    abs(
                        g1RightFlowData$y-flowDataMeans$sdG1Count
                    ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
                g1Mean=G1_1,
                g2Mean=G2_1,
                g1SD=mean(
                    c(
                        rightPeak1$x-flowDataMeans$G1_1,
                        flowDataMeans$G1_1-leftPeak1$x
                    )
                ),
                g2SD=NA,
                numG1=sum(
                    flowData$y[
                        c(
                            which(flowData$x == leftPeak1$x):
                            which(flowData$x == rightPeak1$x)
                        )
                    ]
                ),
                numG2=0
            )
            
            g1Mean <- flowDataMeans01$g1Mean
            g1SD <- flowDataMeans01$g1SD
            
            xpectr::suppress_mw(
                singlePopNLS <- nls(
                    formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                    exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
                    (A + B*x + C*(x^2))*
                    (1/(sqrt(2 * pi)*g1SD*(x/g1Mean))*
                    exp(-((x-g1Mean)^2)/(2*(g1SD*(x/g1Mean))^2))),
                    data=flowData,
                    start=c(
                        N1=flowDataMeans01$numG1,
                        A=0, B=0, C=0
                    ),
                    nls.control(warnOnly=TRUE)
                )
            )
        }
        
        if(saveGraph == TRUE){
            nlsGraph <- ggplot()+
            geom_line(data=flowData, aes(x=x, y=y, col ='Raw Data'), size=1)+
            geom_line(aes(flowData$x,predict(singlePopNLS), col='NLS'))+
            labs(y="Counts", x=xVariable)+
            theme(legend.title=element_blank())
            
            plotDir <- "nlsGraphs"
            dir.create(file.path(dirname(flowDir), plotDir), showWarnings=FALSE)
            plotInitDir <- "nlsSingle_graphs"
            dir.create(
                file.path(dirname(flowDir) ,plotDir, plotInitDir),
                showWarnings=FALSE
            )
            plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
            
            png(
                paste0(plotOutFile, "/", flowNameDs[k], '.png'),
                width=600, height=400
            )
            print(nlsGraph)
            dev.off()
        }
        
        flowDataMeans01$residual <- summary(singlePopNLS)[["sigma"]]
        
        residualDs <- rbind(
            residualDs,
            flowDataMeans01
        )
        
    }
    
    return(residualDs)
}

## popConfidenceDoublet
.popConfidenceDoublet = function(flowDir, ds, xVariable, saveGraph = TRUE){
    ##Removing NOTE 'no visible binding for global variable'
    doublet<-G1_1<-G2_1<-`doublet G1+G2`<-`doublet G2+G2`<-G1Count_1<-NULL
    G2Count_1<-`doublet G1+G2 count`<-`doublet G2+G2 count`<-x<-y<-NULL
    
    ds2 <- ds %>% dplyr::filter(doublet == 1)
    flowNameDs <- unique(ds2$data)
    
    modelData01 <- ds2 %>% dplyr::select(
        data,
        G1_1,
        G2_1,
        `doublet G1+G2`,
        `doublet G2+G2`,
        G1Count_1,
        G2Count_1,
        `doublet G1+G2 count`,
        `doublet G2+G2 count`,
        doublet
    )
    
    modelData02 <- modelData01 %>% dplyr::mutate(
        sdG1Count=G1Count_1*0.6,
        sdG2Count=G2Count_1*0.6,
        sdG1G2Count=`doublet G1+G2 count`*0.6,
        sdG2G2Count=`doublet G2+G2 count`*0.6
    )
    
    residualDoubletDs <- data.frame(
        matrix(nrow=0, ncol=26)
    )
    colnames(residualDoubletDs) <- c(
        "data",
        "G1_1",
        "G2_1",
        "doublet G1+G2",
        "doublet G2+G2",
        "G1Count_1",
        "G2Count_1",
        "doublet G1+G2 count",
        "doublet G2+G2 count",
        "sdG1Count",
        "sdG2Count",
        "sdG1G2Count",
        "sdG2G2Count",
        "g1Mean",
        "g2Mean",
        "G1G2Mean",
        "G2G2Mean",
        "g1SD",
        "g2SD",
        "G1g2SD",
        "G2g2SD",
        "numG1",
        "numG2",
        "numDoublet1",
        "numDoublet2",
        "residualDoublet"
    )
    
    
    for(k in seq_len(length(flowNameDs))){
        flowName <- flowCore::read.FCS(
            paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
        )
        flowData <- .smoothData( flowName, xVariable, 5)
        flowDataMeans <- modelData02%>% dplyr::filter(
            data == flowNameDs[k]
        )
        
        if(0 %in% flowData$x){
            flowData <- flowData %>% dplyr::filter(x != 0)
        }
        
        if(flowDataMeans$G2_1 != 0){
            midPoint <- flowData[
                which(
                    abs(
                        flowData$x-mean(c(flowDataMeans$G1_1,
                                        flowDataMeans$G2_1))
                    ) == min(
                        abs(
                        flowData$x-mean(c(flowDataMeans$G1_1,
                                        flowDataMeans$G2_1))
                        )
                    )
                ),
            ]
            
            if(nrow(midPoint)>1){
                midPoint <- midPoint[nrow(midPoint), ]
            }
            
            #G1 standard deviation 
            G1LeftFlowData <- flowData[
                seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
            G1RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G1_1)+1):
                which(flowData$x == midPoint$x),
            ]
            
            leftPeak1 <- G1LeftFlowData[
                which(
                    abs(G1LeftFlowData$y-flowDataMeans$sdG1Count) == 
                    min(abs(G1LeftFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            if(nrow(leftPeak1)>1 ){
                leftPeak1 <- leftPeak1[1, ]
            }
            
            rightPeak1 <- G1RightFlowData[
                which(
                    abs(G1RightFlowData$y-flowDataMeans$sdG1Count) == 
                    min(abs(G1RightFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            if(nrow(rightPeak1)>1 ){
                rightPeak1 <- rightPeak1[1, ]
            }
            
            #G2 standard deviation
            G2LeftFlowData <- flowData[
                which(flowData$x == midPoint$x):
                (which(flowData$x == flowDataMeans$G2_1)-1),
            ]
            G2RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G2_1)+1):nrow(flowData),
            ]
            
            leftPeak2 <- G2LeftFlowData[
                which(
                    abs(G2LeftFlowData$y-flowDataMeans$sdG2Count) == 
                    min(abs(G2LeftFlowData$y-flowDataMeans$sdG2Count))
                ),
            ]
            
            if(nrow(leftPeak2)>1 ){
                leftPeak2 <- leftPeak2[1, ]
            }
            
            rightPeak2 <- G2RightFlowData[
                which(
                    abs(G2RightFlowData$y-flowDataMeans$sdG2Count) == 
                    min(abs(G2RightFlowData$y-flowDataMeans$sdG2Count))
                ),
            ]
            
            if(nrow(rightPeak2)>1 ){
                rightPeak2 <- rightPeak2[1, ]
            }
            
            #G1+G2 doublet standard deviation
            midPointDoublet <- flowData[
                which(
                    abs(
                        flowData$x-mean(c(flowDataMeans$`doublet G1+G2`,
                                        flowDataMeans$G2_1))
                    ) == min(
                        abs(
                        flowData$x-mean(c(flowDataMeans$`doublet G1+G2`,
                                        flowDataMeans$G2_1))
                        )
                    )
                ),
            ]
            
            if(nrow(midPointDoublet)>1){
                midPointDoublet <- midPointDoublet[nrow(midPointDoublet), ]
            }
            
            doubletG1G2LeftFlowData <- flowData[
                which(flowData$x == midPointDoublet$x):
                (which(flowData$x == flowDataMeans$`doublet G1+G2`)-1),
            ]
            doubletG1G2RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$`doublet G1+G2`)+1):
                nrow(flowData),
            ]
            
            doubletG1G2Left <- doubletG1G2LeftFlowData[
                which(
                    abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count) == 
                        min(
                        abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count)
                    )
                ),
            ]
            
            if(nrow(doubletG1G2Left)>1){
                doubletG1G2Left <- doubletG1G2Left[1, ]
            }
            
            doubletG1G2Right <- doubletG1G2RightFlowData[
                which(
                    abs(
                        doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count
                    ) == min(
                        abs(doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count)
                    )
                ),
            ]
            
            if(nrow(doubletG1G2Right)>1){
                doubletG1G2Right <- doubletG1G2Right[1, ]
            }
            
            if(!is.na(flowDataMeans$`doublet G2+G2`)){
                #G2+G2 doublet standard deviation
                midPointDoublet2 <- flowData[
                    which(
                        abs(
                            flowData$x-mean(
                                c(
                                    flowDataMeans$`doublet G1+G2`,
                                    flowDataMeans$`doublet G2+G2`
                                )
                            )
                        ) == min(
                            abs(
                                flowData$x-mean(
                                    c(
                                        flowDataMeans$`doublet G1+G2`,
                                        flowDataMeans$`doublet G2+G2`
                                    )
                                )
                            )
                        )
                    ),
                ]
                
                if(nrow(midPointDoublet2)>1){
                    midPointDoublet2 <- midPointDoublet2[
                    nrow(midPointDoublet2), ]
                }
                
                doubletG1G2LeftFlowData <- flowData[
                    which(flowData$x == midPointDoublet2$x):
                    (which(flowData$x == flowDataMeans$`doublet G2+G2`)-1),
                ]
                
                doubletG2G2RightFlowData <- flowData[
                    (which(flowData$x == flowDataMeans$`doublet G2+G2`)+1):
                    nrow(flowData),
                ]
                
                doubletG2G2Left <- doubletG1G2LeftFlowData[
                    which(
                        abs(
                            doubletG1G2LeftFlowData$y-
                            flowDataMeans$sdG2G2Count
                        ) == min(
                            abs(
                                doubletG1G2LeftFlowData$y-
                                flowDataMeans$sdG2G2Count
                            )
                        )
                    ),
                ]
                
                if(nrow(doubletG2G2Left)>1){
                    doubletG2G2Left <- doubletG2G2Left[1, ]
                }
                
                doubletG2G2Right <- doubletG2G2RightFlowData[
                    which(
                        abs(
                            doubletG2G2RightFlowData$y-flowDataMeans$sdG2G2Count
                        ) == min(
                            abs(
                                doubletG2G2RightFlowData$y-
                                flowDataMeans$sdG2G2Count
                            )
                        )
                    ),
                ]
                
                if(nrow(doubletG2G2Right)>1){
                    doubletG2G2Right <- doubletG2G2Right[1, ]
                }
                flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
                    g1Mean=G1_1,
                    g2Mean=G2_1,
                    G1G2Mean=`doublet G1+G2`,
                    G2G2Mean=`doublet G2+G2`,
                    g1SD=mean(
                        c(
                            rightPeak1$x-flowDataMeans$G1_1,
                            flowDataMeans$G1_1-leftPeak1$x
                        )
                    ),
                    g2SD=mean(
                        c(
                            rightPeak2$x-flowDataMeans$G2_1,
                            flowDataMeans$G2_1-leftPeak2$x
                        )
                    ),
                    G1g2SD=mean(
                        c(
                            doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
                            flowDataMeans$`doublet G1+G2` - doubletG1G2Left$x
                        )
                    ),
                    G2g2SD=mean(
                        c(
                            doubletG2G2Right$x-flowDataMeans$`doublet G2+G2`,
                            flowDataMeans$`doublet G2+G2` - doubletG2G2Left$x
                        )
                    ),
                    numG1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == leftPeak1$x):
                                which(flowData$x == rightPeak1$x)
                            )
                        ]
                    ),
                    numG2=sum(
                        flowData$y[
                            c(
                                which(flowData$x == leftPeak2$x):
                                which(flowData$x == rightPeak2$x)
                            )
                        ]
                    ),
                    numDoublet1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == doubletG1G2Left$x):
                                which(flowData$x == doubletG1G2Right$x)
                            )
                        ]
                    ),
                    numDoublet2=sum(
                        flowData$y[
                            c(
                                which(flowData$x == doubletG2G2Left$x):
                                which(flowData$x == doubletG2G2Right$x)
                            )
                        ]
                    )
                )
                
                g1Mean <- flowDataMeans01$g1Mean
                g2Mean <- flowDataMeans01$g2Mean
                G1G2Mean <- flowDataMeans01$G1G2Mean
                G2G2Mean <- flowDataMeans01$G2G2Mean
                
                g1SD <- flowDataMeans01$g1SD
                g2SD <- flowDataMeans01$g2SD
                G1g2SD <- flowDataMeans01$G1g2SD
                G2g2SD <- flowDataMeans01$G2g2SD
                
                xpectr::suppress_mw(
                    singlePopNLS <- nls(
                        formula = y ~ (N1/(sqrt(2*pi) * g1SD) * 
                        exp(-((x-g1Mean)^2)/(2 *g1SD^2))) +
                        (A1 + B1*x + C1*(x^2))*
                        (1/(sqrt(2 * pi) * g1SD * (x/g1Mean)) * 
                        exp(-((x - g1Mean)^2)/(2 *(g1SD* (x/g1Mean))^2)))+
                        (N2/(sqrt(2 * pi) * g2SD) * exp(-((x - g2Mean)^2)/
                        (2 *g2SD^2)))+
                        (A2 + B2*x + C2*(x^2))*
                        (1/(sqrt(2 * pi) * g2SD * (x/g2Mean)) * 
                        exp(-((x - g2Mean)^2)/(2 *(g2SD* (x/g2Mean))^2)))+
                        (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
                        (2*G1g2SD^2)))+
                        (A3 + B3*x + C3*(x^2))*
                        (1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*
                        exp(-((x-G1G2Mean)^2)/(2*(G1g2SD*(x/G1G2Mean))^2)))+
                        (numDoublet2/(sqrt(2*pi)*G2g2SD)*exp(-((x-G2G2Mean)^2)/
                        (2 *G2g2SD^2))),
                        data=flowData,
                        start=c(
                            N1=flowDataMeans01$numG1,
                            N2=flowDataMeans01$numG2,
                            numDoublet1=flowDataMeans01$numDoublet1,
                            numDoublet2=flowDataMeans01$numDoublet2,
                            A1=0, B1=0, C1=0,
                            A2=0, B2=0, C2=0,
                            A3=0, B3=0, C3=0
                        ),
                        nls.control(warnOnly=TRUE)
                    )
                )
            }else{
                flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
                    g1Mean=G1_1,
                    g2Mean=G2_1,
                    G1G2Mean=`doublet G1+G2`,
                    G2G2Mean=`doublet G2+G2`,
                    g1SD=mean(
                        c(
                            rightPeak1$x-flowDataMeans$G1_1,
                            flowDataMeans$G1_1-leftPeak1$x
                        )
                    ),
                    g2SD=mean(
                        c(
                            rightPeak2$x-flowDataMeans$G2_1,
                            flowDataMeans$G2_1-leftPeak2$x
                        )
                    ),
                    G1g2SD=mean(
                        c(
                            doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
                            flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x
                        )
                    ),
                    G2g2SD=NA,
                    numG1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == leftPeak1$x):
                                which(flowData$x == rightPeak1$x)
                            )
                        ]
                    ),
                    numG2=sum(
                        flowData$y[
                            c(
                                which(flowData$x == leftPeak2$x):
                                which(flowData$x == rightPeak2$x)
                            )
                        ]
                    ),
                    numDoublet1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == doubletG1G2Left$x):
                                which(flowData$x == doubletG1G2Right$x)
                            )
                        ]
                    ),
                    numDoublet2=0
                )
                
                g1Mean <- flowDataMeans01$g1Mean
                g2Mean <- flowDataMeans01$g2Mean
                G1G2Mean <- flowDataMeans01$G1G2Mean
                g1SD <- flowDataMeans01$g1SD
                g2SD <- flowDataMeans01$g2SD
                G1g2SD <- flowDataMeans01$G1g2SD
                
                xpectr::suppress_mw(
                    singlePopNLS <- nls(
                        formula=y ~ (N1/(sqrt(2*pi)*g1SD)*
                        exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
                        (A1 + B1*x + C1*(x^2))*
                        (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                        (2*(g1SD*(x/g1Mean))^2)))+
                        (N2/(sqrt(2*pi)*g2SD)*exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
                        (A2 + B2*x + C2*(x^2))*
                        (1/(sqrt(2*pi)*g2SD*(x/g2Mean))*exp(-((x-g2Mean)^2)/
                        (2*(g2SD*(x/g2Mean))^2)))+
                        (numDoublet1/(sqrt(2*pi)*G1g2SD)*
                        exp(-((x-G1G2Mean)^2)/(2*G1g2SD^2))),
                        data=flowData,
                        start=c(
                            N1=flowDataMeans01$numG1,
                            N2=flowDataMeans01$numG2,
                            numDoublet1=flowDataMeans01$numDoublet1,
                            A1=0, B1=0, C1=0,
                            A2=0, B2=0, C2=0
                        ),
                        nls.control(warnOnly=TRUE)
                    )
                )
                
            }
        }else{
            #G1 standard deviation 
            G1LeftFlowData <- flowData[
                seq_len(which(flowData$x == flowDataMeans$G1_1)), ]
            G1RightFlowData <- flowData[
                which(flowData$x == flowDataMeans$G1_1):nrow(flowData),
            ]
            
            leftPeak1 <- G1LeftFlowData[
                which(
                    abs(
                        G1LeftFlowData$y-flowDataMeans$sdG1Count
                    ) == min(abs(G1LeftFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            if(nrow(leftPeak1)>1 ){
                leftPeak1 <- leftPeak1[1, ]
            }
            
            rightPeak1 <- G1RightFlowData[
                which(
                    abs(
                        G1RightFlowData$y-flowDataMeans$sdG1Count
                    ) == min(abs(G1RightFlowData$y-flowDataMeans$sdG1Count))
                ),
            ]
            
            if(nrow(rightPeak1)>1 ){
                rightPeak1 <- rightPeak1[1, ]
            }
            
            
            #G1+G2 doublet standard deviation
            midPointDoublet <- flowData[
                which(
                    abs(
                        flowData$x-mean(
                        c(flowDataMeans$`doublet G1+G2`, flowDataMeans$G1_1)
                    )
                    ) == min(
                    abs(
                        flowData$x-mean(
                            c(
                                flowDataMeans$`doublet G1+G2`,
                                flowDataMeans$G1_1
                            )
                        )
                    )
                )
                ),
            ]
            
            if(nrow(midPointDoublet)>1){
                midPointDoublet <- midPointDoublet[nrow(midPointDoublet), ]
            }
            
            doubletG1G2LeftFlowData <- flowData[
                which(flowData$x == midPointDoublet$x):
                which(flowData$x == flowDataMeans$`doublet G1+G2`),
            ]
            
            doubletG1G2RightFlowData <- flowData[
                which(flowData$x == flowDataMeans$`doublet G1+G2`):
                nrow(flowData),
            ]
            
            doubletG1G2Left <- doubletG1G2LeftFlowData[
                which(
                    abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count) == 
                        min(
                        abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count)
                        )
                ),
            ]
            
            if(nrow(doubletG1G2Left)>1){
                doubletG1G2Left <- doubletG1G2Left[1, ]
            }
            
            doubletG1G2Right <- doubletG1G2RightFlowData[
                which(
                    abs(
                        doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count
                    ) == min(
                    abs(
                        doubletG1G2RightFlowData$y-
                        flowDataMeans$sdG1G2Count
                    )
                    )
                ),
            ]
            
            if(nrow(doubletG1G2Right)>1){
                doubletG1G2Right <- doubletG1G2Right[1, ]
            }
            
            if(!is.na(flowDataMeans$`doublet G2+G2`)){
                #G2+G2 doublet standard deviation
                midPointDoublet2 <- flowData[
                    which(
                        abs(
                            flowData$x-mean(
                                c(
                                    flowDataMeans$`doublet G1+G2`,
                                    flowDataMeans$`doublet G2+G2`)
                                )
                            ) == 
                            min(
                                abs(
                                flowData$x-mean(
                                    c(
                                        flowDataMeans$`doublet G1+G2`,
                                        flowDataMeans$`doublet G2+G2`
                                    )
                                    )
                                )
                            )
                    ),
                ]
                
                if(nrow(midPointDoublet2)>1){
                    midPointDoublet2 <- midPointDoublet2[
                        nrow(midPointDoublet2), ]
                }
                
                doubletG1G2LeftFlowData <- flowData[
                    which(flowData$x == midPointDoublet2$x):
                    which(flowData$x == flowDataMeans$`doublet G2+G2`), 
                ]
                doubletG2G2RightFlowData <- flowData[
                    which(flowData$x == flowDataMeans$`doublet G2+G2`):
                    nrow(flowData),
                ]
                
                doubletG2G2Left <- doubletG1G2LeftFlowData[
                    which(
                        abs(
                            doubletG1G2LeftFlowData$y-
                            flowDataMeans$sdG2G2Count
                        ) == min(
                            abs(
                                doubletG1G2LeftFlowData$y-
                                flowDataMeans$sdG2G2Count
                            )
                        )
                    ),
                ]
                
                if(nrow(doubletG2G2Left)>1){
                    doubletG2G2Left <- doubletG2G2Left[1, ]
                }
                
                doubletG2G2Right <- doubletG2G2RightFlowData[
                    which(
                        abs(
                            doubletG2G2RightFlowData$y-
                            flowDataMeans$sdG2G2Count
                        ) == min(
                            abs(
                                doubletG2G2RightFlowData$y-
                                flowDataMeans$sdG2G2Count
                            )
                        )
                    ),
                ]
                
                if(nrow(doubletG2G2Right)>1){
                    doubletG2G2Right <- doubletG2G2Right[1, ]
                }
                flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
                    g1Mean=G1_1,
                    g2Mean=G2_1,
                    G1G2Mean=`doublet G1+G2`,
                    G2G2Mean=`doublet G2+G2`,
                    g1SD=mean(
                        c(
                            rightPeak1$x-flowDataMeans$G1_1,
                            flowDataMeans$G1_1-leftPeak1$x
                        )
                    ),
                    g2SD=NA,
                    G1g2SD=mean(
                        c(
                            doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
                            flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x
                        )
                    ),
                    G2g2SD=mean(
                        c(
                            doubletG2G2Right$x-flowDataMeans$`doublet G2+G2`,
                            flowDataMeans$`doublet G2+G2`-doubletG2G2Left$x
                        )
                    ),
                    numG1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == leftPeak1$x):
                                which(flowData$x == rightPeak1$x)
                            )
                        ]
                    ),
                    numG2=0,
                    numDoublet1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == doubletG1G2Left$x):
                                which(flowData$x == doubletG1G2Right$x)
                            )
                        ]
                    ),
                    numDoublet2=sum(
                        flowData$y[
                            c(
                                which(flowData$x == doubletG2G2Left$x):
                                which(flowData$x == doubletG2G2Right$x)
                            )
                        ]
                    )
                )
                
                g1Mean <- flowDataMeans01$g1Mean
                G1G2Mean <- flowDataMeans01$G1G2Mean
                G2G2Mean <- flowDataMeans01$G2G2Mean
                g1SD <- flowDataMeans01$g1SD
                G1g2SD <- flowDataMeans01$G1g2SD
                G2g2SD <- flowDataMeans01$G2g2SD
                
                xpectr::suppress_mw(
                    singlePopNLS <- nls(
                        formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                        exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
                        (A1 + B1*x + C1*(x^2))*
                        (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x - g1Mean)^2)/
                        (2*(g1SD*(x/g1Mean))^2)))+
                        (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
                        (2*G1g2SD^2)))+
                        (A2 + B2*x + C2*(x^2))*
                        (1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*
                        exp(-((x-G1G2Mean)^2)/(2*(G1g2SD*(x/G1G2Mean))^2)))+
                        (numDoublet2/(sqrt(2*pi)*G2g2SD)*exp(-((x-G2G2Mean)^2)/
                        (2*G2g2SD^2))),
                        data=flowData,
                        start=c(
                            N1=flowDataMeans01$numG1,
                            numDoublet1=flowDataMeans01$numDoublet1,
                            numDoublet2=flowDataMeans01$numDoublet2,
                            A1=0, B1=0, C1=0,
                            A2=0, B2=0, C2=0
                        ),
                        nls.control(warnOnly=TRUE)
                    )
                )
            }else{
                flowDataMeans01 <- flowDataMeans %>% dplyr::mutate(
                    g1Mean=G1_1,
                    g2Mean=G2_1,
                    G1G2Mean=`doublet G1+G2`,
                    G2G2Mean=`doublet G2+G2`,
                    g1SD=mean(
                        c(
                            rightPeak1$x-flowDataMeans$G1_1,
                            flowDataMeans$G1_1-leftPeak1$x
                        )
                    ),
                    g2SD=NA,
                    G1g2SD=mean(
                        c(
                            doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
                            flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x
                        )
                    ),
                    G2g2SD=NA,
                    numG1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == leftPeak1$x):
                                which(flowData$x == rightPeak1$x)
                            )
                        ]
                    ),
                    numG2=0,
                    numDoublet1=sum(
                        flowData$y[
                            c(
                                which(flowData$x == doubletG1G2Left$x):
                                which(flowData$x == doubletG1G2Right$x)
                            )
                        ]
                    ),
                    numDoublet2=0
                )
                
                g1Mean <- flowDataMeans01$g1Mean
                G1G2Mean <- flowDataMeans01$G1G2Mean
                g1SD <- flowDataMeans01$g1SD
                G1g2SD <- flowDataMeans01$G1g2SD
                
                xpectr::suppress_mw(
                    singlePopNLS <- nls(
                        formula = y ~ (N1/(sqrt(2*pi)*g1SD)*
                        exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
                        (A1 + B1*x + C1*(x^2))*
                        (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                        (2*(g1SD*(x/g1Mean))^2)))+
                        (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
                        (2*G1g2SD^2))),
                        data=flowData,
                        start=c(
                            N1=flowDataMeans01$numG1,
                            numDoublet1=flowDataMeans01$numDoublet1,
                            A1=0, B1=0, C1=0
                        ),
                        nls.control(warnOnly=TRUE)
                    )
                )
            }
        }
        
        if(saveGraph == TRUE){
            nlsGraph <- ggplot()+
            geom_line(data=flowData, aes(x=x, y=y, col = 'Raw Data'), size=1)+
            geom_line(aes(flowData$x, predict(singlePopNLS), col ='NLS'))+
            labs(y ="Counts", x =xVariable)+
            theme(legend.title = element_blank())
            
            plotDir <- "nlsGraphs"
            plotInitDir <- "nlsDoublet_graphs"
            dir.create(
                file.path(dirname(flowDir), plotDir, plotInitDir),
                showWarnings=FALSE
            )
            plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
            
            png(
                paste0(plotOutFile, "/", flowNameDs[k], '.png'),
                width=600, height=400
            )
            print(nlsGraph)
            dev.off()
        }
        
        flowDataMeans01$residualDoublet <- summary(singlePopNLS)[["sigma"]]
        
        residualDoubletDs <- rbind(
            residualDoubletDs,
            flowDataMeans01
        )
        
    }
    return(residualDoubletDs)
}

## popConfidence2Pop
.popConfidence2Pop = function(flowDir, ds, xVariable, saveGraph = TRUE){
    
    ##Removing NOTE 'no visible binding for global variable'
    G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
    G1_2<-G2_2<-G1_3<-G1Count_2<-G2Count_2<-NULL
    
    if(TRUE %in% grepl("_3", names(ds))){
        ds1<-ds %>% dplyr::filter(is.na(G1_3))    
    }else{
        ds1<-ds
    }
    
    modelData <- ds1 %>% dplyr::select(
        data,
        G1_1,
        G2_1,
        G1_2,
        G2_2,
        G1Count_1,
        G2Count_1,
        G1Count_2,
        G2Count_2
    )
    
    modelData01<-modelData %>% dplyr::filter(
        !is.na(modelData$G1_2)
    )
    modelData02 <- modelData01 %>% dplyr::mutate(
        sdG1Count=G1Count_1*0.6,
        sdG2Count=G2Count_1*0.6,
        sdG1Count2=G1Count_2*0.6,
        sdG2Count2=G2Count_2*0.6
    )
    
    flowNameDs <- unique(modelData02$data)
    
    residualMultiDs <- data.frame(
        matrix(nrow=0, ncol=26)
    )
    
    colnames(residualMultiDs) <- c(
        "data",
        "G1_1",
        "G2_1",
        "G1_2",
        "G2_2", 
        "G1Count_1",
        "G2Count_1",
        "G1Count_2",
        "G2Count_2", 
        "sdG1Count",
        "sdG2Count",
        "sdG1Count2",
        "sdG2Count2",
        "g1Mean",
        "g2Mean",
        "g1SD",
        "g2SD",
        "numG1",
        "numG2",
        "g1Mean2",
        "g2Mean2",
        "g1SD2",
        "g2SD2",
        "pop2NumG1",
        "pop2NumG2",
        "residual"
    )
    
    
    for(k in seq_len(length(flowNameDs))){
        flowName <- flowCore::read.FCS(
            paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
        )
        flowData <- .smoothData( flowName, xVariable, 5)
        flowDataMeans <- modelData02%>% dplyr::filter(
            data == flowNameDs[k]
        )
        
        if(0 %in% flowData$x){
            flowData <- flowData %>% dplyr::filter(x != 0)
        }
        
        midPoint <- flowData[
            which(
                abs(
                    flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2))
                ) == min(
                    abs(flowData$x-mean(
                    c(flowDataMeans$G1_1,flowDataMeans$G1_2)
                    )
                    )
                )
            ),
        ]
        
        if(nrow(midPoint)>1){
            midPoint <- midPoint[nrow(midPoint), ]
        }
        
        #G1 standard deviation 
        g1LeftFlowData <- flowData[
            seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
        g1RightFlowData <- flowData[
            (which(flowData$x == flowDataMeans$G1_1)+1):
            which(flowData$x == midPoint$x)-1,
        ]
        
        leftPeak1 <- g1LeftFlowData[
            which(
                abs(
                    g1LeftFlowData$y-flowDataMeans$sdG1Count
                ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
            ),
        ]
        
        if(nrow(leftPeak1)>1){
            leftPeak1 <- leftPeak1[1, ]
        }
        
        rightPeak1 <- g1RightFlowData[
            which(
                abs(
                    g1RightFlowData$y-flowDataMeans$sdG1Count
                ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
            ),
        ]
        
        if(nrow(rightPeak1)>1 ){
            rightPeak1 <- rightPeak1[1, ]
        }
        midPointLeft <- flowData[
            which(
                abs(
                    flowData$x-mean(
                        c(flowDataMeans$G2_1, flowDataMeans$G1_2))
                    ) == min(
                        abs(
                        flowData$x-mean(
                            c(flowDataMeans$G2_1,flowDataMeans$G1_2))
                        )
                    )
            ),
        ]
        
        if(nrow(midPointLeft)>1){
            midPointLeft <- midPointLeft[nrow(midPointLeft), ]
        }
        
        midPointRight <- flowData[
            which(
                abs(
                    flowData$x-mean(
                        c(flowDataMeans$G2_1, flowDataMeans$G2_2))
                    ) == min(
                        abs(flowData$x-mean(
                            c(flowDataMeans$G2_1,flowDataMeans$G2_2))
                        )
                    )
            ),
        ]
        
        if(nrow(midPointRight)>1){
            midPointRight <- midPointRight[nrow(midPointRight), ]
        }
        
        #G2 standard deviation
        if(
            length(
                unique(
                    c(
                        flowDataMeans$G1_1,
                        flowDataMeans$G1_2,
                        flowDataMeans$G2_2,
                        flowDataMeans$G2_1 
                    )
                )
            ) == 3
        ){
            g2LeftFlowData <- flowData[
                which(flowData$x == midPoint$x):
                which(flowData$x == flowDataMeans$G2_1)-1,
            ]
            g2RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G2_1)+1):
                nrow(flowData),
            ]
            
            leftPeak2 <- g2LeftFlowData[
                which(
                    abs(g2LeftFlowData$y-flowDataMeans$sdG2Count) == 
                    min(abs(g2LeftFlowData$y-flowDataMeans$sdG2Count))
                ),
            ]
            
            if(nrow(leftPeak2)>1 ){
                leftPeak2 <- leftPeak2[1, ]
            }
        }else{
            g2LeftFlowData <- flowData[
                which(flowData$x == midPointLeft$x):
                which(flowData$x == flowDataMeans$G2_1)-1,
            ]
            g2RightFlowData <- flowData[
                (which(flowData$x == flowDataMeans$G2_1)+1):
                which(flowData$x == midPointRight$x),
            ]
            
            leftPeak2 <- g2LeftFlowData[
                which(
                    abs(g2LeftFlowData$y-flowDataMeans$sdG2Count) == 
                    min(abs(g2LeftFlowData$y-flowDataMeans$sdG2Count))
                ),
            ]
            
            if(nrow(leftPeak2)>1 ){
                leftPeak2 <- leftPeak2[1, ]
            }
        }
        
        rightPeak2 <- g2RightFlowData[
            which(
                abs(g2RightFlowData$y-flowDataMeans$sdG2Count) == 
                min(abs(g2RightFlowData$y-flowDataMeans$sdG2Count))
            ),
        ]
        
        if(nrow(rightPeak2)>1 ){
            rightPeak2 <- rightPeak2[1, ]
        }
        
        #pop2
        
        #G1 standard deviation
        g1LeftFlowData2 <- flowData[
            (which(flowData$x == midPoint$x)+1):
            which(flowData$x == flowDataMeans$G1_2)-1,
        ]
        
        g1RightFlowData2 <- flowData[
            (which(flowData$x == flowDataMeans$G1_2)+1):
            which(flowData$x == midPointLeft$x)-1,
        ]
        
        pop2leftPeak1 <- g1LeftFlowData2[
            which(
                abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2) == 
                min(abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2))
            ),
        ]
        
        if(nrow(pop2leftPeak1)>1){
            pop2leftPeak1 <- pop2leftPeak1[1, ]
        }
        
        pop2rightPeak1 <- g1RightFlowData2[
            which(
                abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2) == 
                min(abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2))
            ),
        ]
        
        if(nrow(pop2rightPeak1)>1 ){
            pop2rightPeak1 <- pop2rightPeak1[1, ]
        }
        
        #G2 standard deviation
        g2LeftFlowData2 <- flowData[
            which(flowData$x == midPointRight$x):
            which(flowData$x == flowDataMeans$G2_2)-1,
        ]
        g2RightFlowData2 <- flowData[
            (which(flowData$x == flowDataMeans$G2_2)+1):nrow(flowData),
        ]
        
        pop2leftPeak2 <- g2LeftFlowData2[
            which(
                abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2) == 
                min(abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2))
            ),
        ]
        
        if(nrow(pop2leftPeak2)>1 ){
            pop2leftPeak2 <- pop2leftPeak2[1, ]
        }
        
        pop2rightPeak2 <- g2RightFlowData2[
            which(
                abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2) == 
                min(abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2))
            ),
        ]
        
        if(nrow(pop2rightPeak2)>1 ){
            pop2rightPeak2 <- pop2rightPeak2[1, ]
        }
        
        flowDataMeans01 = flowDataMeans %>% dplyr::mutate(
            g1Mean=G1_1,
            g2Mean=G2_1,
            g1SD=mean(
                c(
                    rightPeak1$x-flowDataMeans$G1_1,
                    flowDataMeans$G1_1-leftPeak1$x
                )
            ),
            g2SD=mean(
                c(
                    rightPeak2$x-flowDataMeans$G2_1,
                    flowDataMeans$G2_1-leftPeak2$x
                )
            ),
            pop1NumG1=sum(
                flowData$y[
                    c(
                        which(flowData$x == leftPeak1$x):
                        which(flowData$x == rightPeak1$x)
                    )
                ]
            ),
            pop1NumG2=sum(
                flowData$y[
                    c(
                        which(flowData$x == leftPeak2$x):
                        which(flowData$x == rightPeak2$x)
                    )
                ]
            ),
            g1Mean2=G1_2,
            g2Mean2=G2_2,
            g1SD2=mean(
                c(
                    pop2rightPeak1$x-flowDataMeans$G1_2,
                    flowDataMeans$G1_2-pop2leftPeak1$x
                )
            ),
            g2SD2=mean(
                c(
                    pop2rightPeak2$x-flowDataMeans$G2_2,
                    flowDataMeans$G2_2-pop2leftPeak2$x
                )
            ),
            pop2NumG1=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop2leftPeak1$x):
                        which(flowData$x == pop2rightPeak1$x)
                    )
                ]
            ),
            pop2NumG2=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop2leftPeak2$x):
                        which(flowData$x == pop2rightPeak2$x)
                    )
                ]
            )
        )
        
        g1Mean <- flowDataMeans01$g1Mean
        g2Mean <- flowDataMeans01$g2Mean
        g1SD <- flowDataMeans01$g1SD
        g2SD <- flowDataMeans01$g2SD
        g1Mean2 <- flowDataMeans01$g1Mean2
        g2Mean2 <- flowDataMeans01$g2Mean2
        g1SD2 <- flowDataMeans01$g1SD2
        g2SD2 <- flowDataMeans01$g2SD2
        
        if(
            length(
                unique(
                    c(
                        flowDataMeans$G1_1,
                        flowDataMeans$G1_2,
                        flowDataMeans$G2_2,
                        flowDataMeans$G2_1 
                    )
                )
            ) == 3
        ){
        xpectr::suppress_mw(
            multiPopNLS <- nls(
                formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*
                exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
                (A + B*x + C*(x^2))*
                (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                (2*(g1SD*(x/g1Mean))^2)))+
                (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
                (A2 + B2*x + C2*(x^2))*
                (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                (2*(g2SD2*(x/g2Mean2))^2)))+
                (pop2N2/(sqrt(2*pi)*g2SD2)* exp(-((x-g2Mean2)^2)/(2*g2SD2^2))),
                data=flowData,
                start=c(
                    pop1N1=flowDataMeans01$pop1NumG1,
                    pop1N2=flowDataMeans01$pop1NumG2,
                    pop2N2=flowDataMeans01$pop2NumG2,
                    A=0, B=0, C=0,
                    A2=0, B2=0, C2=0
                ),
                nls.control(warnOnly=TRUE)
            )
        )
        }else{
            xpectr::suppress_mw(
                multiPopNLS <- nls(
                    formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*
                    exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
                    (A + B*x + C*(x^2))*
                    (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                    (2*(g1SD*(x/g1Mean))^2)))+
                    (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
                    (A2 + B2*x + C2*(x^2))*
                    (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                    (2*(g2SD2*(x/g2Mean2))^2)))+
                    (pop2N1/(sqrt(2*pi)*g1SD2)*exp(-((x-g1Mean2)^2)/(2*g1SD2^2)))+
                    (pop2N2/(sqrt(2*pi)*g2SD2)*exp(-((x-g2Mean2)^2)/(2*g2SD2^2))),
                    data=flowData,
                    start=c(
                        pop1N1=flowDataMeans01$pop1NumG1,
                        pop1N2=flowDataMeans01$pop1NumG2,
                        pop2N1=flowDataMeans01$pop2NumG1,
                        pop2N2=flowDataMeans01$pop2NumG2,
                        A=0, B=0, C=0,
                        A2=0, B2=0, C2=0
                    ),
                    nls.control(warnOnly=TRUE)
                )
            )
        }
        
        if(saveGraph == TRUE){
            nlsGraph <- ggplot()+
            geom_line(data=flowData, aes(x=x, y=y, col ='Raw Data'), size=1)+
            geom_line(aes(flowData$x,predict(multiPopNLS), col='NLS'))+
            labs(y="Counts", x=xVariable)+
            theme(legend.title=element_blank())
            
            plotDir <- "nlsGraphs"
            dir.create(file.path(dirname(flowDir), plotDir), showWarnings=FALSE)
            plotInitDir <- "nlsMultiple_graphs"
            dir.create(
                file.path(dirname(flowDir) ,plotDir, plotInitDir),
                showWarnings=FALSE
            )
            plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
            
            png(
                paste0(plotOutFile, "/", flowNameDs[k], '.png'),
                width=600, height=400
            )
            print(nlsGraph)
            dev.off()
        }
        
        flowDataMeans01$residual2Pop <- summary(multiPopNLS)[["sigma"]]
        
        residualMultiDs <- rbind(
            residualMultiDs,
            flowDataMeans01
        )
        
    }
    
    return(residualMultiDs)
}

## popConfidence3Pop
.popConfidence3Pop = function(flowDir, ds, xVariable, saveGraph = TRUE){
    
    ##Removing NOTE 'no visible binding for global variable'
    G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
    G1_2<-G2_2<-G1Count_2<-G2Count_2<-NULL
    G1_3<-G2_3<-G1Count_3<-G2Count_3<-NULL
    
    ds1<-ds %>% dplyr::filter(
        !is.na(ds$G1_3)
    )
    
    modelData <- ds1 %>% dplyr::select(
        data,
        G1_1,
        G2_1,
        G1_2,
        G2_2,
        G1_3,
        G2_3,
        G1Count_1,
        G2Count_1,
        G1Count_2,
        G2Count_2,
        G1Count_3,
        G2Count_3
    )
    
    modelData02 <- modelData %>% dplyr::mutate(
        sdG1Count=G1Count_1*0.6,
        sdG2Count=G2Count_1*0.6,
        sdG1Count2=G1Count_2*0.6,
        sdG2Count2=G2Count_2*0.6,
        sdG1Count3=G1Count_3*0.6,
        sdG2Count3=G2Count_3*0.6
    )
    
    flowNameDs <- unique(modelData02$data)
    
    residualMultiDs <- data.frame(
        matrix(nrow=0, ncol=38)
    )
    
    colnames(residualMultiDs) <- c(
        "data",
        "G1_1",
        "G2_1",
        "G1_2",
        "G2_2", 
        "G1_3", 
        "G2_3", 
        "G1Count_1",
        "G2Count_1",
        "G1Count_2",
        "G2Count_2", 
        "G1Count_3",
        "G2Count_3", 
        "sdG1Count",
        "sdG2Count",
        "sdG1Count2",
        "sdG2Count2",
        "sdG1Count3",
        "sdG2Count3",
        "g1Mean",
        "g2Mean",
        "g1SD",
        "g2SD",
        "numG1",
        "numG2",
        "g1Mean2",
        "g2Mean2",
        "g1SD2",
        "g2SD2",
        "pop2NumG1",
        "pop2NumG2",
        "g1Mean3",
        "g2Mean3",
        "g1SD3",
        "g2SD3",
        "pop3NumG1",
        "pop3NumG2",
        "residual"
    )
    
    for(k in seq_len(length(flowNameDs))){
        flowName <- flowCore::read.FCS(
            paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
        )
        flowData <- .smoothData( flowName, xVariable, 5)
        flowDataMeans <- modelData02%>% dplyr::filter(
            data == flowNameDs[k]
        )
        
        if(0 %in% flowData$x){
            flowData <- flowData %>% dplyr::filter(x != 0)
        }
        
        midpoint1<- flowData[
            which(
                abs(
                    flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2))
            ) == min(
                abs(flowData$x-mean(
                        c(flowDataMeans$G1_1, flowDataMeans$G1_2))
                )
                )
            ),
        ]
        
        if(nrow(midpoint1)>1){
            midpoint1 <- midpoint1[nrow(midpoint1), ]
        }
        
        midpoint2<- flowData[
            which(
                abs(
                    flowData$x-mean(c(midpoint1$x, flowDataMeans$G1_3))
                ) == min(
                abs(flowData$x-mean(c(midpoint1$x, flowDataMeans$G1_3)))
                )
            ),
        ]
        
        if(nrow(midpoint2)>1){
            midpoint2 <- midpoint2[nrow(midpoint2), ]
        }
        
        midpoint3<- flowData[
            which(
                abs(
                    flowData$x-mean(c(midpoint2$x, flowDataMeans$G2_2))
                ) == min(
                abs(flowData$x-mean(c(midpoint2$x, flowDataMeans$G2_2)))
                )
            ),
        ]
        
        if(nrow(midpoint3)>1){
            midpoint3 <- midpoint3[nrow(midpoint3), ]
        }
        
        midpoint4<-flowData[
            which(
                abs(
                    flowData$x-mean(
                        c(
                            max(c(midpoint3$x,flowDataMeans$G2_2)),
                            flowDataMeans$G2_3
                        )
                    )
                ) == min(
                abs(
                    flowData$x-mean(
                    c(
                        max(c(midpoint3$x,flowDataMeans$G2_2)),
                        flowDataMeans$G2_3
                    )
                    )
                )
                )
            ),
        ]
        
        if(nrow(midpoint4)>1){
            midpoint4 <- midpoint4[nrow(midpoint4), ]
        }
        
        #pop 1 G1 standard deviation 
        g1LeftFlowData <- flowData[
            seq_len(which(flowData$x == flowDataMeans$G1_1)-1), ]
        g1RightFlowData <- flowData[
            (which(flowData$x == flowDataMeans$G1_1)+1):
            which(flowData$x == midpoint1$x)-1,
        ]
        
        g1LeftPeak1 <- g1LeftFlowData[
            which(
                abs(
                    g1LeftFlowData$y-flowDataMeans$sdG1Count
                ) == min(abs(g1LeftFlowData$y-flowDataMeans$sdG1Count))
            ), 
        ]
        
        if(nrow(g1LeftPeak1)>1){
            g1LeftPeak1 <- g1LeftPeak1[1, ]
        }
        
        g1RightPeak1 <- g1RightFlowData[
            which(
                abs(
                    g1RightFlowData$y-flowDataMeans$sdG1Count
                ) == min(abs(g1RightFlowData$y-flowDataMeans$sdG1Count))
            ),
        ]
        
        if(nrow(g1RightPeak1)>1 ){
            g1RightPeak1 <- g1RightPeak1[1, ]
        }
        
        #pop 2 G1 standard deviation
        g1LeftFlowData2 <- flowData[
            which(flowData$x == midpoint1$x):
            (which(flowData$x == flowDataMeans$G1_2)-1),
        ]
        
        g1RightFlowData2 <- flowData[
            (which(flowData$x == flowDataMeans$G1_2)+1):
            which(flowData$x == midpoint2$x)-1,
        ]
        
        pop2g1LeftPeak1 <- g1LeftFlowData2[
            which(
                abs(
                    g1LeftFlowData2$y-flowDataMeans$sdG1Count2
                ) == min(abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2))
            ), 
        ]
        
        if(nrow(pop2g1LeftPeak1)>1){
            pop2g1LeftPeak1 <- pop2g1LeftPeak1[1, ]
        }
        
        pop2g1RightPeak1 <- g1RightFlowData2[
            which(
                abs(
                    g1RightFlowData2$y-flowDataMeans$sdG1Count2
                ) == min(abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2))
            ),
        ]
        
        if(nrow(pop2g1RightPeak1)>1 ){
            pop2g1RightPeak1 <- pop2g1RightPeak1[1, ]
        }
        
        #pop 3 G1 standard deviation
        g1LeftFlowData3 <- flowData[
            which(flowData$x == midpoint2$x):
            (which(flowData$x == flowDataMeans$G1_3)-1),
        ]
        
        g1RightFlowData3 <- flowData[
            (which(flowData$x == flowDataMeans$G1_3)+1):
            which(flowData$x == midpoint3$x)-1,
        ]
        
        pop3g1LeftPeak1 <- g1LeftFlowData3[
            which(
                abs(
                    g1LeftFlowData3$y-flowDataMeans$sdG1Count3
                ) == min(abs(g1LeftFlowData3$y-flowDataMeans$sdG1Count3))
            ), 
        ]
        
        if(nrow(pop3g1LeftPeak1)>1){
            pop3g1LeftPeak1 <- pop3g1LeftPeak1[1, ]
        }
        
        pop3g1RightPeak1 <- g1RightFlowData3[
            which(
                abs(
                    g1RightFlowData3$y-flowDataMeans$sdG1Count3
                ) == min(abs(g1RightFlowData3$y-flowDataMeans$sdG1Count3))
            ),
        ]
        
        if(nrow(pop3g1RightPeak1)>1 ){
            pop3g1RightPeak1 <- pop3g1RightPeak1[1, ]
        }
        
        #pop 2 G2 standard deviation
        g2LeftFlowData2 <- flowData[
            which(flowData$x == midpoint3$x):
            (which(flowData$x == flowDataMeans$G2_2)-1),
        ]
        
        g2RightFlowData2 <- flowData[
            (which(flowData$x == flowDataMeans$G2_2)+1):
            which(flowData$x == midpoint4$x)-1,
        ]
        
        pop2g2LeftPeak1 <- g2LeftFlowData2[
            which(
                abs(
                    g2LeftFlowData2$y-flowDataMeans$sdG2Count2
                ) == min(abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2))
            ), 
        ]
        
        if(nrow(pop2g2LeftPeak1)>1){
            pop2g2LeftPeak1 <- pop2g2LeftPeak1[1, ]
        }
        
        pop2g2RightPeak1 <- g2RightFlowData2[
            which(
                abs(
                    g2RightFlowData2$y-flowDataMeans$sdG2Count2
                ) == min(abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2))
            ),
        ]
        
        if(nrow(pop2g2RightPeak1)>1 ){
            pop2g2RightPeak1 <- pop2g2RightPeak1[1, ]
        }
        
        #pop 2 G3 standard deviation
        g2LeftFlowData3 <- flowData[
            which(flowData$x == midpoint4$x):
            (which(flowData$x == flowDataMeans$G2_3)-1),
        ]
        
        g2RightFlowData3 <- flowData[
            (which(flowData$x == flowDataMeans$G2_3)+1):
            nrow(flowData),
        ]
        
        pop3g2LeftPeak1 <- g2LeftFlowData3[
            which(
                abs(
                    g2LeftFlowData3$y-flowDataMeans$sdG2Count3
                ) == min(abs(g2LeftFlowData3$y-flowDataMeans$sdG2Count3))
            ), 
        ]
        
        if(nrow(pop3g2LeftPeak1)>1){
            pop3g2LeftPeak1 <- pop3g2LeftPeak1[1, ]
        }
        
        pop3g2RightPeak1 <- g2RightFlowData3[
            which(
                abs(
                    g2RightFlowData3$y-flowDataMeans$sdG2Count3
                ) == min(abs(g2RightFlowData3$y-flowDataMeans$sdG2Count3))
            ),
        ]
        
        if(nrow(pop3g2RightPeak1)>1 ){
            pop3g2RightPeak1 <- pop3g2RightPeak1[1, ]
        }
        
        flowDataMeans01 = flowDataMeans %>% dplyr::mutate(
            g1Mean=G1_1,
            g2Mean=G2_1,
            g1SD=mean(
                c(
                    g1RightPeak1$x-flowDataMeans$G1_1,
                    flowDataMeans$G1_1-g1LeftPeak1$x
                )
            ),
            g2SD=mean(
                c(
                    pop3g1RightPeak1$x-flowDataMeans$G2_1,
                    flowDataMeans$G2_1-pop3g1LeftPeak1$x
                )
            ),
            pop1NumG1=sum(
                flowData$y[
                    c(
                        which(flowData$x == g1LeftPeak1$x):
                        which(flowData$x == pop3g1RightPeak1$x)
                    )
                ]
            ),
            pop1NumG2=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop3g1LeftPeak1$x):
                        which(flowData$x == pop3g1RightPeak1$x)
                    )
                ]
            ),
            g1Mean2=G1_2,
            g2Mean2=G2_2,
            g1SD2=mean(
                c(
                    pop2g1RightPeak1$x-flowDataMeans$G1_2,
                    flowDataMeans$G1_2-pop2g1LeftPeak1$x
                )
            ),
            g2SD2=mean(
                c(
                pop2g2RightPeak1$x-flowDataMeans$G2_2,
                flowDataMeans$G2_2-pop2g2LeftPeak1$x
                )
            ),
            pop2NumG1=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop2g1LeftPeak1$x):
                        which(flowData$x == pop2g1RightPeak1$x)
                    )
                ]
            ),
            pop2NumG2=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop2g2LeftPeak1$x):
                        which(flowData$x == pop2g2RightPeak1$x)
                    )
                ]
            ),
            
            g1Mean3=G1_3,
            g2Mean3=G2_3,
            g1SD3=mean(
                c(
                pop3g1RightPeak1$x-flowDataMeans$G1_3,
                flowDataMeans$G1_3-pop3g1LeftPeak1$x
                )
            ),
            g2SD3=mean(
                c(
                pop3g2RightPeak1$x-flowDataMeans$G2_3,
                flowDataMeans$G2_3-pop3g2LeftPeak1$x
                )
            ),
            pop3NumG1=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop3g1LeftPeak1$x):
                        which(flowData$x == pop3g2RightPeak1$x)
                    )
                ]
            ),
            pop3NumG2=sum(
                flowData$y[
                    c(
                        which(flowData$x == pop3g2LeftPeak1$x):
                        which(flowData$x == pop3g2RightPeak1$x)
                    )
                ]
            )
        )
        
        g1Mean <- flowDataMeans01$g1Mean
        g2Mean <- flowDataMeans01$g2Mean
        g1SD <- flowDataMeans01$g1SD
        g2SD <- flowDataMeans01$g2SD
        g1Mean2 <- flowDataMeans01$g1Mean2
        g2Mean2 <- flowDataMeans01$g2Mean2
        g1SD2 <- flowDataMeans01$g1SD2
        g2SD2 <- flowDataMeans01$g2SD2
        g1Mean3 <- flowDataMeans01$g1Mean3
        g2Mean3 <- flowDataMeans01$g2Mean3
        g1SD3 <- flowDataMeans01$g1SD3
        g2SD3 <- flowDataMeans01$g2SD3
        
        xpectr::suppress_mw(
            multiPopNLS <- nls(
                formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*
                exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
                (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
                (pop2N1/(sqrt(2*pi)*g1SD2)*exp(-((x-g1Mean2)^2)/(2*g1SD2^2)))+
                (pop2N2/(sqrt(2*pi)*g2SD2)* exp(-((x-g2Mean2)^2)/(2*g2SD2^2)))+
                (pop3N2/(sqrt(2*pi)*g2SD3)* exp(-((x-g2Mean3)^2)/(2*g2SD3^2)))+
                (A + B*x + C*(x^2))*
                (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                (2*(g1SD*(x/g1Mean))^2)))+
                (A2 + B2*x + C2*(x^2))*
                (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                (2*(g2SD2*(x/g2Mean2))^2)))+
                (A3 + B3*x + C3*(x^2))*
                (1/(sqrt(2*pi)*g2SD3*(x/g2Mean3))*exp(-((x-g2Mean3)^2)/
                (2*(g2SD3*(x/g2Mean3))^2)))+
                (A4 + B4*x + C4*(x^2))*
                (1/(sqrt(2*pi)*g2SD*(x/g2Mean))*exp(-((x-g2Mean)^2)/
                (2*(g2SD*(x/g2Mean))^2))),
                data=flowData,
                start=c(
                    pop1N1=flowDataMeans01$pop1NumG1,
                    pop1N2=flowDataMeans01$pop1NumG2,
                    pop2N1=flowDataMeans01$pop2NumG1,
                    pop2N2=flowDataMeans01$pop2NumG2,
                    pop3N2=flowDataMeans01$pop3NumG2,
                    A=0, B=0, C=0,
                    A2=0, B2=0, C2=0,
                    A3=0, B3=0, C3=0,
                    A4=0, B4=0, C4=0
                ),
                nls.control(warnOnly=TRUE)
            )
        )
        
        if(saveGraph == TRUE){
            nlsGraph <- ggplot()+
            geom_line(data=flowData, aes(x=x, y=y, col ='Raw Data'), size=1)+
            geom_line(aes(flowData$x,predict(multiPopNLS), col='NLS'))+
            labs(y="Counts", x=xVariable)+
            theme(legend.title=element_blank())
            
            plotDir <- "nlsGraphs"
            dir.create(
                file.path(dirname(flowDir), plotDir),
                showWarnings=FALSE
            )
            plotInitDir <- "nlsMultiple_graphs"
            dir.create(
                file.path(dirname(flowDir) ,plotDir, plotInitDir),
                showWarnings=FALSE
            )
            plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
            
            png(
                paste0(plotOutFile, "/", flowNameDs[k], '.png'),
                width=600, height=400
            )
          print(nlsGraph)
          dev.off()
        }
        
        flowDataMeans01$residual3Pop <- summary(multiPopNLS)[["sigma"]]
        
        residualMultiDs <- rbind(
            residualMultiDs,
            flowDataMeans01
        )
      }
      
      return(residualMultiDs)
}
