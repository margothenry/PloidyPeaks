#' outputData
#'
#' outputData creates a .csv and saves it in a folder called ‘analysis’ for the user.
#' The dataset includes the location of all peaks identified, the height of the peaks,
#' an option to have the information on doublets, a flag if the data was messy,
#' a confidence interval for the population.
#'
#' @param flowDir The data set to analyze for doublets
#' @param singleDs The data set for the single populations found in PeakAlgorithm1
#' @param finishedDs The data set for all other flow frames that have been analyzed
#' @param messyDs The data set that got identified as messy
#' @param xVariable The fluorescence channel on the x axis
#' @param doubletFlag TRUE/FALSE
#' @param saveGraph TRUE/FALSE
#' @export
#'
#' @examples
#' outputData(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  singleDs = single_data,
#'  finishedDs = finished_data,
#'  messyDs = messy_data,
#'  xVariable = "FITC-A",
#'  doubletFlag = FALSE,
#'  saveGraph = TRUE
#'  )

outputData = function(
  flowDir,
  singleDs,
  finishedDs,
  messyDs,
  xVariable,
  doubletFlag,
  saveGraph
  ){
  ##Removing NOTE 'no visible binding for global variable'
  x<-y<-.<-possiblePairX<-possiblePairY<-G1<-G1Count<-G2<-G2Count<-id<-NULL
  g1G2Doublet<-g1G2DoubletCount<-g2G2Doublet<-g2G2DoubletCount<-NULL
  residual<-residualDoublet<-Success<-Algorithm<-Data<-`doublet G1+G2`<-NULL
  `doublet G1+G2 count`<-`doublet G2+G2`<- `doublet G2+G2 count`<-NULL
  
  ##If the algorithm classified each sub population in the first algorithm
  if(
    purrr::is_empty(finishedDs) &
    purrr::is_empty(messyDs) &
    !purrr::is_empty(singleDs)
    ){

    ##Formatting the single pop data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
      dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
      dplyr::mutate(
        messy=0,
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
        "messy"
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
      tidyr::pivot_wider(names_from=id, values_from=c(G1, G1Count, G2, G2Count))

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
    initialRSE=popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
      )

    doubletRSE=popConfidenceDoublet(
      flowDir, ds=finalData3, xVariable, saveGraph
      )

    finalData4=sqldf::sqldf(
      "select ds.*,
              ds2.residual,
              ds3.residualDoublet
       from finalData3 ds
       left join initialRSE ds2
        on ds.data = ds2.data
      left join doubletRSE ds3
        on ds.data = ds3.data"
    )

    finalData5=finalData4 %>% dplyr::mutate(
      finalResidual=do.call(
        pmin, 
        c(subset(., select = c(residual, residualDoublet)), na.rm=TRUE))
    ) %>% dplyr::rename(
      initialRSE=residual
    ) %>% dplyr::select(
      -residualDoublet
    )
    

    ##adding which algorithm analyzed the flow frame
    logDs2=logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3=logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4=logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5=logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)

    ##merging the algorithm used with the final dataset
    finalData6=merge(finalData5, logDs5, by="data")

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
      paste0(file.path(dirname(getwd()), subDir), "/flow_analysis.csv")
      )

  }else if(
    purrr::is_empty(messyDs) &
    !purrr::is_empty(finishedDs) &
    !purrr::is_empty(singleDs)
    ){
    ##If the algorithm classified each sub population before prior the 4th
    ##algorithm

    ##Formatting the diploid data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
      dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
      dplyr::mutate(
        messy=0,
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
        messy=0,
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
        "messy"
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
      tidyr::pivot_wider(names_from=id, values_from=c(G1, G1Count, G2, G2Count))

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
    initialRSE=popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
      )

    doubletRSE=popConfidenceDoublet(
      flowDir, ds=finalData3, xVariable, saveGraph
      )

    finalData4=sqldf::sqldf(
      "select ds.*,
              ds2.residual,
              ds3.residualDoublet
       from finalData3 ds
       left join initialRSE ds2
        on ds.data = ds2.data
      left join doubletRSE ds3
        on ds.data = ds3.data"
    )

    finalData5=finalData4 %>% dplyr::mutate(
      finalResidual=do.call(
        pmin, 
        c(subset(., select=c(residual, residualDoublet)), na.rm=TRUE))
    ) %>% dplyr::rename(
      initialRSE=residual
    ) %>% dplyr::select(
      -residualDoublet
    )
    

    ##adding which algorithm analyzed the flow frame
    logDs2=logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))
    
    logDs3=logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))
    
    logDs4=logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )

    logDs5=logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)

    ##merging the algorithm used with the final dataset
    finalData6=merge(finalData5, logDs5, by="data")

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
      paste0(file.path(dirname(getwd()), subDir), "/flow_analysis.csv"
             )
      )

  }else if(
    purrr::is_empty(singleDs) &
    !purrr::is_empty(finishedDs) &
    !purrr::is_empty(messyDs)
    ){
    ##If the algorithm classified each sub population before prior the 4th
    ##algorithm

    ##Formatting the diploid data from the first peak algorithm
    finalPart1=messyDs %>% data.frame() %>%
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
        messy=0,
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
        "messy"
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
      tidyr::pivot_wider(names_from=id, values_from=c(G1, G1Count, G2, G2Count))
    
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
    initialRSE=popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    doubletRSE=popConfidenceDoublet(
      flowDir, ds=finalData3, xVariable, saveGraph
    )
    
    finalData4=sqldf::sqldf(
      "select ds.*,
              ds2.residual,
              ds3.residualDoublet
       from finalData3 ds
       left join initialRSE ds2
        on ds.data = ds2.data
      left join doubletRSE ds3
        on ds.data = ds3.data"
    )

    finalData5=finalData4 %>% dplyr::mutate(
      finalResidual=do.call(
        pmin, 
        c(subset(., select=c(residual, residualDoublet)), na.rm=TRUE))
    ) %>% dplyr::rename(
      initialRSE=residual
    ) %>% dplyr::select(
      -residualDoublet
    )
    

    ##adding which algorithm analyzed the flow frame
    logDs2=logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))

    logDs3=logDs2 %>%
      dplyr::mutate(id=rowid(Data)) %>%
      tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))

    logDs4=logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )

    logDs5=logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)

    ##merging the algorithm used with the final dataset
    finalData6=merge(finalData5, logDs5, by="data")
    
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
      paste0(file.path(dirname(getwd()), subDir), "/flow_analysis.csv"
      )
    )
    
  }else{

    ##Formatting the diploid data from the first peak algorithm
    finalPart1=singleDs %>% data.frame() %>%
      dplyr::select(-c("g3LL", "g3UL", "g4LL", "g4UL")) %>%
      dplyr::mutate(
        messy=0,
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
        messy=0,
        G1=x,
        G1Count=y,
        G2=possiblePairX,
        G2Count=possiblePairY
      )

    ##Formatting the data from the last peak algorthm
    finalPart3=messyDs %>% data.frame() %>%
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
        "messy"
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
      tidyr::pivot_wider(names_from=id, values_from=c(G1, G2,G1Count, G2Count))

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
    initialRSE=popConfidenceInitial(
      flowDir, ds=finalData3, xVariable, saveGraph
      )

    doubletRSE=popConfidenceDoublet(
      flowDir, ds=finalData3, xVariable, saveGraph
      )

    finalData4=sqldf::sqldf(
      "select ds.*,
              ds2.residual,
              ds3.residualDoublet
       from finalData3 ds
       left join initialRSE ds2
        on ds.data = ds2.data
      left join doubletRSE ds3
        on ds.data = ds3.data"
    )
    
    finalData5=finalData4 %>% dplyr::mutate(
      finalResidual=do.call(
        pmin, 
        c(subset(., select=c(residual, residualDoublet)), na.rm=TRUE))
    ) %>% dplyr::rename(
      initialRSE=residual
    ) %>% dplyr::select(
      -residualDoublet
    )
    

    ##adding which algorithm analyzed the flow frame
    logDs2=logDs %>% 
      dplyr::select(-Success) %>% 
      dplyr::mutate(Algorithm=as.numeric(Algorithm))

    logDs3=logDs2 %>%
    dplyr::mutate(id=rowid(Data)) %>%
    tidyr::pivot_wider(names_from=id, values_from=c(Algorithm))

    logDs4=logDs3 %>% dplyr::mutate(
      Algorithm=do.call(
        pmax, 
        c(subset(., select=2:ncol(logDs3)), na.rm=TRUE))
    )
    
    logDs5=logDs4 %>%
      dplyr::select(Data, Algorithm) %>%
      dplyr::rename(data=Data)

    ##merging the algorithm used with the final dataset
    finalData6=merge(finalData5, logDs5, by="data")

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
      paste0(file.path(dirname(getwd()), subDir), "/flow_analysis.csv")
      )

  }

}



