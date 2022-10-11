#' popConfidence2Pop
#'
#'
#' @param flowDir The directory of the gated .fcs data
#' @param ds dataset to be analyzed
#' @param xVariable The fluorescence channel on the x axis
#' @param saveGraph  T/F for saving the graphs as an output of the NLS
#'
#' @export
#'
#' @examples
#' \dontrun{
#' popConfidence2Pop(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  ds = data,
#'  xVariable = "FITC-A",
#'  saveGraph = TRUE
#'  )
#' }
#'
popConfidence2Pop = function(flowDir, ds, xVariable, saveGraph = TRUE){
  
  ##Removing NOTE 'no visible binding for global variable'
  G1_1<-G2_1<-G1Count_1<-G2Count_1<-x<-y<-NULL
  G1_2<-G2_2<-G1Count_2<-G2Count_2<-NULL
  
  if(TRUE %in% grepl("_3", names(ds))){
    ds1<-ds %>% dplyr::filter(
      is.na(G1_3)
    )    
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
  
  
  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
    )
    flowData <- smoothData( flowName, xVariable, 5)
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
          abs(flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2)))
        )
      ),
    ]

    if(nrow(midPoint)>1){
      midPoint <- midPoint[nrow(midPoint), ]
    }

    #G1 standard deviation
    g1LeftFlowData <- flowData[1:which(flowData$x == flowDataMeans$G1_1)-1, ]
    g1RightFlowData <- flowData[
      (
        which(flowData$x == flowDataMeans$G1_1)+1
      ):which(flowData$x == midPoint$x)-1,
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
          flowData$x-mean(c(flowDataMeans$G2_1, flowDataMeans$G1_2))
        ) == min(
          abs(flowData$x-mean(c(flowDataMeans$G2_1, flowDataMeans$G1_2)))
        )
      ),
    ]

    if(nrow(midPointLeft)>1){
      midPointLeft <- midPointLeft[nrow(midPointLeft), ]
    }

    midPointRight <- flowData[
      which(
        abs(
          flowData$x-mean(c(flowDataMeans$G2_1, flowDataMeans$G2_2))
        ) == min(
          abs(flowData$x-mean(c(flowDataMeans$G2_1, flowDataMeans$G2_2)))
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
            flowDataMeans$G2_1 )
        )
      ) == 3
    ){
      g2LeftFlowData <- flowData[
        which(
          flowData$x == midPoint$x
        ):which(flowData$x == flowDataMeans$G2_1)-1,
      ]
      g2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G2_1)+1):nrow(flowData),
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
    }else{
      g2LeftFlowData <- flowData[
        which(
          flowData$x == midPointLeft$x
        ):which(flowData$x == flowDataMeans$G2_1)-1,
      ]
      g2RightFlowData <- flowData[
        (
          which(
            flowData$x == flowDataMeans$G2_1)+1
        ):which(flowData$x == midPointRight$x),
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

    #pop2

    #G1 standard deviation
    g1LeftFlowData2 <- flowData[
      (
        which(flowData$x == midPoint$x)+1
      ):which(flowData$x == flowDataMeans$G1_2)-1,
    ]

    g1RightFlowData2 <- flowData[
      (
        which(flowData$x == flowDataMeans$G1_2)+1
      ):which(flowData$x == midPointLeft$x)-1,
    ]

    pop2leftPeak1 <- g1LeftFlowData2[
      which(
        abs(
          g1LeftFlowData2$y-flowDataMeans$sdG1Count2
        ) == min(abs(g1LeftFlowData2$y-flowDataMeans$sdG1Count2))
      ),
    ]

    if(nrow(pop2leftPeak1)>1){
      pop2leftPeak1 <- pop2leftPeak1[1, ]
    }

    pop2rightPeak1 <- g1RightFlowData2[
      which(
        abs(
          g1RightFlowData2$y-flowDataMeans$sdG1Count2
        ) == min(abs(g1RightFlowData2$y-flowDataMeans$sdG1Count2))
      ),
    ]

    if(nrow(pop2rightPeak1)>1 ){
      pop2rightPeak1 <- pop2rightPeak1[1, ]
    }

    #G2 standard deviation
    g2LeftFlowData2 <- flowData[
      which(
        flowData$x == midPointRight$x
      ):which(flowData$x == flowDataMeans$G2_2)-1,
    ]
    g2RightFlowData2 <- flowData[
      (which(flowData$x == flowDataMeans$G2_2)+1):nrow(flowData),
    ]

    pop2leftPeak2 <- g2LeftFlowData2[
      which(
        abs(
          g2LeftFlowData2$y-flowDataMeans$sdG2Count2
        ) == min(abs(g2LeftFlowData2$y-flowDataMeans$sdG2Count2)
        )
      ),
    ]

    if(nrow(pop2leftPeak2)>1 ){
      pop2leftPeak2 <- pop2leftPeak2[1, ]
    }

    pop2rightPeak2 <- g2RightFlowData2[
      which(
        abs(
          g2RightFlowData2$y-flowDataMeans$sdG2Count2
        ) == min(abs(g2RightFlowData2$y-flowDataMeans$sdG2Count2)
        )
      ),
    ]

    if(nrow(pop2rightPeak2)>1 ){
      pop2rightPeak2 <- pop2rightPeak2[1, ]
    }

    flowDataMeans01 = flowDataMeans %>% dplyr::mutate(
      g1Mean=G1_1,
      g2Mean=G2_1,
      g1SD=mean(
        c(rightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-leftPeak1$x)
      ),
      g2SD=mean(
        c(rightPeak2$x-flowDataMeans$G2_1, flowDataMeans$G2_1-leftPeak2$x)
      ),
      pop1NumG1=sum(
        flowData$y[
          c(which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x))
        ]
      ),
      pop1NumG2=sum(
        flowData$y[
          c(which(flowData$x == leftPeak2$x):which(flowData$x == rightPeak2$x))
        ]
      ),
      g1Mean2=G1_2,
      g2Mean2=G2_2,
      g1SD2=mean(
        c(pop2rightPeak1$x-flowDataMeans$G1_2, flowDataMeans$G1_2-pop2leftPeak1$x)
      ),
      g2SD2=mean(
        c(pop2rightPeak2$x-flowDataMeans$G2_2, flowDataMeans$G2_2-pop2leftPeak2$x)
      ),
      pop2NumG1=sum(
        flowData$y[
          c(which(flowData$x == pop2leftPeak1$x):which(flowData$x == pop2rightPeak1$x))
        ]
      ),
      pop2NumG2=sum(
        flowData$y[
          c(which(flowData$x == pop2leftPeak2$x):which(flowData$x == pop2rightPeak2$x))
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
            flowDataMeans$G2_1 )
        )
      ) == 3
    ){
      xpectr::suppress_mw(
        multiPopNLS <- nls(
          formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
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
          formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
            (A + B*x + C*(x^2))*
            (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
                                                  (2*(g1SD*(x/g1Mean))^2)))+
            (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
            (A2 + B2*x + C2*(x^2))*
            (1/(sqrt(2*pi)*g2SD2*(x/g1Mean2))*exp(-((x-g2Mean2)^2)/
                                                    (2*(g2SD2*(x/g2Mean2))^2)))+
            (pop2N1/(sqrt(2*pi)*g1SD2)*exp(-((x-g1Mean2)^2)/(2*g1SD2^2))) +
            (pop2N2/(sqrt(2*pi)*g2SD2)* exp(-((x-g2Mean2)^2)/(2*g2SD2^2))),
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
      plotInitDir <- "NLSMultiple_graphs"
      dir.create(
        file.path(dirname(flowDir) ,plotDir, plotInitDir), showWarnings=FALSE
      )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)
      
      png(paste0(plotOutFile, "/", flowNameDs[k], '.png'), width=600, height=400)
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



