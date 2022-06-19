#' popConfidenceInitial
#'
#'
#' @param flowDir The directory of the gated .fcs data
#' @param ds dataset to be analyzed
#' @param x_variable The fluorescence channel on the x axis
#' @param saveGraph  T/F for saving the graphs as an output of the NLS
#'
#' @export
#'
#' @examples
#' popConfidenceInitial(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  ds = data,
#'  xVariable = "FITC-A",
#'  saveGraph = TRUE
#'  )
#'
#'
popConfidenceInitial = function(flowDir, ds, xVariable, saveGraph = TRUE){

  modelData <- ds %>% dplyr::select(
    data,
    G1_1,
    G2_1,
    G1Count_1,
    G2Count_1
  )

  modelData01 <- modelData %>% dplyr::mutate(
    sdG1Count=G1Count_1*0.6,
    sdG2Count=G2Count_1*0.6
  )

  flowNameDs <- unique(modelData01$data)

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

  for(k in 1:length(flowNameDs)){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowNameDs[k]), transformation=FALSE
      )
    flowData <- smoothData( flowName, xVariable, 5)
    flowDataMeans <- modelData01%>% dplyr::filter(
      data == flowNameDs[k]
    )

    if(0 %in% flowData$x){
      flowData <- flowData %>% dplyr::filter(x != 0)
    }

    if(flowDataMeans$G2_1 != 0){
      midPoint <- flowData[
        which(
          abs(
            flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G2_1))
            ) == min(
              abs(flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G2_1)))
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

      #G2 standard deviation
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
          c(rightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-leftPeak1$x)
          ),
        g2SD=mean(
          c(rightPeak2$x-flowDataMeans$G2_1, flowDataMeans$G2_1-leftPeak2$x)
          ),
        numG1=sum(
          flowData$y[
            c(which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x))
            ]
          ),
        numG2=sum(
          flowData$y[
            c(which(flowData$x == leftPeak2$x):which(flowData$x == rightPeak2$x))
            ]
          )
      )
      g1Mean <- flowDataMeans01$g1Mean
      g2Mean <- flowDataMeans01$g2Mean
      g1SD <- flowDataMeans01$g1SD
      g2SD <- flowDataMeans01$g2SD

      xpectr::suppress_mw(
      singlePopNLS <- nls(
          formula = y ~ (N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
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
        1:which(flowData$x == flowDataMeans$G1_1)-1,
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
          c(rightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-leftPeak1$x)
          ),
        g2SD=NA,
        numG1=sum(
          flowData$y[
            c(which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x))
            ]
          ),
        numG2=0
      )

      g1Mean <- flowDataMeans01$g1Mean
      g1SD <- flowDataMeans01$g1SD
      
      xpectr::suppress_mw(
      singlePopNLS <- nls(
        formula = y ~ (N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
          (A + B*x + C*(x^2))*
          (1/(sqrt(2 * pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
            (2*(g1SD*(x/g1Mean))^2))),
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
      plotInitDir <- "NLSInitial_graphs"
      dir.create(
        file.path(dirname(flowDir) ,plotDir, plotInitDir), showWarnings=FALSE
        )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)

      png(paste0(plotOutFile, "/", flowNameDs[k], '.png'), width=600, height=400)
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



