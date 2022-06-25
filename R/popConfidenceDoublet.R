#' popConfidenceDoublet
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
#' popConfidenceDoublet(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  ds = data,
#'  xVariable = "FITC-A",
#'  saveGraph = TRUE
#'  )
#'
#'
popConfidenceDoublet = function(flowDir, ds, xVariable, saveGraph = TRUE){

  ds2 <- ds %>% dplyr::filter(doublet == 1)
  flowNameDs <- unique(ds2$data)

  modelData <- ds2 %>% dplyr::select(
    data,
    G1_1,
    G2_1,
    `doublet G1+G2`,
    `doublet G2+G2`,
    G1Count_1,
    G2Count_1,
    `doublet G1+G2 count`,
    `doublet G2+G2 count`
  )

  modelData01 <- modelData %>% dplyr::mutate(
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
      G1LeftFlowData <- flowData[1:which(flowData$x == flowDataMeans$G1_1)-1, ]
      G1RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G1_1)+1):which(flowData$x == midPoint$x),
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

      #G2 standard deviation
      G2LeftFlowData <- flowData[
        which(
          flowData$x == midPoint$x
          ):(which(flowData$x == flowDataMeans$G2_1)-1),
        ]
      G2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$G2_1)+1):nrow(flowData),
        ]

      leftPeak2 <- G2LeftFlowData[
        which(
          abs(
            G2LeftFlowData$y-flowDataMeans$sdG2Count
            ) == min(abs(G2LeftFlowData$y-flowDataMeans$sdG2Count))
        ),
      ]

      if(nrow(leftPeak2)>1 ){
        leftPeak2 <- leftPeak2[1, ]
      }

      rightPeak2 <- G2RightFlowData[
        which(
          abs(
            G2RightFlowData$y-flowDataMeans$sdG2Count
            ) == min(abs(G2RightFlowData$y-flowDataMeans$sdG2Count))
        ),
      ]

      if(nrow(rightPeak2)>1 ){
        rightPeak2 <- rightPeak2[1, ]
      }

      #G1+G2 doublet standard deviation
      midPointDoublet <- flowData[
        which(
          abs(
            flowData$x-mean(
              c(flowDataMeans$`doublet G1+G2`, flowDataMeans$G2_1)
              )
            ) == min(
              abs(
                flowData$x-mean(
                  c(flowDataMeans$`doublet G1+G2`, flowDataMeans$G2_1)
                  )
                )
              )
        ),
      ]

      if(nrow(midPointDoublet)>1){
        midPointDoublet <- midPointDoublet[nrow(midPointDoublet), ]
      }

      doubletG1G2LeftFlowData <- flowData[
        which(
          flowData$x == midPointDoublet$x
          ):(which(flowData$x == flowDataMeans$`doublet G1+G2`)-1),
        ]
      doubletG1G2RightFlowData <- flowData[
        (which(flowData$x == flowDataMeans$`doublet G1+G2`)+1):nrow(flowData),
        ]

      doubletG1G2Left <- doubletG1G2LeftFlowData[
        which(
          abs(
            doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count
            ) == min(abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count))
        ),
      ]

      if(nrow(doubletG1G2Left)>1){
        doubletG1G2Left <- doubletG1G2Left[1, ]
      }

      doubletG1G2Right <- doubletG1G2RightFlowData[
        which(
          abs(
            doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count
            ) == min(abs(doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count))
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
                c(flowDataMeans$`doublet G1+G2`, flowDataMeans$`doublet G2+G2`)
                )
              ) == min(
                abs(
                  flowData$x-mean(
                    c(flowDataMeans$`doublet G1+G2`, flowDataMeans$`doublet G2+G2`)
                    )
                  )
                )
          ),
        ]

        if(nrow(midPointDoublet2)>1){
          midPointDoublet2 <- midPointDoublet2[nrow(midPointDoublet2), ]
        }

        doubletG1G2LeftFlowData <- flowData[
          which(
            flowData$x == midPointDoublet2$x
            ):(which(flowData$x == flowDataMeans$`doublet G2+G2`)-1),
          ]
        
        doubletG2G2RightFlowData <- flowData[
          (which(flowData$x == flowDataMeans$`doublet G2+G2`)+1):nrow(flowData),
          ]

        doubletG2G2Left <- doubletG1G2LeftFlowData[
          which(
            abs(
              doubletG1G2LeftFlowData$y-flowDataMeans$sdG2G2Count
              ) == min(abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG2G2Count))
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
                abs(doubletG2G2RightFlowData$y-flowDataMeans$sdG2G2Count)
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
            c(rightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-leftPeak1$x)
            ),
          g2SD=mean(
            c(rightPeak2$x-flowDataMeans$G2_1, flowDataMeans$G2_1-leftPeak2$x)
            ),
          G1g2SD=mean(
            c(doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2` - doubletG1G2Left$x)
            ),
          G2g2SD=mean(
            c(doubletG2G2Right$x-flowDataMeans$`doublet G2+G2`,
              flowDataMeans$`doublet G2+G2` - doubletG2G2Left$x)
            ),
          numG1=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x)
                )
              ]
            ),
          numG2=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak2$x):which(flowData$x == rightPeak2$x)
                )
              ]
            ),
          numDoublet1=sum(
            flowData$y[
              c(
                which(
                  flowData$x == doubletG1G2Left$x
                  ):which(flowData$x == doubletG1G2Right$x)
                )
              ]
            ),
          numDoublet2=sum(
            flowData$y[
              c(
                which(
                  flowData$x == doubletG2G2Left$x
                  ):which(flowData$x == doubletG2G2Right$x)
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
          formula = y ~ (
            N1/(sqrt(2*pi) * g1SD) * exp(-((x-g1Mean)^2)/(2 *g1SD^2))) +
            (A1 + B1*x + C1*(x^2))*
            (1/(sqrt(2 * pi) * g1SD * (x/g1Mean)) * exp(-((x - g1Mean)^2)/
              (2 *(g1SD* (x/g1Mean))^2)))+
            (N2/(sqrt(2 * pi) * g2SD) * exp(-((x - g2Mean)^2)/(2 *g2SD^2)))+
            (A2 + B2*x + C2*(x^2))*
            (1/(sqrt(2 * pi) * g2SD * (x/g2Mean)) * exp(-((x - g2Mean)^2)/
              (2 *(g2SD* (x/g2Mean))^2)))+
            (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
              (2*G1g2SD^2)))+
            (A3 + B3*x + C3*(x^2))*
            (1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*exp(-((x-G1G2Mean)^2)/
              (2*(G1g2SD*(x/G1G2Mean))^2)))+
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
              flowDataMeans$G1_1-leftPeak1$x)
            ),
          g2SD=mean(
            c(
              rightPeak2$x-flowDataMeans$G2_1,
              flowDataMeans$G2_1-leftPeak2$x)
            ),
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x)
            ),
          G2g2SD=NA,
          numG1=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x)
                )
              ]
            ),
          numG2=sum(
            flowData$y[
              c(
                which(flowData$x == leftPeak2$x):which(flowData$x == rightPeak2$x)
                )
              ]
            ),
          numDoublet1=sum(
            flowData$y[
              c(
                which(
                  flowData$x == doubletG1G2Left$x
                  ):which(flowData$x == doubletG1G2Right$x)
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
          formula=y ~ (N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
            (A1 + B1*x + C1*(x^2))*
            (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/
              (2*(g1SD*(x/g1Mean))^2)))+
            (N2/(sqrt(2*pi)*g2SD)*exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
            (A2 + B2*x + C2*(x^2))*
            (1/(sqrt(2*pi)*g2SD*(x/g2Mean))*exp(-((x-g2Mean)^2)/
              (2*(g2SD*(x/g2Mean))^2)))+
            (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
              (2*G1g2SD^2))),
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
      G1LeftFlowData <- flowData[1:which(flowData$x == flowDataMeans$G1_1), ]
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
            flowData$x-mean(c(flowDataMeans$`doublet G1+G2`, flowDataMeans$G1_1))
            ) == min(
              abs(
                flowData$x-mean(
                  c(flowDataMeans$`doublet G1+G2`, flowDataMeans$G1_1))
                )
            )
        ),
      ]

      if(nrow(midPointDoublet)>1){
        midPointDoublet <- midPointDoublet[nrow(midPointDoublet), ]
      }

      doubletG1G2LeftFlowData <- flowData[
        which(
          flowData$x == midPointDoublet$x
          ):which(flowData$x == flowDataMeans$`doublet G1+G2`),
        ]
      
      doubletG1G2RightFlowData <- flowData[
        which(flowData$x == flowDataMeans$`doublet G1+G2`):nrow(flowData),
        ]

      doubletG1G2Left <- doubletG1G2LeftFlowData[
        which(
          abs(
            doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count
            ) == min(abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG1G2Count))
        ),
      ]

      if(nrow(doubletG1G2Left)>1){
        doubletG1G2Left <- doubletG1G2Left[1, ]
      }

      doubletG1G2Right <- doubletG1G2RightFlowData[
        which(
          abs(
            doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count
            ) == min(abs(doubletG1G2RightFlowData$y-flowDataMeans$sdG1G2Count))
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
                c(flowDataMeans$`doublet G1+G2`, flowDataMeans$`doublet G2+G2`)
                )
              ) == min(
                abs(
                  flowData$x-mean(
                    c(flowDataMeans$`doublet G1+G2`, flowDataMeans$`doublet G2+G2`)
                    )
                  )
                )
          ),
        ]

        if(nrow(midPointDoublet2)>1){
          midPointDoublet2 <- midPointDoublet2[nrow(midPointDoublet2), ]
        }

        doubletG1G2LeftFlowData <- flowData[
          which(
            flowData$x == midPointDoublet2$x
            ):which(flowData$x == flowDataMeans$`doublet G2+G2`), 
          ]
        doubletG2G2RightFlowData <- flowData[
          which(flowData$x == flowDataMeans$`doublet G2+G2`):nrow(flowData),
          ]

        doubletG2G2Left <- doubletG1G2LeftFlowData[
          which(
            abs(
              doubletG1G2LeftFlowData$y-flowDataMeans$sdG2G2Count
              ) == min(abs(doubletG1G2LeftFlowData$y-flowDataMeans$sdG2G2Count))
          ),
        ]

        if(nrow(doubletG2G2Left)>1){
          doubletG2G2Left <- doubletG2G2Left[1, ]
        }

        doubletG2G2Right <- doubletG2G2RightFlowData[
          which(
            abs(
              doubletG2G2RightFlowData$y-flowDataMeans$sdG2G2Count
              ) == min(abs(doubletG2G2RightFlowData$y-flowDataMeans$sdG2G2Count))
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
            c(rightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-leftPeak1$x)
            ),
          g2SD=NA,
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x)
            ),
          G2g2SD=mean(
            c(
              doubletG2G2Right$x-flowDataMeans$`doublet G2+G2`,
              flowDataMeans$`doublet G2+G2`-doubletG2G2Left$x)
            ),
          numG1=sum(
            flowData$y[
              c(which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x))
              ]
            ),
          numG2=0,
          numDoublet1=sum(
            flowData$y[
              c(
                which(
                  flowData$x == doubletG1G2Left$x
                  ):which(flowData$x == doubletG1G2Right$x)
                )
              ]
            ),
          numDoublet2=sum(
            flowData$y[
              c(
                which(
                  flowData$x == doubletG2G2Left$x
                  ):which(flowData$x == doubletG2G2Right$x)
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
          formula = y ~ (N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
            (A1 + B1*x + C1*(x^2))*
            (1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x - g1Mean)^2)/
              (2*(g1SD*(x/g1Mean))^2)))+
            (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/
              (2*G1g2SD^2)))+
            (A2 + B2*x + C2*(x^2))*
            (1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*exp(-((x-G1G2Mean)^2)/
              (2*(G1g2SD*(x/G1G2Mean))^2)))+
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
            c(rightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-leftPeak1$x)
            ),
          g2SD=NA,
          G1g2SD=mean(
            c(
              doubletG1G2Right$x-flowDataMeans$`doublet G1+G2`,
              flowDataMeans$`doublet G1+G2`-doubletG1G2Left$x)
            ),
          G2g2SD=NA,
          numG1=sum(
            flowData$y[
              c(which(flowData$x == leftPeak1$x):which(flowData$x == rightPeak1$x))
              ]
            ),
          numG2=0,
          numDoublet1=sum(
            flowData$y[
              c(
                which(
                  flowData$x == doubletG1G2Left$x
                  ):which(flowData$x == doubletG1G2Right$x)
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
          formula = y ~ (N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2))) +
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
        ))

      }

    }

    if(saveGraph == TRUE){
      nlsGraph <- ggplot()+
        geom_line(data=flowData, aes(x=x, y=y, col = 'Raw Data'), size=1)+
        geom_line(aes(flowData$x, predict(singlePopNLS), col ='NLS'))+
        labs(y ="Counts", x =xVariable)+
        theme(legend.title = element_blank())

      plotDir <- "nlsGraphs"
      plotInitDir <- "NLSDoublet_graphs"
      dir.create(
        file.path(dirname(flowDir), plotDir, plotInitDir), showWarnings=FALSE
        )
      plotOutFile <- file.path(dirname(flowDir), plotDir, plotInitDir)

      png(
        paste0(plotOutFile, "/", flowNameDs[k], '.png'), width=600, height=400
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



