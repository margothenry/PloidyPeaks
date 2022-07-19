#' popConfidence3Pop
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
#' popConfidence3Pop(
#'  flowDir = "FlowData/T10_FLC/gated_data",
#'  ds = data,
#'  xVariable = "FITC-A",
#'  saveGraph = TRUE
#'  )
#'
#'
popConfidence3Pop = function(flowDir, ds, xVariable, saveGraph = TRUE){
  
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
    
    midpoint1<- flowData[
      which(
        abs(
          flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2))
        ) == min(
          abs(flowData$x-mean(c(flowDataMeans$G1_1, flowDataMeans$G1_2)))
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
            c(max(c(midpoint3$x,flowDataMeans$G2_2)), flowDataMeans$G2_3))
        ) == min(
          abs(
            flowData$x-mean(
              c(max(c(midpoint3$x,flowDataMeans$G2_2)), flowDataMeans$G2_3)
              )
            )
        )
      ),
    ]
    
    if(nrow(midpoint4)>1){
      midpoint4 <- midpoint4[nrow(midpoint4), ]
    }
    
    #pop 1 G1 standard deviation
    g1LeftFlowData <- flowData[1:which(flowData$x == flowDataMeans$G1_1)-1, ]
    g1RightFlowData <- flowData[
      (
        which(flowData$x == flowDataMeans$G1_1)+1
      ):which(flowData$x == midpoint1$x)-1,
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
      which(flowData$x == midpoint1$x):(which(flowData$x == flowDataMeans$G1_2)-1), ]
    
    g1RightFlowData2 <- flowData[
      (
        which(flowData$x == flowDataMeans$G1_2)+1
      ):which(flowData$x == midpoint2$x)-1,
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
      which(flowData$x == midpoint2$x):(which(flowData$x == flowDataMeans$G1_3)-1), ]
    
    g1RightFlowData3 <- flowData[
      (
        which(flowData$x == flowDataMeans$G1_3)+1
      ):which(flowData$x == midpoint3$x)-1,
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
      which(flowData$x == midpoint3$x):(which(flowData$x == flowDataMeans$G2_2)-1), ]
    
    g2RightFlowData2 <- flowData[
      (
        which(flowData$x == flowDataMeans$G2_2)+1
      ):which(flowData$x == midpoint4$x)-1,
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
      which(flowData$x == midpoint4$x):(which(flowData$x == flowDataMeans$G2_3)-1), ]
    
    g2RightFlowData3 <- flowData[
      (
        which(flowData$x == flowDataMeans$G2_3)+1
      ):nrow(flowData),
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
        c(g1RightPeak1$x-flowDataMeans$G1_1, flowDataMeans$G1_1-g1LeftPeak1$x)
      ),
      g2SD=mean(
        c(pop3g1RightPeak1$x-flowDataMeans$G2_1, flowDataMeans$G2_1-pop3g1LeftPeak1$x)
      ),
      pop1NumG1=sum(
        flowData$y[
          c(which(flowData$x == g1LeftPeak1$x):which(flowData$x == pop3g1RightPeak1$x))
        ]
      ),
      pop1NumG2=sum(
        flowData$y[
          c(which(flowData$x == pop3g1LeftPeak1$x):which(flowData$x == pop3g1RightPeak1$x))
        ]
      ),
      g1Mean2=G1_2,
      g2Mean2=G2_2,
      g1SD2=mean(
        c(pop2g1RightPeak1$x-flowDataMeans$G1_2, flowDataMeans$G1_2-pop2g1LeftPeak1$x)
      ),
      g2SD2=mean(
        c(pop2g2RightPeak1$x-flowDataMeans$G2_2, flowDataMeans$G2_2-pop2g2LeftPeak1$x)
      ),
      pop2NumG1=sum(
        flowData$y[
          c(which(flowData$x == pop2g1LeftPeak1$x):which(flowData$x == pop2g1RightPeak1$x))
        ]
      ),
      pop2NumG2=sum(
        flowData$y[
          c(which(flowData$x == pop2g2LeftPeak1$x):which(flowData$x == pop2g2RightPeak1$x))
        ]
      ),
      
      g1Mean3=G1_3,
      g2Mean3=G2_3,
      g1SD3=mean(
        c(pop3g1RightPeak1$x-flowDataMeans$G1_3, flowDataMeans$G1_3-pop3g1LeftPeak1$x)
      ),
      g2SD3=mean(
        c(pop3g2RightPeak1$x-flowDataMeans$G2_3, flowDataMeans$G2_3-pop3g2LeftPeak1$x)
      ),
      pop3NumG1=sum(
        flowData$y[
          c(which(flowData$x == pop3g1LeftPeak1$x):which(flowData$x == pop3g2RightPeak1$x))
        ]
      ),
      pop3NumG2=sum(
        flowData$y[
          c(which(flowData$x == pop3g2LeftPeak1$x):which(flowData$x == pop3g2RightPeak1$x))
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
        formula = y ~ (pop1N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
          (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
          
          (pop2N1/(sqrt(2*pi)*g1SD2)*exp(-((x-g1Mean2)^2)/(2*g1SD2^2))) +
          
          (pop2N2/(sqrt(2*pi)*g2SD2)* exp(-((x-g2Mean2)^2)/(2*g2SD2^2)))+
          (pop3N2/(sqrt(2*pi)*g2SD3)* exp(-((x-g2Mean3)^2)/(2*g2SD3^2))) +
          
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
                                                (2*(g2SD*(x/g2Mean))^2)))
        ,
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
    
    flowDataMeans01$residual3Pop <- summary(multiPopNLS)[["sigma"]]
    
    residualMultiDs <- rbind(
      residualMultiDs,
      flowDataMeans01
    )
    }
  
  return(residualMultiDs)
}



