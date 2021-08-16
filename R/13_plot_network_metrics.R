rm(list=ls())
library(ggplot2)

results_dir <- file.path(getwd(),"results")
networks_dir <- file.path(results_dir,"network_construction")
sim_dir <- file.path(results_dir,"simulations")
networkMetrics <- readRDS(file = file.path(networks_dir,"networkMetrics.rds"))



networkMetrics <- networkMetrics[networkMetrics$Connected,]
#networkMetrics <- networkMetrics[which((networkMetrics$MGenes + networkMetrics$Egenes)>4),]

ggplot(networkMetrics) +
  geom_point(aes(x=Nodes,y=Accuracy), alpha=0.5, size = 0.5) +
  theme_bw()+
  xlab("Number of TFs") +
  #  facet_wrap(~Dataset, ncol = 3) +
  theme(text=element_text(size=20),
        axis.text.x=element_text(angle=90, hjust=1))

ggplot(networkMetrics) +
  geom_point(aes(x=Interactions,y=Accuracy),  alpha=0.5, size = 0.5) +
  theme_bw()+
  xlab("Number of Interactions") +
  # facet_wrap(~Dataset, ncol = 3) +
  theme(text=element_text(size=20),
        axis.text.x=element_text(angle=90, hjust=1))

plotData <- networkMetrics
plotData$EMNodes <- (networkMetrics$MGenes + networkMetrics$Egenes)/networkMetrics$Nodes


ggplot(plotData) +
  geom_boxplot(aes(x=EMNodes ,y=Accuracy), alpha=0.5, size = 0.5) +
  theme_bw()+
  xlab("Fraction of EM TFs") +
  #  facet_wrap(~Dataset, ncol = 3) +
  theme(text=element_text(size=20),
        axis.text.x=element_text(angle=90, hjust=1))


plotData <- networkMetrics 
plotData <- plotData[plotData$Connected == TRUE,]
plotData$DiffFraction <- (plotData$Interactions - plotData$PosInt)/plotData$Interactions

ggplot(plotData, aes(x=Interactions,y=Accuracy)) +
  geom_point() +
  xlab("Number of Interactions") +
  theme_bw() +
  theme(text=element_text(size=20)
  )

ggplot(plotData, aes(x=Interactions,y=Accuracy)) +
  geom_point() +
  xlab("Number of Interactions") +
  theme_bw() +
  theme(text=element_text(size=20)
  )
plotData$EMFraction <- (plotData$Egenes+ plotData$MGenes)/plotData$Nodes
ggplot(plotData, aes(x=EMFraction,y=Accuracy)) +
  geom_point() +
  xlab("Fraction of EM TFs") +
  theme_bw() +
  theme(text=element_text(size=20)
  )


ggplot(plotData, aes(x=EMFraction,y=Accuracy)) +
  geom_point() +
  xlab("Fraction of EM TFs") +
  theme_bw() +
  facet_wrap(~Dataset) +
  theme(text=element_text(size=20)
  )

ggplot(plotData, aes(x=EMFraction,y=Accuracy)) +
  geom_point() +
  xlab("Fraction of EM TFs") +
  theme_bw() +
  facet_wrap(~Dataset) +
  theme(text=element_text(size=20)
  )



plotData$posIntFrac <- plotData$PosInt/plotData$Interactions
ggplot(plotData, aes(x=posIntFrac,y=Accuracy)) +
  geom_point() +
  xlab("Fraction of Excitatory Interactions") +
  theme_bw() +
  facet_wrap(~Dataset) +
  theme(text=element_text(size=20)
  )

plotData$NegIntFrac <- 1- plotData$posIntFrac
ggplot(plotData, aes(x=NegIntFrac,y=Accuracy)) +
  geom_point() +
  xlab("Fraction of Inhibitory Interactions") +
  theme_bw() +
  facet_wrap(~Dataset) +
  theme(text=element_text(size=20)
  )
ggplot(plotData, aes(x=EMFraction,y=Accuracy)) +
  geom_point() +
  xlab("Fraction of EM TFs") +
  theme_bw() +
  facet_wrap(~Dataset) +
  theme(text=element_text(size=20)
  )
