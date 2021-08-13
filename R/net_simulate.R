rm(list = ls())
library(sRACIPE)
library(stringr)

working_dir <- "output/network_construction/"
simulation_dir <-"simAll/"

dir.create(file.path("simAll/"))

createNetworkMiValue <- function(tmpNetwork=tmpNetwork,posMi = posMi,
                                 negMi = negMi, name = "A549", simulate = FALSE){
  tmpNetwork <- tmpNetwork[order(tmpNetwork$Mi, decreasing = TRUE),]
  tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork[c("Source","Target","Type")]),]
  posNetwork <- tmpNetwork[tmpNetwork$Cor >0,]
  posNetwork <- posNetwork[which(posNetwork$Mi > posMi),]
  negNetwork <- tmpNetwork[tmpNetwork$Cor < 0,]
  negNetwork <- negNetwork[which(negNetwork$Mi > negMi),]
  if(dim(posNetwork)[1] == 0 | dim(negNetwork)[1] == 0) return()
  
  tmpNetwork <- rbind(posNetwork,negNetwork)
  sum(duplicated(tmpNetwork[,c("Source","Target")]))
  sum(duplicated(tmpNetwork[,c("Source","Target","Type")]))
  tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork[,c("Source","Target","Type")]),]
  tmpNetwork[duplicated(tmpNetwork[,c("Source","Target")]),]
  tmpNetwork[duplicated(tmpNetwork[,c("Source","Target")], fromLast = TRUE),]
  delRows <- anyDuplicated(tmpNetwork[,c("Source","Target")])
  delRows <- c(delRows, anyDuplicated(tmpNetwork[,c("Source","Target")], fromLast = TRUE))
  delRows <- delRows[delRows>0]
  if(length(delRows>0)){ tmpNetwork <- tmpNetwork[-delRows,]}
  
  interactionRemoved = length(tmpNetwork$Source)
  # If signal nodes are to be removed iteratively, uncomment the next line
  # while(interactionRemoved>0)
  {
    tmpVar = length(tmpNetwork$Source)
    
    tmpNetwork <- tmpNetwork[(which(tmpNetwork$Source %in% tmpNetwork$Target )),]
    #     tmpNetwork <- tmpNetwork[which(tmpNetwork$Target %in% tmpNetwork$Source),]
    interactionRemoved = tmpVar - length(tmpNetwork$Source)
  }
  
  tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork),]
  message(length(tmpNetwork$Source))
  networkTfs <- union(tmpNetwork$Source,tmpNetwork$Target)
  coreTfGenes <- substr(coreTf,0,nchar(coreTf)-3)
  print(networkTfs[which(networkTfs %in% coreTfGenes)])
  
  if(length(union(tmpNetwork$Source,tmpNetwork$Target)) <5) simulate = FALSE
  #  if(length(which(tmpNetwork$Type ==2)) <2) simulate = FALSE
  #  if(length(which(tmpNetwork$Type ==1)) <2) simulate = FALSE
  
  if(simulate){
    tmp <- sracipeSimulate(circuit = tmpNetwork[,c("Source", "Target", "Type")], plotToFile = TRUE,
                           numModels = 2000, plots = FALSE, stepper = "RK4", integrateStepSize = 0.02)
    annotation(tmp) <- paste0(name,"_",as.character(posMi), "_",as.character(negMi),"Interactions")
    # tmp <- sracipePlotData(tmp, plotToFile = T)
    saveRDS(tmp, file = paste0(simulation_dir,name,"_",as.character(posMi), "_",as.character(negMi),".rds"))
  }
  # return(tmpNetwork)
}

coreTf <- readRDS(paste0(working_dir,"coreTf.rds"))
allNetworkFull <- readRDS(paste0(working_dir,"allNetworkFullCoreMiTop23AllRegAllInt.rds"))

for(posMi in seq(1,0.05,-0.05)){
  for(negMi in seq(0.5,0.05,-0.05)){
    cellLine = "DU145"
    treatment <- "TNF"
    print(paste0(cellLine,"_",treatment,"_",as.character(posMi),"_",as.character(negMi)))
    tmpNetwork <- allNetworkFull[(allNetworkFull$CellLine == cellLine & allNetworkFull$Signal == treatment),]
    
    createNetworkMiValue(tmpNetwork,posMi=posMi,negMi = negMi,
                         name=paste0(cellLine,"_",treatment), simulate = TRUE)
  }
}


#### check results ##### SAMPLE  ##########
tmp <- readRDS(paste0(simulation_dir,"DU145_TNF_0.6_0.05.rds"))
tmp <- sracipePlotData(tmp)

