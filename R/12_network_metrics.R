# 12 - network metrics: some computations to evaluate the generated networks and simulations
# First, we binarize the simulated expression data for each condition, then count up the expression of E and M genes to identify the number of E and M models
# We also compute connectivity, transitivity, and mean distance of each network using the igraph package

rm(list=ls())
library(Binarize)
library(sRACIPE)

results_dir <- file.path(getwd(),"results")
networks_dir <- file.path(results_dir,"network_construction")
sim_dir <- file.path(results_dir,"simulations")

EMClassCondition <- readRDS(file.path(results_dir,"EMClassCondition.rds"))

fileNames <- list.files(sim_dir)
networkMetrics <- data.frame(Nodes=integer(),Interactions=integer(),PosInt = integer(),
                             EModels = integer(), MModels = integer(),Connected = logical(),
                             MGenes = integer(), Egenes = integer(), ModelCount = integer(), 
                             CellLine = character(), Signal = character(), Dataset = character(), 
                             transitivity = numeric(), MeanDistance = numeric(), Accuracy = numeric(), 
                             Filename = character(), MiPos = numeric(), MiNeg = numeric())

hammingTable <- c(0,7,12,17,22,28,35,42,49,56,64,71,79)
for(fileCounter in seq_along(fileNames)){
  
  fileNameTmp <- fileNames[fileCounter]
  print(fileNameTmp)
  
  cellLine <- str_match(fileNameTmp, "(.*?)_")[,2]
  signal <- str_match(fileNameTmp, "_(.*?)_")[,2]
  posMi <- as.numeric(gsub("^(?:[^_]+_){2}([^_]+).*", "\\1", fileNameTmp))
  negMi <- gsub("^(?:[^_]+_){3}([^_]+).*", "\\1", fileNameTmp)
  negMi <- as.numeric(substr(negMi,0,nchar(negMi)-4))
  
  condition <- paste0(cellLine, "_",signal)
  
  racipe <- readRDS(file.path(sim_dir,fileNameTmp))
  racipe <- sracipeNormalize(racipe)
  expData <- assay(racipe,1)
  
  EMGenes <- EMClassCondition[,condition]
  MState <- EMGenes[!is.na(EMGenes)]
  circuit <-  as.data.frame(sracipeCircuit(racipe))
  g <- graph_from_data_frame(circuit, directed = TRUE, vertices = NULL)
  
  commonGenes <- rownames(expData)
  commonGenes <- commonGenes[which(commonGenes %in% names(MState))]
  expData <- expData[commonGenes,]
  expData <- Binarize::binarizeMatrix(expData, method = "kMeans")[,1:dim(expData)[2]]
  
  MState <- MState[commonGenes]
  EState <- !MState
  expDataTmp <- as.matrix(expData)
  storage.mode(expDataTmp) <- "logical"
  if(length(MState)>85) message("Extend hamming table cutoffs")
  hammingCutoff <- findInterval(length(MState),hammingTable)
  hamDistE <- apply((expDataTmp), 2, function(x) sum(x != EState))
  hamDistM <- apply((expDataTmp), 2, function(x) sum(x != MState))
  networkMetricsTmp <- data.frame(Nodes= length(union(circuit$Source,circuit$Target)),
                                  Interactions=length(circuit$Source),PosInt = length(which(circuit$Type ==1)),
                                  EModels = length(which(hamDistE < hammingCutoff)), MModels = length(which(hamDistM < hammingCutoff)),
                                  Connected = igraph::is_connected(g),MGenes = sum(MState), Egenes = sum(EState), 
                                  ModelCount = sracipeConfig(racipe)$simParams["numModels"], CellLine = cellLine, 
                                  Signal = signal, Dataset = condition, Transitivity = igraph::transitivity(g),
                                  MeanDistance = igraph::mean_distance(g, directed = TRUE, unconnected = FALSE), 
                                  Filename=fileNameTmp, MiPos = posMi, MiNeg = negMi)
  networkMetricsTmp$Accuracy <- sum(networkMetricsTmp$EModels + networkMetricsTmp$MModels)/networkMetricsTmp$ModelCount
  rownames(networkMetricsTmp) <- condition
  networkMetrics <- rbind(networkMetrics, networkMetricsTmp)
}
saveRDS(networkMetrics, file = file.path(networks_dir,"networkMetrics.rds"))
