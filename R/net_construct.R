rm(list = ls())
treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
directions <- c("fwd","rev")

working_dir <- "output/network_construction/"

allRegulons <- list()
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    for(dCounter in seq_along(directions)){
      
      ## Select a treatment to work with (or iterate thru t in a loop)
      t <- treatments[tCounter]
      c <- cancer_types[cCounter]
      direc <- directions[dCounter]
      treatmentName <- paste0(c,"_",t,"_",direc)
      output_dir <- paste0("output/auc/")
      
      regulon_fname <- paste0("EmtRds/SCENIC/results/regulons_",c,"_",t,"_",direc,".json")
      regs <- read_json(regulon_fname, simplify = T)
      network <- data.frame(Source = character(), Target = character())
      for(src in seq_along(regs)){
        targets <- unlist(regs[src])
        networkTmp <- data.frame(Source = rep(names(regs)[src],length(targets)), Target = targets)
        network <- rbind(network, networkTmp)
      }
      network$Dataset <- treatmentName
      allRegulons[[treatmentName]] <- network
      
    }
  }
}

saveRDS(allRegulons, file = "allRegulons.rds")

library(infotheo)
#' @param data expression matrix. Genes in rows and samples/cells in columns
#' @param nbins number of bins
#' @value A mutual information matrix
calculateMI <- function(data = actMat){
#  data <- t(data)
  nGenes <- dim(data)[2]
  miMat <- matrix(0,nrow = nGenes,ncol = nGenes)
  geneNames <- colnames(data)
  rownames(miMat) <- geneNames
  colnames(miMat) <- geneNames
  data <- infotheo::discretize(data)
  #i=1
  for(i in 1:(nGenes-1))
  {
    for(j in (i+1):(nGenes))
    {
      gene1 <- geneNames[i]
      gene2 <- geneNames[j]
      miMat[gene1, gene2] <- infotheo::mutinformation(data[,gene1], data[,gene2], method = "mm")
      miMat[gene2, gene1] <- miMat[gene1, gene2]
    }
  }
  return(miMat)
}

allRegulons <- readRDS(paste0(working_dir,"allRegulons.rds"))
coreTf <- readRDS(paste0(working_dir,"coreTf.rds"))
fullNetworkList <- list()
counter <- 1
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    for(dCounter in seq_along(directions)){
      t <- treatments[tCounter]
      c <- cancer_types[cCounter]
      direc <- directions[dCounter]
      treatmentName <- paste0(c,"_",t,"_",direc)
      print(treatmentName)
      tfInNetwork <- coreTf
      tfInNetwork <- substr(tfInNetwork, 0, nchar(as.character(tfInNetwork)) -3)
      geneNetwork <- data.frame(allRegulons[treatmentName])
      colnames(geneNetwork) <- c("Source", "Target", "Dataset")
      geneNetwork$Source <- substr(geneNetwork$Source, 0, nchar(as.character(geneNetwork$Source)) -3)
      
      # Include any interactions where TFs from core network are source or targets
      coreInteractions <- geneNetwork[which(geneNetwork$Source %in% tfInNetwork),]
      coreInteractionsTmp <- geneNetwork[which(geneNetwork$Target %in% tfInNetwork),]
      
      coreInteractions <- rbind(coreInteractions, coreInteractionsTmp)
      colnames(coreInteractions) <- c("Source", "Target", "Dataset")
      
      coreInteractions <- coreInteractions[!duplicated(coreInteractions),]
      
      seurat_regs <- readRDS(paste0("output/auc/seurat_regs",c,"_",t,"_",direc,".Rds"))
      networkGenes <- union(coreInteractions$Source, coreInteractions$Target)
      tmp <- seurat_regs@assays$RNA@scale.data
      regNames <- substr(rownames(tmp), 0, nchar(as.character(rownames(tmp))) -3)
      networkGenes <- networkGenes[which(networkGenes %in% regNames)]
      coreInteractions <- coreInteractions[which(coreInteractions$Source %in% networkGenes),]
      coreInteractions <- coreInteractions[which(coreInteractions$Target %in% networkGenes),]
      networkGenes <- union(coreInteractions$Source, coreInteractions$Target)
      rownames(tmp) <- regNames
      tmp <- tmp[which(rownames(tmp) %in% networkGenes),]
      tmp <- t(tmp)

      actMi <- calculateMI(data = tmp)

      mi <- integer(length = length(coreInteractions$Source))
      for(int in seq_along(coreInteractions$Source)){
        mi[int] <- actMi[as.character(coreInteractions$Source[int]), as.character(coreInteractions$Target[int])]
      }
      coreInteractions$Mi <- mi
      
      actCor <- cor(tmp, method = "s")
      type <- integer(length = length(coreInteractions$Source))
      for(int in seq_along(coreInteractions$Source)){
        type[int] <- actCor[as.character(coreInteractions$Source[int]), as.character(coreInteractions$Target[int])]
      }
      coreInteractions$Cor <- type
      type[type>0] <- 1
      type[type<=0] <- 2
      coreInteractions$Type <- type
      coreInteractions$Dataset <- treatmentName
      #  }
      fullNetworkList[[treatmentName]] <- coreInteractions
      counter <- counter+1
    }
  }
}
saveRDS(fullNetworkList, file = paste0(working_dir,"fullNetworkList_top23RegAllRegAllInt.rds"))
        
srcGenes <- lapply(fullNetworkList, function(x)x$Source)
tgtGenes <- lapply(fullNetworkList, function(x)x$Target)
dataset <- lapply(fullNetworkList, function(x)x$Dataset)
types <- lapply(fullNetworkList, function(x)x$Type)
actCor <- lapply(fullNetworkList, function(x)x$Cor)
mi <- lapply(fullNetworkList, function(x)x$Mi)
cellLine <- str_match(unlist(dataset), "(.*?)_")[,2]
signal <- str_match(unlist(dataset), "_(.*?)_")[,2]

allNetworkFull <- data.frame(Source = unlist(srcGenes), Target = unlist(tgtGenes), Type = unlist(types), Dataset = unlist(dataset), Cor = unlist(actCor), Mi = unlist(mi), CellLine = cellLine,Signal = signal, stringsAsFactors = FALSE)
saveRDS(allNetworkFull, paste0(working_dir,"allNetworkFullCoreMiTop23AllRegAllInt.rds"))
