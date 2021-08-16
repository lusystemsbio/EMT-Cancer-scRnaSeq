# 11 - Classify DATFs as E or M markers
# DATFs are considered M markers if they have a positive fold change in the fwd transition and negative in the rev transition
# E markers, conversely, should decrease in the fwd transition and subsequently increase in the rev transition
rm(list=ls())
library(igraph)
library(stringr)

options(stringsAsFactors = FALSE)
results_dir <- file.path(getwd(),"results")
networks_dir <- file.path(results_dir,"network_construction")
sim_dir <- file.path(results_dir,"simulations")

treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
directions <- c("fwd","rev")

allNetworkFull <- readRDS(file.path(networks_dir,"allNetworkFullCoreMiTop23AllRegAllInt.rds"))
allGenes <- union(allNetworkFull$Source,allNetworkFull$Target)
degs <- readRDS(file = file.path(results_dir,"degRegsFwd0dRev7d.rds"))
degs$Signal <- str_match(degs$Comparison, "_(.*?)_")[,2]


counter <- 1
treatmentNames <- c()
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    for(dCounter in seq_along(directions)){
      t <- treatments[tCounter]
      c <- cancer_types[cCounter]
      direc <- directions[dCounter]
      treatmentNames <- c(treatmentNames,paste0(c,"_",t,"_",direc))
      
    }
  }
}


# Classify DATFs as E or M based on whether they increase/decrease in the fwd and rev transitions, with p-value < 0.05
EMClass <- matrix(data = NA,ncol = length(treatmentNames),nrow = length(allGenes),dimnames = list(allGenes, treatmentNames))
for(tCounter in seq_along(treatmentNames)){
  degTmp <- degs[which((degs$Comparison == treatmentNames[tCounter]) & (degs$Time == "7d") & (degs$p_val_adj < 0.05)),]
  Mgenes <- degTmp$gene[which(degTmp$avg_log2FC>0)]
  Mgenes <- substr(Mgenes,0,nchar(Mgenes)-3)
  Mgenes <- Mgenes[which(Mgenes %in% allGenes)]
  Egenes <- degTmp$gene[which(degTmp$avg_log2FC<0)]
  Egenes <- substr(Egenes,0,nchar(Egenes)-3)
  Egenes <- Egenes[which(Egenes %in% allGenes)]
  EMClass[Egenes,treatmentNames[tCounter]] <- FALSE
  EMClass[Mgenes,treatmentNames[tCounter]] <- TRUE
  
  degTmp <- degs[which((degs$Comparison == treatmentNames[tCounter]) & (degs$Time == "3d_rm") & (degs$p_val_adj < 0.05)),]
  Mgenes <- degTmp$gene[which(degTmp$avg_log2FC<0)]
  Mgenes <- substr(Mgenes,0,nchar(Mgenes)-3)
  Mgenes <- Mgenes[which(Mgenes %in% allGenes)]
  Egenes <- degTmp$gene[which(degTmp$avg_log2FC>0)]
  Egenes <- substr(Egenes,0,nchar(Egenes)-3)
  Egenes <- Egenes[which(Egenes %in% allGenes)]
  EMClass[Egenes,treatmentNames[tCounter]] <- FALSE
  EMClass[Mgenes,treatmentNames[tCounter]] <- TRUE
}
saveRDS(EMClass, file = file.path(results_dir,"EMClassP05.rds"))


counter <- 1
conditionNames <- c()
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    t <- treatments[tCounter]
    c <- cancer_types[cCounter]
    conditionNames <- c(conditionNames,paste0(c,"_",t))
  }
}

## Unify DATF classifications across fwd and rev direction for each condition
# If classifications disagree, leave DATFs as NA; if a classification is only present for one direction, use that for the overall condition
EMClassCondition <- matrix(data = NA,ncol = length(conditionNames),nrow = length(allGenes),dimnames = list(allGenes, conditionNames))
for(tCounter in seq_along(conditionNames)){
  for(geneCounter in seq_along(allGenes)){
    if(is.na(EMClass[allGenes[geneCounter],2*tCounter-1])){
      if(!is.na(EMClass[allGenes[geneCounter],2*tCounter])){
        EMClassCondition[allGenes[geneCounter],conditionNames[tCounter]] <- EMClass[allGenes[geneCounter],2*tCounter]
      }
    }
    if(!is.na(EMClass[allGenes[geneCounter],2*tCounter-1])){
      if(is.na(EMClass[allGenes[geneCounter],2*tCounter])){
        EMClassCondition[allGenes[geneCounter],conditionNames[tCounter]] <- EMClass[allGenes[geneCounter],2*tCounter-1]
      }
      if(!is.na(EMClass[allGenes[geneCounter],2*tCounter])){
        if(EMClass[allGenes[geneCounter],2*tCounter-1] == EMClass[allGenes[geneCounter],2*tCounter]){
          EMClassCondition[allGenes[geneCounter],conditionNames[tCounter]] <- EMClass[allGenes[geneCounter],2*tCounter-1]}
        if(!(EMClass[allGenes[geneCounter],2*tCounter-1] == EMClass[allGenes[geneCounter],2*tCounter])){
          EMClassCondition[allGenes[geneCounter],conditionNames[tCounter]] <- NA}
      }
    }
  }
}
saveRDS(EMClassCondition, file.path(results_dir,"EMClassCondition.rds"))


