## 5 - Identifying differentially active TFs across several experimental conditions

rm(list = ls()) # clean workspace
library(Seurat)
library(jsonlite)
library(sjmisc)
library(stringr)
library(MASS)
library(UpSetR)
library(ggplot2)
library(dplyr)

options(stringsAsFactors = FALSE)

## Global variables
treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
direction <- c("fwd", "rev")
results_dir <- file.path(getwd(),"results")

# Make a list of all conditions
counter <- 1
treatmentNames <- c()
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    t <- treatments[tCounter]
    c <- cancer_types[cCounter]
    treatmentNames <- c(treatmentNames,paste0(c,"_",t))
    
  }
}

# Read in complete set of DEGs
degs <- readRDS(file = file.path(results_dir, "degRegsFwd0dRev7d.rds"))
degs$Signal <- str_match(degs$Comparison, "_(.*?)_")[,2]

# Loop through all conditions 
diffRegulonList <- list()
nDiffRegulon <- 100
set.seed(123)
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    for(dCounter in seq_along(direction)){
      t <- treatments[tCounter]
      c <- cancer_types[cCounter]
      direc <- direction[dCounter]
      treatmentName <- paste0(c,"_",t,"_",direc)
      print(treatmentName)
      
      # Identify DEGs from this condition 
      deg_reg <- degs[degs$Comparison == paste0(c,"_",t,"_",direc),] 
      seurat_regs <- readRDS(file.path(results_dir, paste0("seurat_regs",c,"_",t,"_",direc,".Rds")))
      
      # Read regulons from SCENIC for this condition
      regulon_fname <- file.path(results_dir,paste0("regulons_",c,"_",t,"_",direc,".json"))
      regs <- read_json(regulon_fname, simplify = T)
      
      # Take the top N regulons (100 by default) by adjusted p-value
      networkReg <- deg_reg$gene[rank(deg_reg$p_val_adj, ties.method = "random")]
      networkReg <- unique(networkReg)
      networkReg <- networkReg[1:min(nDiffRegulon,length(networkReg))]
      regs <- regs[networkReg]
      diffRegulonList[[treatmentName]]  <-  regs
    }
  }
}
saveRDS(diffRegulonList, file = file.path(results_dir, "diffRegulonTop100Deg0d7d.rds"))

# Combine regulon sets for fwd and rev into one list for each condition
diffRegulonNames <- lapply(diffRegulonList, function(x) names(x))
diffRegulonFwdRev <- list()
for(i in 1:(length(diffRegulonList)/2)){
  diffRegulonFwdRev[[i]] <- union(diffRegulonNames[[2*i-1]],diffRegulonNames[[2*i]])
}

# Save combined regulon lists to file
names(diffRegulonFwdRev) <- treatmentNames
saveRDS(diffRegulonFwdRev, file = file.path(results_dir, "diffRegulonFwdRev.rds"))

# Visualize overlap in differentially active regulons with upSet plots
diffRegulonFwdRev = diffRegulonFwdRev[order(names(diffRegulonFwdRev))]
upset(fromList(diffRegulonFwdRev),keep.order = F, group.by = "degree", text.scale = 2, nsets = length(diffRegulonFwdRev)) 

# Histogram of how many comparisons each TF is differentially active across
degs <- readRDS(file = file.path(results_dir, "degRegsFwd0dRev7d.rds"))
degsFiltered <- degs[degs$p_val_adj < 0.05,]
TfFrequency <- degsFiltered %>%
  group_by(gene) %>%
  summarize(n())
hist(TfFrequency$`n()`, breaks = 20)

# Select all DATFs which appear in at least 24/84 comparisons as the core network
coreTf <- TfFrequency$gene[TfFrequency$`n()` > 23]
saveRDS(coreTf, file = file.path(results_dir,"coreTf.rds"))
