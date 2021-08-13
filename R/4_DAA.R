## 4 - Differential activity analysis: this script will identify regulons with differential levels of activity between timepoints.
## Specifically, it makes 7 comparisons for each of 12 conditions (84 total): 0d vs 8h, 1d, 3d, and 7d; and 7d vs 8h_rm, 1d_rm, 3d_rm

rm(list = ls()) # clean workspace
library(Seurat)

options(stringsAsFactors = FALSE)

# Global variables
treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
direction <- c("fwd", "rev")
results_dir <- file.path(getwd(),"results")

# Compile a set of all DARs across all conditions
# Note: since we placed the regulon activities in the "counts" slot, this analysis will yield differentially active regulons, 
# not differentially expressed genes (DEGs). 
degs <- data.frame()
for(tCounter in seq_along(treatments)){
  for(cCounter in seq_along(cancer_types)){
    for(dCounter in seq_along(direction)){
      
      ## Select a treatment to work with (default is iterate thru all conditions in a loop)
      t <- treatments[tCounter]
      c <- cancer_types[cCounter]
      direc <- direction[dCounter]
      treatmentName <- paste0(c,"_",t,"_",direc)
      print(treatmentName)
      
      seurat_regs <- readRDS(file.path(results_dir,paste0("seurat_regs",c,"_",t,"_",direc,".Rds")))
      
      # For fwd, identify DEGs between 0d and each later timepoint
      if(direc == "fwd"){
        for(degTime in c("8h", "1d", "3d", "7d")){
          deg <- FindMarkers(seurat_regs, ident.1 = degTime, ident.2 = "0d", logfc.threshold = 0.0)
          deg$CellLine <- c
          deg$gene <- rownames(deg)
          deg$Comparison <- paste0(c,"_",t,"_",direc)
          deg$Time <-  degTime
          degs <- rbind(degs,deg)
        }
      }
      # For rev, identify DEGs between 7d and each subsequent timepoint after signal removal
      if(direc == "rev"){
        for(degTime in c("8h_rm", "1d_rm", "3d_rm")){
          deg <- FindMarkers(seurat_regs, ident.1 = degTime, ident.2 = "7d", logfc.threshold = 0.0)
          deg$CellLine <- c
          deg$gene <- rownames(deg)
          deg$Comparison <- paste0(c,"_",t,"_",direc)
          deg$Time <-  degTime
          degs <- rbind(degs,deg)
        }
      }
    }
  }
}

# Save compiled results to Rds file
saveRDS(degs, file = file.path(results_dir,"degRegsFwd0dRev7d.rds"))
