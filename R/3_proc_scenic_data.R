## 3 - processing pySCENIC results: this script reads as input the results from pySCENIC and puts them into R-friendly data formats for downstream analysis. 
## It also conducts some basic data processing: regulon activity is scaled and highly variable TFs are identified
## Finally, it performs dimensionality reduction via PCA and UMAP

rm(list = ls()) # clean workspace
library(Seurat)

options(stringsAsFactors = FALSE)
## Network Construction
treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
direction <- c("fwd", "rev")

data_dir <- file.path(getwd(),"data")
results_dir <- file.path(getwd(),"results")

# Iterate through each condition
for(tCounter in 1:length(treatments)){
  t = treatments[tCounter]
  for(c in 1:length(cancer_types)) {
    
    ## Read in seurat object and regulon data
    ## We will use regulon activities here in place of counts and perform routine Seurat data processing
    # First, process data from fwd direction 
    seurat <- readRDS(file.path(data_dir,paste0(cancer_types[c],"_",t,".rds")))
    auc_fname <- file.path(results_dir,paste0("AUC_",cancer_types[c],"_",t,"_fwd.csv"))
    auc_matrix <- read.table(auc_fname, header=TRUE, sep = ",", row.names = 1, check.names = FALSE)
    auc_matrix <- t(auc_matrix)
    seuratFwd <- subset(seurat, Time %in% c("0d", "1d", "3d", "7d", "8h"))
    auc_matrix <- auc_matrix[,colnames(seuratFwd)]
    seurat_regs <- CreateSeuratObject(counts = auc_matrix)
    seurat_regs <- AddMetaData(seurat_regs, metadata = FetchData(seurat, "Time"))
    Idents(seurat_regs) <- "Time"
    allGenes <- rownames(seurat_regs)
    
    # Scale data, find highly variable genes
    seurat_regs <- ScaleData(seurat_regs, features = allGenes)
    seurat_regs <- FindVariableFeatures(seurat_regs)
    
    # Dimensionality reduction w/ PCA and UMAP
    seurat_regs <- RunPCA(seurat_regs)
    seurat_regs <- RunUMAP(seurat_regs, dims = 1:50)
   
    # Save Rds files
    saveRDS(seurat_regs, file = file.path(results_dir,paste0("seurat_regs",cancer_types[c],"_",t,"_fwd.Rds")))
    
    ## Next, process the rev direction the same way
    auc_fname <- file.path(results_dir,paste0("AUC_",cancer_types[c],"_",t,"_rev.CSV"))
    auc_matrix <- read.table(auc_fname, header=TRUE, sep = ",", row.names = 1, check.names = FALSE)
    auc_matrix <- t(auc_matrix)
    seuratRev <- subset(seurat, Time %in% c("1d_rm","3d_rm",  "7d", "8h_rm" ))
    auc_matrix <- auc_matrix[,colnames(seuratRev)]
    
    seurat_regs <- CreateSeuratObject(counts = auc_matrix)
    seurat_regs <- AddMetaData(seurat_regs, metadata = FetchData(seurat, "Time"))
    Idents(seurat_regs) <- "Time"
    allGenes <- rownames(seurat_regs)
    
    seurat_regs <- ScaleData(seurat_regs, features = allGenes)
    seurat_regs <- FindVariableFeatures(seurat_regs)
    
    seurat_regs <- RunPCA(seurat_regs)
    seurat_regs <- RunUMAP(seurat_regs, dims = 1:50)
    
    
    saveRDS(seurat_regs, file = file.path(results_dir,paste0("seurat_regs",cancer_types[c],"_",t,"_rev.Rds")))
  }
}
