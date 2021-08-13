## 6 - plot activities: this script visualizes the regulon activity across conditions using UMAP
# It will create a 12-page pdf with UMAP plots of regulatory activity for each condition
# Each timepoint will be marked by a black circle, which increase in opacity and size as time progresses, and are connected by a dotted line.
## Expected runtime: 


### Setup
rm(list = ls())
library(Seurat)
library(jsonlite)
library(sjmisc)
library(stringr)
library(varhandle)
library(plyr)
library(ggplot2)
library(ggpmisc)
library(spatstat.utils)
library(umap)
library(ggthemes)
library(RColorBrewer)

# Install a forked version of the LSD package, for heatscatter plots compatible with ggplot
dir.create(file.path(getwd(),"LSD_dir"))
remotes::install_github("stineb/LSD", lib=file.path(getwd(),"LSD_dir"))
library(LSD, lib.loc = file.path(getwd(),"LSD_dir"))


## Initialize some values
treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
directions <- c("fwd","rev")
results_dir <- file.path(getwd(),"results")
fig_dir <- file.path(getwd(),"figs")
if(!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

plotFName <- file.path(fig_dir,"UMAP_density_alpha.pdf")
working_directory <- getwd()
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E84855")

## Color setup for plots
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot_color <- rf(20)

## Make plot list
umap_density_list <- list()
iter <- 1

## Loop thru cell lines, treatments, directions
for(cancer in seq_along(cancer_types)) {
  c <- cancer_types[cancer]
  for(treatment in seq_along(treatments)) {
    t <- treatments[treatment]
    shared_regs <- c()
    print(paste0("Starting work on ",c,"_",t))
    for(dir in seq_along(directions)) {
      d <- directions[dir]
      if(d == "fwd") {
        seurat_regs_fwd <- readRDS(file.path(results_dir,paste0("seurat_regs",c,"_",t,"_",d,".Rds")))
        auc_mtx_fwd <- as.data.frame(t(GetAssayData(seurat_regs_fwd, slot = "scale.data")))
        auc_mtx_fwd$Time <- FetchData(seurat_regs_fwd, "Time")[,"Time"]
        auc_mtx_fwd$Direction <- d
        auc_mtx_fwd$Cell <- rownames(auc_mtx_fwd)
        
      } else {
        seurat_regs_rev <- readRDS(file.path(results_dir,paste0("seurat_regs",c,"_",t,"_",d,".Rds")))
        auc_mtx_rev <- as.data.frame(t(GetAssayData(seurat_regs_rev, slot = "scale.data")))
        auc_mtx_rev$Time <- FetchData(seurat_regs_rev, "Time")[,"Time"]
        auc_mtx_rev$Direction <- d
        auc_mtx_rev$Cell <- rownames(auc_mtx_rev)
        
      }
    }
    
    ## Scale using 7d data
    shared_regs <- colnames(auc_mtx_fwd)[which(colnames(auc_mtx_fwd) %in% colnames(auc_mtx_rev))]
    shared_regs_short <- shared_regs[-((length(shared_regs)-2):length(shared_regs))]
    auc_mtx_fwd_scale <- auc_mtx_fwd[,shared_regs]
    auc_mtx_rev_scale <- auc_mtx_rev[,shared_regs]
    for(reg in shared_regs_short) {
      fwd_7d_reg <- auc_mtx_fwd_scale[which(auc_mtx_fwd_scale$Direction == "fwd" & 
                                              auc_mtx_fwd_scale$Time == "7d"),reg]
      rev_7d_reg <- auc_mtx_rev_scale[which(auc_mtx_rev_scale$Direction == "rev"& 
                                              auc_mtx_rev_scale$Time == "7d"),reg]
      ## Generate linear models
      x_lm <- lm(fwd_7d_reg ~ rev_7d_reg)
      new_rev_reg <- (auc_mtx_rev_scale[which(auc_mtx_rev_scale$Direction == "rev" ),reg] * 
                        unname(x_lm$coefficients["rev_7d_reg"])) + 
        unname(x_lm$coefficients["(Intercept)"])
      
      ## Replace reverse data
      auc_mtx_rev_scale[,reg] <- new_rev_reg
      
    }
    
    ## Combine fwd and rev data
    auc_mtx_fwd_scale <- auc_mtx_fwd_scale[which(!auc_mtx_fwd_scale$Time == "7d"),]
    combined_mtx <- rbind(auc_mtx_fwd_scale, auc_mtx_rev_scale)
    na_regs <- colnames(combined_mtx[,colSums(is.na(combined_mtx)) != 0])
    combined_mtx_noNA <- combined_mtx[ , colSums(is.na(combined_mtx)) == 0]
    pca_regs <- shared_regs_short[which(shared_regs_short %in% colnames(combined_mtx_noNA))]
    combined_cellInfo <- combined_mtx[,c("Time","Direction",'Cell')]
    
    ## Progress report
    print(paste0("Finished scaling data for ",c,"_",t,". Starting dimensional reduction."))
    
    ## Run UMAP
    umap_rs <- umap(combined_mtx_noNA[,pca_regs])
    umap_data <- as.data.frame(umap_rs$layout)
    colnames(umap_data) <- c("x","y")
    if(all.equal(combined_cellInfo$Cell, rownames(umap_data))) {
      umap_data$Time <- combined_cellInfo$Time
    }
    umap_data$Time <- factor(umap_data$Time, levels = c("0d", "8h", "1d", "3d", "7d", "8h_rm", "1d_rm", "3d_rm"))
    
    ## Progress report
    print(paste0("Done with dimensional reductions for ",c,"_",t,". Creating plots now."))
    
    ## Calculate centroids for each timepoint
    centroid_df_umap <- aggregate(umap_data[,c("x","y")], list(Time=umap_data$Time), mean)
    
    ## Plot UMAP density with timepoint data
    ggscatter <- heatscatter(umap_data[,1], umap_data[,2], 
                             xlab = "UMAP1", 
                             ylab= "UMAP2", 
                             ggplot = TRUE)
    image <- ggscatter + 
      geom_point(data=centroid_df_umap, mapping = aes(x=x, y=y, size=Time, alpha=Time)) +
      geom_path(data = centroid_df_umap, arrow = arrow(), linetype="dashed") +
      guides(alpha=FALSE, size=FALSE, color=FALSE, shape=FALSE) +
      scale_shape_manual(values = c(15,16,17,18,8,21,24,23)) +
      scale_alpha_discrete(breaks=c(0.4,0.485,0.57,0.655,0.74,0.825,0.91,1)) +
      ggtitle(paste(c," ",t)) +
      theme(title = element_text(size=20),
            axis.title = element_text(size = 20),
            axis.text = element_text(size=18),
            plot.subtitle = element_text(size=16),
            legend.title = element_text(size=20),
            legend.title.align = 0,
            legend.text = element_text(size=18),
            legend.key.size = unit(20, units = "pt")) +
      labs(subtitle = paste0("UMAP on ",length(pca_regs), " regulons in ", nrow(umap_data), " cells"))
    if(t == "TNF") {
      image <- image + 
        guides(size=guide_legend(title = "Time"), 
               color=guide_legend(override.aes = list(size = 3.5)), 
               shape=guide_legend(override.aes = list(size = 3.5)))
    }
    umap_density_list[[iter]] <- image
    
    
    ## Increase iterator
    iter <- iter+1
  }
}

## Save all plots
pdf(file = plotFName, width=11.5, height=8)
for(plot in umap_density_list) {
  print(plot, newpage = TRUE)
}
graphics.off()

