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
#library(rowr)
library(spatstat.utils)
library(umap)
library(ggthemes)
library(RColorBrewer)

## ggplot theme for graphics
ggtheme <- function (base_size = 11, base_family = "") {
  theme_get() %+replace%
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size=16),
          axis.line = element_line(color = "black", linetype = "solid"),
          title = element_text(size = 20),
          legend.title = element_text(size=18),
          legend.title.align = 0,
          legend.text = element_text(size=16),
          legend.key.size = unit(20, units = "pt"),
          panel.grid = element_blank(),
          panel.background = element_blank())
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


## Initialize some values
treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
directions <- c("fwd","rev")
results_dir <- file.path(getwd(),"results")
fig_dir <- file.path(getwd(),"figs")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E84855")

## Parameters
genepair_lists <- list(c("RELB+NFKB2","FOSL1+JUNB"),c("FOSL1","RELB"),c("JUNB","RELB"))

## make dummy plot for null cases
dummy_df <- data.frame(x=1,y=1)
dummy_plot <- ggplot(data = dummy_df) +
  geom_point(aes(x=x, y=y, alpha=0.0001)) +
  guides(alpha=FALSE) +
  ggtheme() +
  geom_text(mapping = aes(x=1, y=1, label="NULL"))

for(genes in genepair_lists) {
  trans_centroid_plotlist <- list()
  iter <- 1
  
  g1_expanded <- strsplit(genes[1],"\\+")[[1]]
  g2_expanded <- strsplit(genes[2],"\\+")[[1]]
  full_genelist <- c(g1_expanded, g2_expanded)
  
  regs_g1 <-paste0(g1_expanded, "(+)") 
  regs_g2 <-paste0(g2_expanded, "(+)") 
  regs <- paste0(full_genelist, "(+)")
  
  ## Loop thru cancer types and treatments (12 conditions)
  for(c in cancer_types) {
    for(t in treatments) {
      activity_df_full <- data.frame(x=numeric(),y=numeric(),Time=character(),Direction=character())
      # Combine data for fwd and rev
      for(d in directions) {
        seurat_regs <- readRDS(file.path(results_dir,paste0("seurat_regs",c,"_",t,"_",d,".Rds")))
        if(length(regs[which(regs %in% rownames(seurat_regs))]) == length(regs)) {
          auc_mtx <- as.data.frame(t(GetAssayData(seurat_regs, slot = "scale.data")[regs,]))
          if(length(regs) == 2) {
            colnames(auc_mtx) <- c("x","y")
            auc_mtx$Time <- FetchData(seurat_regs, "Time")[,"Time"]
            auc_mtx$Direction <- d
          } 
          else if(length(regs) == 4){
            auc_mtx$x <- rowSums(auc_mtx[,c(regs_g1)], dims = 1)
            auc_mtx$y <- rowSums(auc_mtx[,c(regs_g2)], dims = 1)
            auc_mtx <- auc_mtx[,c("x","y")]
            auc_mtx$Time <- FetchData(seurat_regs, "Time")[,"Time"]
            auc_mtx$Direction <- d
          } 
          else {
            auc_mtx <- data.frame()
          }
          
        } else if(c == "MCF7" & t == "EGF") {
          ## Handle an exception case where data is missing for one regulon
          auc_mtx <- as.data.frame(t(GetAssayData(seurat_regs, slot = "scale.data")[c("RELB(+)",regs_g2),]))
          auc_mtx$x <- auc_mtx$`RELB(+)`
          auc_mtx$y <- rowSums(auc_mtx[,c(regs_g2)], dims = 1)
          auc_mtx <- auc_mtx[,c("x","y")]
          auc_mtx$Time <- FetchData(seurat_regs, "Time")[,"Time"]
          auc_mtx$Direction <- d
        }
        else {
          auc_mtx <- data.frame()
        }
        activity_df_full <- rbind(activity_df_full, auc_mtx)
        
      }
      ## If data for both directions is present, scale rev data to match fwd
      if(nrow(activity_df_full) > 0 & length(unique(activity_df_full$Direction)) == 2) {  
        ## Set colors 
        activity_df_full$Time <- factor(activity_df_full$Time, 
                                        levels = c("0d", "8h", "1d", "3d", "7d", "8h_rm", "1d_rm", "3d_rm"))
        plotColors <- cbPalette
        
        #### Scaled plots
        if(length(unique(activity_df_full$Direction)) == 2) {
          fwd_7d_x <- activity_df_full[which(activity_df_full$Direction == "fwd" & activity_df_full$Time == "7d"),"x"]
          fwd_7d_y <- activity_df_full[which(activity_df_full$Direction == "fwd"& activity_df_full$Time == "7d"),"y"]
          
          rev_7d_x <- activity_df_full[which(activity_df_full$Direction == "rev"& activity_df_full$Time == "7d"),"x"]
          rev_7d_y <- activity_df_full[which(activity_df_full$Direction == "rev"& activity_df_full$Time == "7d"),"y"]
          
          ## Generate linear models
          x_lm <- lm(fwd_7d_x ~ rev_7d_x)
          y_lm <- lm(fwd_7d_y ~ rev_7d_y)
          
          ## Rescale rev data
          new_rev_x <- (activity_df_full[which(activity_df_full$Direction == "rev" ),"x"] * 
                          unname(x_lm$coefficients["rev_7d_x"])) + 
            unname(x_lm$coefficients["(Intercept)"])
          new_rev_y <- (activity_df_full[which(activity_df_full$Direction == "rev" ),"y"] * 
                          unname(y_lm$coefficients["rev_7d_y"])) + 
            unname(y_lm$coefficients["(Intercept)"])
          
          ## Add to df
          activity_df_full$transformed_x <- activity_df_full$x
          activity_df_full[which(activity_df_full$Direction == "rev"),"transformed_x"] <- new_rev_x
          
          activity_df_full$transformed_y <- activity_df_full$y
          activity_df_full[which(activity_df_full$Direction == "rev"),"transformed_y"] <- new_rev_y
        }
        
        ## Make centroid df
        centroid_df <- aggregate(activity_df_full[,c("transformed_x","transformed_y")], 
                                 list(Time=activity_df_full$Time), mean)
        
        ## Plot centroids
        rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
        plot_color <- rf(20)
        image <- ggplot(centroid_df, aes(x=transformed_x,y=transformed_y)) +
          geom_point(mapping=aes(size=1,  color=Time)) +
          geom_path(aes(alpha=0.5), arrow = arrow(), linetype="dashed") +
          xlab(genes[1]) +
          ylab(genes[2]) +
          scale_size(range=c(1.75, 3)) +
          scale_shape_manual(values=c(16,17,3,4,8,18,15,13)) +
          guides(alpha=FALSE, size=FALSE, color=FALSE) + 
          scale_color_manual(values=plot_color) +
          ggtheme() +
          ggtitle(paste(c," ",t))
        trans_centroid_plotlist[[iter]] <- image
        iter <- iter+1
        
        
      } 
      else {
        trans_centroid_plotlist[[iter]] <- dummy_plot
        iter <- iter+1
      }
    }
  }
  
  title <- file.path(fig_dir,paste0(genes[1],"_",genes[2],"_activity_centroids_SCALED.pdf"))
  pdf(file = title, width = 22, height = 17)
  multiplot(plotlist = trans_centroid_plotlist, cols = 4)
  dev.off()
}
