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

## ggplot theme for graphics
ggtheme <- function (base_size = 11, base_family = "") {
  theme_get() %+replace%
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size=14),
          axis.line = element_line(color = "black", linetype = "solid"),
          title = element_text(size = 16),
          legend.title = element_text(size=16),
          legend.title.align = 0,
          legend.text = element_text(size=14),
          legend.key.size = unit(20, units = "pt"),
          panel.grid = element_blank(),
          panel.background = element_blank())
}

# Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
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
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E84855")
fig_dir <- file.path(getwd(),"figs")
results_dir <- file.path(getwd(),"results")


## Get reg list
reg_list <- c("BHLHE40(+)", "CEBPB(+)", "CREB3(+)", "CREM(+)", "ELK3(+)",
              'ESRRA(+)', 'ETS1(+)', 'ETS2(+)', 'FOS(+)', 'FOSB(+)',
              'FOSL1(+)', 'IRF1(+)', 'IRF2(+)', 'IRF3(+)', 'IRF7(+)', 
              'JUN(+)', 'JUNB(+)', 'KLF6(+)', 'MAFG(+)', 'MYC(+)', 
              'NFE2L3(+)', 'NFKB2(+)', 'RELB(+)', 'SOX4(+)', 'SPDEF(+)',
              'STAT1(+)', 'STAT3(+)', 'XBP1(+)')
gene_list <- reg_list

## make dummy plot for null cases
dummy_df <- data.frame(x=1,y=1)
dummy_plot <- ggplot(data = dummy_df) +
  geom_point(aes(x=x, y=y, alpha=0.0001)) +
  guides(alpha=FALSE) +
  ggtheme() +
  geom_text(mapping = aes(x=1, y=1, label="NULL"))

for(g in seq_along(gene_list)) {
  gene_plotlist <- list()
  iter <- 1
  for(c in cancer_types) {
    for(t in treatments) {
      gene_df_full <- data.frame(Activity=numeric(), Time=character(),Direction=character())
      for(d in directions) {
        ## Read auc matrix, seurat object, de regs
        treatmentName <- paste0(c,"_",t,"_",d)
        seurat_regs <- readRDS(file.path(results_dir,paste0("seurat_regs",c,"_",t,"_",d,".Rds")))
        auc_mtx <- GetAssayData(seurat_regs, slot = "scale.data")
        
        ## Select gene of interest
        if(gene_list[g] %in% rownames(auc_mtx)) {
          gene_df <- as.data.frame(auc_mtx[gene_list[g],])
          colnames(gene_df)[1] <- "Activity"
          gene_df$Time <- FetchData(seurat_regs, "Time")[,"Time"]
          gene_df$Direction <- d
          gene_df_full <- rbind(gene_df_full, gene_df)
        }
      }
      ## Only proceed if there is data for the regulon - also correct color scheme if data for one direction is missing
      if(nrow(gene_df_full) > 0) {        
        gene_df_full$Time <- factor(gene_df_full$Time, levels = c("0d", "8h", "1d", "3d", "7d", "8h_rm", "1d_rm", "3d_rm"))
        
        if(length(unique(gene_df_full$Direction)) == 1) {
          if(unique(gene_df_full$Direction) == "fwd") {
            plotColors <- cbPalette
          } else {
            plotColors <- cbPalette[5:8]
          }
        } else {
          plotColors <- cbPalette
        }
        
        ## Plot
        image <- ggplot(gene_df_full, aes(x=Time,y=Activity, fill=Time)) +
          geom_boxplot() +
          ggtheme() +
          ylab(paste0(gene_list[g])) +
          scale_fill_manual(values=plotColors) +
          guides(fill=FALSE) +
          facet_wrap(~Direction, scales = "free") +
          theme(axis.text.y = element_text(size = 14), 
                title = element_text(size = 20), 
                axis.title.y = element_text(size = 16),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                strip.text = element_text(size = 16),
                axis.text.x = element_blank(),
                axis.title.x = element_blank()) +
          ggtitle(paste0(c," ", t))
        
        if(t == "TNF") {
          image <- image +
            theme(axis.text.x = element_text(size = 16, angle = 90), 
                  axis.title.x = element_text(size=20))
        }
        gene_plotlist[[iter]] <- image
        iter <- iter+1
      } 
      else {
        gene_plotlist[[iter]] <- dummy_plot
        iter <- iter+1
        
      }
    }
  }
  
  ## Save plot for this regulon
  title <- file.path(fig_dir,paste0(gene_list[g],"_time_activity.pdf"))
  pdf(file = title, width = 22, height = 17)
  multiplot(plotlist = gene_plotlist, cols = 4)
  dev.off()
  
}