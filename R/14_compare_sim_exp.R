# Compare experimental & simulated data - here we take the example of one generated network for 
# OVCA420_TGFB1 with posMi and negMi cutoffs 0.4 and 0.1, respectively

rm(list=ls())
library(sRACIPE)

results_dir <- file.path(getwd(),"results")
networks_dir <- file.path(results_dir,"network_construction")
sim_dir <- file.path(results_dir,"simulations")

treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
directions <- c("fwd","rev")


c <- cancer_types[4]
t <- treatments[2]
posMi <- 0.4
negMi <- 0.1

# Read the simulated data from the simulations above.
expRacipe <- readRDS(file.path(sim_dir,paste0(c,"_",t,"_",posMi,"_",negMi,".rds")))
circuit <- sracipeCircuit(expRacipe)
networkTFs <- union(circuit$Source,circuit$Target)
expRacipeNorm <- sracipeNormalize(expRacipe)
expData <- assay(expRacipeNorm,1)

expDataOrig <- assay(expRacipe,1)
expDataOrig <- log2(1+expDataOrig)
#expDataOrig <- expDataOrig[commonTfs,]

tmpMeans <- rowMeans(expDataOrig)
tmpSds <- apply(expDataOrig,1,sd)

seuratFwd <- readRDS(file = file.path(results_dir,paste0("seurat_regs",c,"_",t,"_fwd.Rds")))
actFwd <- seuratFwd@assays$RNA@scale.data
tfNames <- rownames(actFwd)
tfNames <- substr(tfNames, 0, nchar(tfNames)-3)
rownames(actFwd) <- substr(rownames(actFwd), 0, nchar(rownames(actFwd))-3)
commonTfs <- intersect(networkTFs,rownames(actFwd))

seuratRev <- readRDS(file =file.path(results_dir,paste0("seurat_regs",c,"_",t,"_rev.Rds")))
actRev <- seuratRev@assays$RNA@scale.data
tfNames <- rownames(actRev)
tfNames <- substr(tfNames, 0, nchar(tfNames)-3)
rownames(actRev) <- substr(rownames(actRev), 0, nchar(rownames(actRev))-3)
commonTfs <- intersect(commonTfs,rownames(actRev))

expData <- expData[commonTfs,]
actFwd <- actFwd[commonTfs,]
actRev <- actRev[commonTfs,]

plotColor <- c("#5E4FA2", "#4F61AA", "#4173B3", "#3386BC", "#4198B6",
               "#51ABAE", "#62BEA6", "#77C8A4", "#8ED1A4", "#A4DAA4",
               "#B8E2A1", "#CBEA9D", "#DEF199", "#EAF69F", "#F2FAAC",
               "#FAFDB8", "#FEFAB6", "#FEF0A5", "#FEE695", "#FDD985",
               "#FDC978", "#FDB96A", "#FCA75E", "#F99254", "#F67D4A",
               "#F26943", "#E85A47", "#DE4B4B", "#D33C4E", "#C1284A",
               "#AF1446", "#9E0142")

pca = prcomp(t(expData), scale. = T, center = T)
pcaData <- pca$x
pcaDataFwd <- (scale(t(actFwd), pca$center, pca$scale) %*%
                 pca$rotation)
pcaDataFwd <- data.frame(PC1=pcaDataFwd[,1],PC2=pcaDataFwd[,2], Time = seuratFwd$Time)
pcaDataRev <- (scale(t(actRev), pca$center, pca$scale) %*%
                 pca$rotation)
pcaDataRev <- data.frame(PC1=pcaDataRev[,1],PC2=pcaDataRev[,2], Time = seuratRev$Time)
expDataSim <- data.frame(pca$x[,1:2])
expDataSim$Time <- "Sim"
pcaDataFwd$Dir <- "Fwd"
pcaDataRev$Dir <- "Rev"
expDataSim$Dir <- "None"

# Fit a linear model to activities from cells from day 7 from forward and
# reverse directions.

sevenD <- rownames(pcaDataFwd[which(pcaDataFwd$Time == "7d"),])
sevenDFwd <- pcaDataFwd[sevenD,]$PC1
sevenDRev  <- pcaDataRev[sevenD,]$PC1
pred <- pcaDataRev$PC1
pcaDataRevScaled <- data.frame(PC1 = predict(lm(sevenDFwd ~ sevenDRev),data.frame(sevenDRev = pred)))

sevenDFwd <- pcaDataFwd[sevenD,]$PC2
sevenDRev  <- pcaDataRev[sevenD,]$PC2
pred <- pcaDataRev$PC2
pcaDataRevScaled$PC2 <- predict(lm(sevenDFwd ~ sevenDRev),data.frame(sevenDRev = pred))

pcaDataRevScaled$Time <- pcaDataRev$Time
pcaDataRevScaled$Dir <- pcaDataRev$Dir
allData <- rbind(expDataSim, pcaDataFwd,pcaDataRevScaled)

p <- ggplot2::ggplot(allData) +
  #    geom_point(aes(x = allData[,1], y =allData[,2], color = Time, shape = Dir), alpha = 0.5) +
  xlab(paste0("PC1(",100*summary(pca)$importance[2,1],"%)")) +
  ylab(paste0("PC2(",100*summary(pca)$importance[2,2],"%)")) +
  #    stat_summary(aes(x = allData[,1], y =allData[,2], group=Dir, color = Dir), fun.y=mean, geom="line") +
  geom_point(aes(x = allData[,1], y =allData[,2], group=Dir, color = Time, shape = Dir), alpha = 0.5) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  xlim(-5,11) +
  ylim(-5,4)
p


