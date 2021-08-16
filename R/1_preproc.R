# Prepare data for processing
# First, install the needed packages from CRAN
install.packages(c("stringr", "jsonlite", "sjmisc", "MASS", "UpSetR", "ggplot2", "googledrive","dplyr","Binarize"))
install.packages(c("varhandle", "plyr", "ggpmisc", "spatstat.utils", "umap", "ggthemes","RColorBrewer","infotheo", "igraph")) # rowr?

install.packages("BiocManager")
BiocManager::install("Seurat")

library(remotes)
remotes::install_version("spatstat", version = "1.64-1")

rm(list = ls()) # clean workspace

library(Seurat)
library(stringr)
library(sjmisc)
library(MASS)
library(googledrive)

cancers <- c("A549","DU145","MCF7","OVCA420")
treatments <- c("EGF","TGFB1","TNF")

fIDs <- c("1BJePXMkn_-iOsanOkx5fJO5vDVpLgfLe",       # A549_EGF.rds
          "19aRaQ29uWNxQOSlyuzxEychHNseFK-kR",      # A549_TGFB1.rds
          "1ffnfyx1siBFzgSJENx98xvHdeqOifXH5",     # A549_TNF.rds
          "1mLIBRFrLvvcjBFoX6vBZ4ehiJ1zxIaVF",     # DU145_EGF.rds
          "1uMNae00bJ8AfP4u9IwZYF29S67tAz0Tj",     # DU145_TGFB1.rds
          "1cKo8xXB_rEG6JI3rMXDIQHtDhAgfns92",     # DU145_TNF.rds
          "130Mao1De6lGaQSw9_hzLXuW9cfPROX7z",     # MCF7_EGF.rds
          "1Jd8JvfaIvq4yqLkgbNRPj5OBGBd0J4Vu",     # MCF7_TGFB1.rds
          "1K43TgiAh_Chql6XPKud0oOhr4c_s_sxU",     # MCF7_TNF.rds
          "1J2oUOJlTeZDCX-1A5akGCA6jofUj1Tfn",     # OVCA420_EGF.rds
          "1dnN6BoZP3GZj62_GAEWEqwgohpRVbO3k",     # OVCA420_TGFB1.rds
          "1O4Yt-qaZFX1is3uT_OZqYJmrKukJrJUA"      # OVCA420_EGF.rds
)

dir <- "data" # Directory where the downloaded rds files are stored. From: https://drive.google.com/drive/folders/1lZ38Uj2ZjmFus7XbHGTATh8f9MqXLAf_
if(!dir.exists(dir)) {
  dir.create(dir)
}

# Loop through google drive IDs above & download data
# You must run the function drive_auth() to login before downloading files
treatmentIDs <- rep(c(1:3),4)
for(fID in seq_along(fIDs)) {
  cID <- ceiling(fID / 3)
  tID <- treatmentIDs[fID]
  destName <- paste0(cancers[cID],"_",treatments[tID],".rds")
  drive_download(as_id(fIDs[fID]), path = file.path(dir,destName), overwrite = T)
}

files <- list.files(dir, pattern = ".\\.rds")

for(i in seq_along(files)){
  filename <- file.path(dir,files[i])
  data <- readRDS(filename)
  fwd <- subset(data, Time %in% c("0d", "1d", "3d", "7d", "8h"))
  tmp <- as.matrix(Seurat::GetAssayData(fwd[["RNA"]], slot = "data"))
  write.table(tmp,row.names = TRUE,col.names = TRUE, sep = "\t", file = file.path(dir,paste0(tools::file_path_sans_ext(files[i]),"_fwd.tsv")),quote = FALSE)
  rm(tmp,fwd)
  rev <- subset(data, Time %in% c("8h_rm", "1d_rm","3d_rm", "7d"))
  tmp <- as.matrix(Seurat::GetAssayData(rev[["RNA"]], slot = "data"))
  write.table(tmp,row.names = TRUE,col.names = TRUE, sep = "\t", file = file.path(dir,paste0(tools::file_path_sans_ext(files[i]),"_rev.tsv")),quote = FALSE)
  
}

## Next, we'll set up the directory structure for the python script
resultsDir <- file.path(getwd(),"results")
dir.create(resultsDir)
resourceDir <- file.path(getwd(),"resources")
dir.create(resourceDir)
databaseDir <- file.path(getwd(),"databases")
dir.create(databaseDir)

# download a list of known TFs in the human genome
download.file("https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_curated_tfs.txt",
              destfile = file.path(resourceDir,"hs_hgnc_curated_tfs.txt"))

# download motif databases from CisTarget (resources.aertslab.org/cistarget)
# here, we use the 7-species databases, using 500bpUp, TSS+/-10kbp, and TSS+/-5kbp
# the options command can adjust the default time-out parameter - this may need to be increased w/ slow internet connection
options(timeout = 1800)
download.file("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
              destfile = file.path(databaseDir,"hg19-500bp-upstream-7species.mc9nr.feather"))

download.file("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather",
              destfile = file.path(databaseDir,"hg19-tss-centered-10kb-7species.mc9nr.feather"))

download.file("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-7species.mc9nr.feather",
              destfile = file.path(databaseDir,"hg19-tss-centered-5kb-7species.mc9nr.feather"))

# download motif annotations from CisTarget
download.file("https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
              destfile = file.path(databaseDir,"motifs-v9-nr.hgnc-m0.001-o0.0.tbl"))
options(timeout = 60)



