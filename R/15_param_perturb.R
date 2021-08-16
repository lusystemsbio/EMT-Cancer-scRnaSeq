rm(list=ls())
library(sRACIPE)
library(ggplot2)

results_dir <- file.path(getwd(),"results")
networks_dir <- file.path(results_dir,"network_construction")
sim_dir <- file.path(results_dir,"simulations")
perturb_dir <- file.path(results_dir,"perturbation_simulations")
if(!dir.exists(perturb_dir)) {
  dir.create(perturb_dir)
}

treatments <- c("EGF","TGFB1", "TNF")
cancer_types <- c("A549","DU145","MCF7","OVCA420")
directions <- c("fwd","rev")


c <- cancer_types[4]
t <- treatments[2]
posMi <- 0.4
negMi <- 0.1

# Read the simulated data from the simulations above.
newRacipe <- readRDS(file.path(sim_dir,paste0(c,"_",t,"_",posMi,"_",negMi,".rds")))
#newRacipe <- readRDS(file.path("eModelsRacipe.rds"))
tmp <- sracipeParams(newRacipe)
tmp[,"G_JUN"] <- 10000*tmp[,"G_JUN"]
tmp[,"G_RELB"] <- 10000*tmp[,"G_RELB"]
sracipeParams(newRacipe) <- tmp
tmp <- assay(newRacipe,1)
sracipeIC(newRacipe) <- (tmp)
set.seed(123)
newRacipe2_D05_T2000 <-  sracipeSimulate(
  newRacipe, numModels = 2000, initialNoise = 0.05,simDet = FALSE,
  nNoise=0,stepper = "EM", plots = FALSE,integrateStepSize = 0.01,
  genParams = FALSE,genIC = FALSE, simulationTime = 2000)

saveRDS(newRacipe2_D05_T2000, file = file.path(perturb_dir,"eModelsRacipe_D05_T2000.rds"))

# We can instead also record the state of the models at multiple time points,
# for example, here we record the state at every 0.02 time. 

newRacipeTime <-  sracipeSimulate(
  newRacipe, numModels = 2000, initialNoise = 0.05,simDet = FALSE,
  nNoise=0,stepper = "EM", plots = FALSE,integrateStepSize = 0.01,
  genParams = FALSE,genIC = FALSE, simulationTime = 50, printInterval = 0.02, printStart = 0)

saveRDS(newRacipeTime, file = file.path(perturb_dir,"eModelsRacipe_D05_T50_interval02.rds"))



