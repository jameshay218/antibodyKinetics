######################
## JAMES HAY 11.02.2019 - jameshay218@gmail.com
## This script simulates antibody trajectories from the "best" fitting models,
## simulating data using the MLE parameters and different sampling strategies

library(ggplot2)
library(coda)
setwd("~/net/home/ferret")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
save_wd <- "~/net/home/ferret/inputs_Jan2019/sim_best_model/"

## Which model variants to simulate?
sim_runs <- c("CNAY6BY", "CYTY6BN", "CYAY6BY", "CNTY6BY",  
              "CNTY6BN", "CNTY6MY", "CYTY6BY", "CYAY6BN",  
              "CNAY6MY", "CYAY6MY",  "CNAY6BN", "CYTY6MY",
              "CYAY6MN")

## Where to read in parameter and exposure tables
runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs$parTab_file <- paste0("inputs/parTabs/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")

## Subset runs only by those we want to simulate from
runs <- runs[runs$runName %in% sim_runs, ]
all_runs_sim <- NULL
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

adaptive <- adaptive_period <- 1000000

## For each of these models, simulate data with each sampling protocol
for(i in 1:nrow(runs)){
  print(runs$runID[i])
  run <- runs$runName[i]
  times <- c(0,21,37,49,70)
  
  parTab_file <- runs$parTab_file[i]
  exposureTab_file <- runs$exposureTab_file[i]
  created_data_runs <- runs[runs$runName == run,]
  save_wd1 <- paste0(save_wd, run)
  if(!dir.exists(save_wd1)) dir.create(save_wd1,recursive=TRUE)
 
  parTab <- read.csv(runs$parTab_file[i],stringsAsFactors=FALSE)
  chain <- quiet(as.data.frame(load_mcmc_chains(paste0("/media/james/Storage 2/ferret_results/rerun_Jan2019/outputs_real/",
                                                       runs$runID,"_",runs$runName[i]), 
                                                parTab, FALSE, 1, adaptive, FALSE, TRUE,TRUE)[["chain"]]))
  pars <- get_best_pars(chain)
  pars <- pars[1:(length(pars)-1)]
  
  #print(pars)
  pars[which(names(pars) == "m")] <- sapply(pars[which(names(pars) == "m")], function(x) min(1, x))
  
  all_data <- create_data(run, "protocol",
                          paste0("~/net/home/ferret/",parTab_file), 
                          paste0("~/net/home/ferret/",exposureTab_file), 
                          ngroup=5, nstrain=5,nindiv=3,times=times,
                          wd=save_wd1, normal=TRUE,
                          pars=pars)
  
  created_data_runs["data_file"] <- strsplit(all_data$filename,split="~/net/home/ferret/")[[1]][2]
  created_data_runs["parTab_file"] <- strsplit(all_data$parTab_filename,split="~/net/home/ferret/")[[1]][2]
  
  created_data_runs["runID"] <- "protocol"
  all_runs_sim <- rbind(all_runs_sim, created_data_runs)
  
  times <- seq(0,70,by=7)
  all_data <- create_data(run, "weekly",
                          paste0("~/net/home/ferret/",parTab_file), 
                          paste0("~/net/home/ferret/",exposureTab_file), 
                          ngroup=5, nstrain=5,nindiv=3,times=times,
                          wd=save_wd1, normal=TRUE,
                          pars=pars)
  created_data_runs["data_file"] <- strsplit(all_data$filename,split="~/net/home/ferret/")[[1]][2]
  created_data_runs["parTab_file"] <- strsplit(all_data$parTab_filename,split="~/net/home/ferret/")[[1]][2]
  created_data_runs["runID"] <- "weekly"
  all_runs_sim <- rbind(all_runs_sim, created_data_runs)
  times <- seq(0,70,by=3)
  
  all_data <- create_data(run, "3days",
                          paste0("~/net/home/ferret/",parTab_file), 
                          paste0("~/net/home/ferret/",exposureTab_file), 
                          ngroup=5, nstrain=5,nindiv=3,times=times,
                          wd=save_wd1, normal=TRUE,
                          pars=pars)
  
  created_data_runs["data_file"] <- strsplit(all_data$filename,split="~/net/home/ferret/")[[1]][2]
  created_data_runs["parTab_file"] <- strsplit(all_data$parTab_filename,split="~/net/home/ferret/")[[1]][2]
  created_data_runs["runID"] <- "3days"
  all_runs_sim <- rbind(all_runs_sim, created_data_runs)
}

## Create a run tracker for these simulated runs
write.table(all_runs_sim,"~/net/home/ferret/inputs_Jan2019/run_tracker_sim_best_models.csv",sep=",",row.names=FALSE)
write.table(all_runs_sim,"~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker_sim_best_models.csv",sep=",",row.names=FALSE)
