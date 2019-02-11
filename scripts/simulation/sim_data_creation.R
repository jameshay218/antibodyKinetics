library(ggplot2)
library(coda)
setwd("~/net/home/ferret")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")

runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs$parTab_file <- paste0("inputs/parTabs/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")

times <- c(0,21,37,49,70)
all_data <- NULL
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
#runs <- runs[runs$form == "I",]
#runs <- runs[runs$runID %in% results[results$elpd_loos > (max(results$elpd_loos) - 10),"runID"],]
adaptive <- 1000000
for(i in 1:nrow(runs)){
  print(runs$runID[i])
  if(runs$form[i] == "I") save_wd <- "~/net/home/ferret/inputs_Jan2019/sim_all/isolated"
  if(runs$form[i] == "C") save_wd <- "~/net/home/ferret/inputs_Jan2019/sim_all/competitive"
  parTab <- read.csv(runs$parTab_file[i],stringsAsFactors=FALSE)
  chain <- quiet(as.data.frame(load_mcmc_chains(paste0("/media/james/Storage 2/ferret_results/rerun_Jan2019/outputs_real/",
                                                 runs$runID,"_",runs$runName[i]), 
                                          parTab, FALSE, 1, adaptive, FALSE, TRUE,TRUE)[["chain"]]))
  pars <- get_best_pars(chain)
  pars <- pars[1:(length(pars)-1)]
  
  #print(pars)
  pars[which(names(pars) == "m")] <- sapply(pars[which(names(pars) == "m")], function(x) min(1, x))
  all_data[[i]] <- create_data(runs$runName[i], runs$runID[i],
                                 runs$parTab_file[i], runs$exposureTab_file[i], 
                                 ngroup=5, nstrain=5,nindiv=3,times=times,
                                wd=save_wd, normal=TRUE,pars=pars)
  runs[i,"parTab_file"] <- all_data[[i]]$parTab_filename
}
runs$dat_file <- paste0("inputs_Jan2019/sim_all/isolated/",runs$runID, "_",runs$runName,"_data.csv")
write.table(runs, "~/net/home/ferret/inputs_Jan2019/run_tracker_sim_best_models.csv",sep=",",row.names=FALSE)


