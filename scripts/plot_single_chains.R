## Plots the univariate chains
source("~/Documents/Ferret_Model/analysis/convergence_check_funcs.R")
library(coda)

runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runID <- 53
runName <- runs$runName[which(runs$runID == runID)]
wd <- paste0(runID,"_",runName)
wd <- paste0("~/net/home/ferret/outputs_real/",wd)

univ_adaptive <- 1000000

parTab_loc <- "~/net/home/ferret/inputs/parTabs/"
parTab_file <- runs$parTab_file[which(runs$runID == runID)]
parTab_file <- paste0(parTab_loc,parTab_file,".csv")
parTab <- read.csv(parTab_file,stringsAsFactors = FALSE)
chains <- load_mcmc_chains(wd,parTab,FALSE,1,univ_adaptive,FALSE,FALSE)

chains <- chains[["list"]]
for(i in 1:length(chains)){
  filename <- paste0(wd,"/univariate_chain_plot_",i,".pdf")
  pdf(filename)
  plot(as.mcmc(chains[[i]]))
  dev.off()
}


#tmp_indices <- read.csv("~/Documents/Ferret_Model/tmp_transfer_indices.csv")
tmp_indices <- tmp[tmp$`Antigenic Seniority`=="Absent",c("runID","runName")]

wd_new <- "~/Documents/Ferret_Model/raw_results_22092017/outputs_real/"
wd_old <- "~/Documents/Ferret_Model/raw_results_08092017/outputs_real/"
tmp_wd <- "~/Documents/Ferret_Model/raw_results_08092017/tmp/"

tmp_indices$name <- paste0(tmp_indices$runID,"_",tmp_indices$runName)
file.rename(paste0(wd_old,tmp_indices$name),paste0(wd_new,tmp_indices$name))

unlink(paste0(wd_old,tmp_indices$name))
file.rename(paste0(wd_new,tmp_indices$name),paste0(wd_old,tmp_indices$name))

