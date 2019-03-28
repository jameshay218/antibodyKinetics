setwd("~/net/home/ferret")
source("scripts_Jan2019/full_run_paralleltemp.R")
source("scripts_Jan2019/cluster_setup_pt.R")


## MCMC control parameters
n_temperatures <- 30
mcmcPars <- list("iterations"=4000000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=200,"adaptive_period"=1000000,"save_block"=500,"temperature" = seq(1,301,length.out=n_temperatures),
                 "parallel_tempering_iter" = 5,"max_adaptive_period" = 1000000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 2000000)


## Run 3 chains per model
maxChain <- 3

## Some universal options
sim <- FALSE
fixed_S <- FALSE
typing <- TRUE

times <- c(0,21,37,49,70)  

## Read in the run tracker table
runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)

#######################################
#######################################
## Enumerate run IDs for each chain
chainNos <- rep(1:maxChain,each=nrow(runs))
runs <- mefa:::rep.data.frame(runs,maxChain)
runs <- cbind(runs,"chainNo"=chainNos)
runs <- runs[runs$runID %in% c(21L, 53L, 62L),]

## Correct classes for runs table
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$runName <- as.character(runs$runName)
runs$runID <- as.character(runs$runID)

####
########
runs$dat_file <- "inputs/real_data_simple.csv"
runs$dat_file <- as.character(runs$dat_file)
runs <- runs[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file")]
runs$parTab_file <- paste0("inputs_Jan2019/free_tp/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")
runs$runID <- paste0(runs$runID, "_free_tp")

## Submit jobs to cluster
jobs_free_tp <- queuer::enqueue_bulk(obj1,runs,"full_run_paralleltemp",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=5,nstrain=5,nindiv=3,
                             times=times,
                             mcmcPars=mcmcPars,
                             sim=sim,
                             do_call=TRUE,timeout=0)
jobs_free_tp$status()
jobs_free_tp$ids[which(jobs_free_tp$status() == "ERROR")]

#i <- 29
#print(runs[i,])
#devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics/")
#full_run_paralleltemp(runs$runName[i],
#                      runs$runID[i],
#                      runs$chainNo[i],
#                      runs$parTab_file[i],
#                      runs$exposureTab_file[i],
#                      runs$dat_file[i],
#                      fixed_S=FALSE,
#                      typing=TRUE,
#                      ngroup=5,
#                      nstrain=5,
 #                     nindiv=3,
 #                     times=times,
 #                     mcmcPars=mcmcPars,
 #                     sim=TRUE)
