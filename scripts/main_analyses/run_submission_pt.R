setwd("~/net/home/ferret")
source("scripts_Jan2019/full_run_paralleltemp.R")
source("scripts_Jan2019/cluster_setup_pt.R")

## Run 3 chains per model
maxChain <- 3

## MCMC control parameters
n_temperatures <- 40
mcmcPars <- list("iterations"=4000000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=200,"adaptive_period"=1000000,"save_block"=500,"temperature" = seq(1,251,length.out=n_temperatures),
                 "parallel_tempering_iter" = 5,"max_adaptive_period" = 1000000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 2000000)

## Some universal options
sim <- FALSE
fixed_S <- FALSE
typing <- TRUE
times <- c(0,21,37,49,70)  

#######################################
## REAL DATA FITTING
#######################################
## Read in the run tracker table
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_all.csv",stringsAsFactors=FALSE)
runs$parTab_file <- paste0("inputs/parTabs/",runs$parTab_file, ".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")
## Enumerate run IDs for each chain
chainNos <- rep(1:maxChain,each=nrow(runs))
runs <- mefa:::rep.data.frame(runs,maxChain)
runs <- cbind(runs,"chainNo"=chainNos)
## Correct classes for runs table
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$runName <- as.character(runs$runName)
runs$runID <- as.character(runs$runID)
runs$dat_file <- "inputs/real_data_simple.csv"
runs$dat_file <- as.character(runs$dat_file)
runs <- runs[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file")]
## Submit jobs to cluster
jobs <- queuer::enqueue_bulk(obj1,runs,"full_run_paralleltemp",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=5,nstrain=5,nindiv=3,
                             times=times,
                             mcmcPars=mcmcPars,
                             sim=sim,
                             do_call=TRUE,timeout=0)


#######################################
## SIM DATA FITTING
#######################################
sim <- TRUE
runs_sim <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_all_sim.csv",stringsAsFactors=FALSE)
runs_sim$parTab_file <- paste0("inputs/parTabs/",runs_sim$parTab_file, ".csv")
runs_sim$exposureTab_file <- paste0("inputs/exposureTabs/",runs_sim$exposureTab_file,".csv")
## Enumerate run IDs for each chain
chainNos <- rep(1:maxChain,each=nrow(runs_sim))
runs_sim <- mefa:::rep.data.frame(runs_sim,maxChain)
runs_sim <- cbind(runs_sim,"chainNo"=chainNos)
## Correct classes for runs table
runs_sim$parTab_file <- as.character(runs_sim$parTab_file)
runs_sim$exposureTab_file <- as.character(runs_sim$exposureTab_file)
runs_sim$runName <- as.character(runs_sim$runName)
runs_sim$runID <- as.character(runs_sim$runID)
runs_sim$dat_file <- as.character(runs_sim$dat_file)
runs_sim <- runs_sim[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file")]
## Submit jobs to cluster
jobs_sim <- queuer::enqueue_bulk(obj1,runs_sim,"full_run_paralleltemp",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=5,nstrain=5,nindiv=3,
                             times=times,
                             mcmcPars=mcmcPars,
                             sim=sim,
                             do_call=TRUE,timeout=0)
