setwd("~/net/home/ferret")
source("scripts_Jan2019/full_run_paralleltemp.R")
source("scripts_Jan2019/cluster_setup_pt.R")

## Run 3 chains per model
maxChain <- 5

## Some universal options
sim <- TRUE
fixed_S <- FALSE
typing <- TRUE

times <- c(0,21,37,49,70)  
times1 <-seq(0,70,by=7)
times2 <- seq(0,72,by=3)

## Read in the run tracker table
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_sim_best_models.csv",stringsAsFactors=FALSE)


#######################################
#######################################
## Enumerate run IDs for each chain
chainNos <- rep(1:maxChain,each=nrow(runs))
runs <- mefa:::rep.data.frame(runs,maxChain)
runs <- cbind(runs,"chainNo"=chainNos)

## Correct classes for runs table
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$runName <- as.character(runs$runName)
runs$runID <- as.character(runs$runID)

runs$dat_file <- as.character(runs$dat_file)
runs <- runs[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file")]

## MCMC control parameters
n_temperatures <- 25
mcmcPars <- list("iterations"=2000000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=100,"adaptive_period"=1000000,"save_block"=500,"temperature" = seq(1,101,by=n_temperatures),
                 "parallel_tempering_iter" = 5,"max_adaptive_period" = 1000000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 2000000)


## Submit jobs to cluster
jobs_protocol <- queuer::enqueue_bulk(obj1,runs[runs$runID == "protocol",],"full_run_paralleltemp",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=5,nstrain=5,nindiv=3,
                             times=times,
                             mcmcPars=mcmcPars,
                             sim=sim,
                             do_call=TRUE,timeout=0)

jobs_weekly <- queuer::enqueue_bulk(obj1,
                                    runs[runs$runID == "weekly" &
                                         runs$runName %in% c("CNAY6BY","CYAY6BY","CYTY6BN"),],
                                    "full_run_paralleltemp",
                                    fixed_S=fixed_S,
                                    typing=typing,
                                    ngroup=5,nstrain=5,nindiv=3,
                                    times=times1,
                                    mcmcPars=mcmcPars,
                                    sim=sim,
                                    do_call=TRUE,timeout=0)

jobs_3days <- queuer::enqueue_bulk(obj1,runs[runs$runID == "3days" &
                                             runs$runName %in% c("CNAY6BY","CYAY6BY","CYTY6BN"),],
                                   "full_run_paralleltemp",
                                   fixed_S=fixed_S,
                                   typing=typing,
                                   ngroup=5,nstrain=5,nindiv=3,
                                   times=times2,
                                   mcmcPars=mcmcPars,
                                   sim=sim,
                                   do_call=TRUE,timeout=0)

