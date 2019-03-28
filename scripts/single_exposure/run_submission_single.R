######################
## JAMES HAY 11.02.2019 - jameshay218@gmail.com
## Script to submit single exposure analysis runs
## to the DIDE cluster. File paths would need to be changed
## for use on a different computer.

setwd("~/net/home/ferret")
source("scripts_Jan2019/full_run_paralleltemp.R")
source("scripts_Jan2019/cluster_setup_pt.R")

## Some universal options
## Run 3 chains per model
maxChain <- 3
## MCMC control parameters
n_temperatures <- 25
mcmcPars <- list("iterations"=2000000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=200,"adaptive_period"=1000000,"save_block"=500,
                 "temperature" = seq(1,251,length.out=n_temperatures),
                 "parallel_tempering_iter" = 5,"max_adaptive_period" = 1000000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 2000000)

fixed_S <- FALSE
typing <- TRUE

times <- c(0,21,37,49,70)  
times_days <- seq(0,72,by=3)
times_weekly <- seq(0,70,by=7)

#######################################
## REAL SINGLE EXPOSURE MAIN MODEL
#######################################
sim <- FALSE
## Read in the run tracker table
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_single.csv",stringsAsFactors=FALSE)

## Enumerate run IDs for each chain
chainNos <- rep(1:maxChain,each=nrow(runs))
runs <- mefa:::rep.data.frame(runs,maxChain)
runs <- cbind(runs,"chainNo"=chainNos)
runs$parTab_file <- paste0("inputs_Jan2019/single_exposure/",runs$parTab_file, ".csv")
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- paste0("inputs_Jan2019/single_exposure/",runs$exposureTab_file, ".csv")
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$runName <- as.character(runs$runName)
runs$runID <- as.character(runs$runID)
runs$dat_file <- "inputs_Jan2019/single_exposure/single_real_data.csv"
runs <- runs[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file")]
## Submit jobs to cluster
jobs_real <- queuer::enqueue_bulk(obj1,runs,"full_run_paralleltemp",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=1,nstrain=1,nindiv=3,
                             times=times,
                             mcmcPars=mcmcPars,
                             sim=sim,
                             do_call=TRUE,timeout=0)


#######################################
## SIM MAIN MODEL
#######################################
sim <- TRUE
## Read in the run tracker table
runs_sim <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_single_sim.csv",stringsAsFactors=FALSE)
## Enumerate run IDs for each chain
chainNos <- rep(1:maxChain,each=nrow(runs_sim))
runs_sim <- mefa:::rep.data.frame(runs_sim,maxChain)
runs_sim <- cbind(runs_sim,"chainNo"=chainNos)
runs_sim$parTab_file <- paste0("inputs_Jan2019/single_exposure/",runs_sim$parTab_file, ".csv")
runs_sim$parTab_file <- as.character(runs_sim$parTab_file)
runs_sim$exposureTab_file <- paste0("inputs_Jan2019/single_exposure/",runs_sim$exposureTab_file, ".csv")
runs_sim$exposureTab_file <- as.character(runs_sim$exposureTab_file)
runs_sim$runName <- as.character(runs_sim$runName)
runs_sim$runID <- as.character(runs_sim$runID)
runs_sim$dat_file <- paste0("inputs_Jan2019/", runs_sim$dat_file)
runs_sim$sim <- as.logical(runs_sim$sim)
runs_sim$run_times <- as.character(runs_sim$run_times)
runs_sim <- runs_sim[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file","sim","run_times")]
## Submit jobs to cluster
jobs_sim_protocol <- queuer::enqueue_bulk(obj1,runs_sim[runs_sim$runID == "sim_protocol_3",],
                                          "full_run_paralleltemp",
                                          fixed_S=fixed_S,
                                          typing=typing,
                                          ngroup=1,nstrain=1,nindiv=3,
                                          times=times,
                                          mcmcPars=mcmcPars,
                                          do_call=TRUE,timeout=0)

jobs_sim_days <- queuer::enqueue_bulk(obj1,runs_sim[runs_sim$runID == "sim_frequent_3",],
                                      "full_run_paralleltemp",
                                      fixed_S=fixed_S,
                                      typing=typing,
                                      ngroup=1,nstrain=1,nindiv=3,
                                      times=times,
                                      mcmcPars=mcmcPars,
                                      do_call=TRUE,timeout=0)

jobs_sim_days_10 <- queuer::enqueue_bulk(obj1,runs_sim[runs_sim$runID == "sim_frequent_10",],
                                         "full_run_paralleltemp",
                                         fixed_S=fixed_S,
                                         typing=typing,
                                         ngroup=1,nstrain=1,nindiv=10,
                                         times=times,
                                         mcmcPars=mcmcPars,
                                         do_call=TRUE,timeout=0)

jobs_sim_weekly <- queuer::enqueue_bulk(obj1,runs_sim[runs_sim$runID == "sim_weekly_3",],
                                        "full_run_paralleltemp",
                                        fixed_S=fixed_S,
                                        typing=typing,
                                        ngroup=1,nstrain=1,nindiv=3,
                                        times=times,
                                        mcmcPars=mcmcPars,
                                        do_call=TRUE,timeout=0)

jobs_sim_versions <-  queuer::enqueue_bulk(obj1,runs_sim[runs_sim$runID %in%
                                                         c("biphasic_fixed","biphasic_full","monophasic_fixed",
                                                           "monophasic_full"),],
                                           "full_run_paralleltemp",
                                           fixed_S=fixed_S,
                                           typing=typing,
                                           ngroup=1,nstrain=1,nindiv=3,
                                           times=times,
                                           mcmcPars=mcmcPars,
                                           do_call=TRUE,timeout=0)

