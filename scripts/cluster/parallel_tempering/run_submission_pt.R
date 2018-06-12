setwd("~/net/home/ferret")
source("scripts/parallel_tempering/full_run_paralleltemp.R")
source("scripts/parallel_tempering/cluster_setup_pt.R")

## Run 3 chains per model
maxChain <- 3

## Some universal options
sim <- FALSE
fixed_S <- FALSE
typing <- TRUE

times <- c(0,21,36,49,70)   


## Read in the run tracker table
runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs$parTab_file <- paste0("inputs/parTabs/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")


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

####
### Sim or real data?
if(sim){
    runs$dat_file <- paste0("inputs/simulation_data/",runs$runName,"_data.csv")
} else {
    runs$dat_file <- "inputs/real_data_simple.csv"
}
########

runs <- runs[,c("runName","runID","chainNo","parTab_file","exposureTab_file","dat_file")]


## MCMC control parameters
n_temperatures <- 10
mcmcPars <- list("iterations"=5000000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=100,"adaptive_period"=1000000,"save_block"=500,"temperature" = seq(1,101,by=n_temperatures),
                 "parallel_tempering_iter" = 10,"max_adaptive_period" = 1000000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 5000000)


## Submit jobs to cluster
jobs <- queuer::enqueue_bulk(obj1,runs,"full_run_paralleltemp",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=5,nstrain=5,nindiv=3,
                             times=times,
                             mcmcPars=mcmcPars,
                             sim=sim,
                             do_call=TRUE,timeout=0)

