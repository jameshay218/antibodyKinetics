setwd("~/net/home/ferret")
source("scripts/full_run.R")
source("scripts/cluster_setup.R")

## Run 3 chains per model
maxChain <- 3

## Some universal options
sim <- FALSE
fixed_S <- FALSE
typing <- TRUE

times <- c(0,21,36,49,70)   

## Read in the run tracker table
runs <- read.csv("inputs/run_tracker.csv",stringsAsFactors=FALSE)
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
mcmcPars1 <- c("adaptive_period"=1000000,"iterations"=2000000,"opt_freq"=1000,"thin"=1000,"save_block"=1000,"popt"=0.1)
mcmcPars2 <- c("adaptive_period"=1000000,"iterations"=2000000,"opt_freq"=1000,"thin"=1000,"save_block"=1000,"popt"=0.1)

jobs <- queuer::enqueue_bulk(obj1,runs,"full_run",
                             fixed_S=fixed_S,
                             typing=typing,
                             ngroup=5,nstrain=5,nindiv=3,
                             times=times,
                             mcmcPars1=mcmcPars1,mcmcPars2=mcmcPars2,
                             sim=sim,
                             do_call=TRUE,timeout=0)
