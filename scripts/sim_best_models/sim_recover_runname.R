######################
## JAMES HAY 11.02.2019 - jameshay218@gmail.com
## This script simulates antibody trajectories for a given model run idenfifier, and
## then runs the MCMC framework for this simultaed data

library(ggplot2)
setwd("~/net/home/ferret/inputs_Jan2019/")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
runName <- "CYAY6BY"
times <- seq(0,70,by=1)
times <- c(0,21,37,49,70)
parTab <- read.csv("~/net/home/ferret/inputs_Jan2019/sim_best_model/CYAY6BY/protocol_CYAY6BY_parTab.csv", stringsAsFactors=FALSE)
options <- convert_runName_to_options(runName)
parTab <- parTab_modification(parTab, options,FALSE)

parTab[parTab$names == "mod","fixed"] <- 1
parTab[parTab$names == "mod","values"] <- 1
parTab[parTab$names == "sigma" & parTab$type == "infection1","fixed"] <- 1
parTab[parTab$names == "m","upper_bound"] <- 10
#parTab[parTab$names == "m","values"] <- 0
#parTab[parTab$names == "m","fixed"] <- 1
parTab[parTab$names %in% c("sigma","beta"),"upper_bound"] <- 100
parTab[parTab$names == "ts","upper_bound"] <- 30

exposureTab <- read.csv("~/net/home/ferret/inputs_Jan2019/sim_new_version/IYTY6BN/D_exposureTab.csv", stringsAsFactors=FALSE)
all_data <- create_data(runName, 3,
                        "~/net/home/ferret/inputs_Jan2019/sim_new_version/IYTY6BN/parTab.csv", 
                        "~/net/home/ferret/inputs_Jan2019/sim_new_version/IYTY6BN/D_exposureTab.csv", 
                        ngroup=5, nstrain=5,nindiv=3,times=times,
                        wd=save_wd, normal=TRUE)

dat <- as.matrix(read.csv(all_data$filename))
dat <- rbind(times, dat)
mcmcPars <- c("adaptive_period"=20000,"iterations"=100000,"opt_freq"=2000,"thin"=1,"save_block"=1000,"popt"=0.44)
individuals <- rep(3,5)

data(ferret_titres)
ferret_titres <- as.matrix(ferret_titres[,4:ncol(ferret_titres)])
ferret_titres <- rbind(c(0,21,37,49,70),ferret_titres)
rownames(ferret_titres) <- NULL
parTab[parTab$names == "boost_limit","fixed"] <- 0
run_without_point <- antibodyKinetics::run_MCMC(parTab=parTab, data=ferret_titres, mcmcPars=mcmcPars, 
                                                filename="wow3",CREATE_POSTERIOR_FUNC=create_model_group_func_cpp,
                                                mvrPars=NULL,PRIOR_FUNC=NULL, 
                                                version="posterior",form="isolated",
                                                individuals=individuals,exposureTab=exposureTab,
                                                cross_reactivity=TRUE,typing=TRUE)
chain <- read.csv(run_without_point$file)
