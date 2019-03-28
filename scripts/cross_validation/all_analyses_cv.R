## Load required packages and auxiliary functions
library(coda)
library(loo)
#library(antibodyKinetics)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")

## Read in PSIS-LOO estimates from main runs
load("~/Drive/Influenza/Ferret/PLOS Comp Bio/main_results/all_loo_estimates.RData")
load("~/Drive/Influenza/Ferret/PLOS Comp Bio/main_results/all_loo_ks.RData")

old_loo_estimates <- all_loo_estimates[1:64]
tmp_loo_labels <- all_loo_labels[1:64]

## Input area
## Full file path to giving location to save all results
res_wd <- top_wd <- "~/Drive/Influenza/Ferret/PLOS Comp Bio/cross_validation_results"
## Full file path to folder containing the MCMC chains
chain_wd_base <- "/media/james/Storage 2/ferret_results/cross_validation/outputs_reruns/"

## Location of the data file used
if(!dir.exists(res_wd)) dir.create(res_wd,recursive=TRUE)

## How many initial iterations of the MCMC chains should be disgarded
## as burn in/adaptive period?
adaptive <- adaptive_period <- 1000000

## Where are the input parameter and exposure tables saved?
parTab_loc <- "~/net/home/ferret/inputs/parTabs/"
exposureTab_loc <- "~/net/home/ferret/inputs/exposureTabs/"

## Table describing which files and options are used for which runs
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_cv.csv",stringsAsFactors=FALSE)
all_runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_all.csv",stringsAsFactors = FALSE)
all_runs <- all_runs[all_runs$form=="C",]
## For each LOO run, which PSIS-LOO estimate in the list does this correspond to?
loo_indices <- match(runs$runName, all_runs$runName)
###
times <- c(0,21,37,49,70)
## How many posterior samples to take to compute ELPD
n1 <- 1000
thin <- 1

attach(ferret_titres)
ferret_titres_no_indiv <- unique(ferret_titres[,c("group","strain")])
tmp_labels <- tmp_loo_labels[[1]][,c("groups","strains","indiv","all_times")]

## For each LOO run
for(i in 1:nrow(runs)){
    runName <- runs$runName[i]
    runID <- runs$runID[i]
    
    ## Which data point was left out?
    left_strain <- runs$cv_strain[i]
    left_group <- runs$cv_group[i]
    left_indiv <- runs$cv_indiv[i]
    left_time <- runs$cv_time[i]
    
    ## Get old estimate
    tmp_loo_estimates <- old_loo_estimates[[loo_indices[i]]]
    
    print(paste0(runID, "_", runName))
    parTab_file <- paste0(parTab_loc,runs$parTab_file[i],".csv")
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
    parTab[parTab$names == "m","upper_bound"] <- 12
    
    ## What options were set for this model?
    options <- antibodyKinetics::convert_runName_to_options(runName)
    exposureTab <- read.csv(paste0(exposureTab_loc,runs$exposureTab_file[i],".csv"),stringsAsFactors=FALSE)
    parTab <- antibodyKinetics::parTab_modification(parTab,options,FALSE)
    setwd(top_wd)
    chain_wd <- paste0(chain_wd_base,"/",runID,"_",runName)
    
    ## Get ferret titre data without labels (ie. as in fitted data)
    tmp_dat <- ferret_titres
    tmp_dat <- as.matrix(tmp_dat[,4:ncol(tmp_dat)])
    dat <- rbind(times, tmp_dat)
    
    ## Load in MCMC chains
    chain <- as.data.frame(load_mcmc_chains(chain_wd,parTab,FALSE,
                                            thin,adaptive_period,multi_chain,FALSE,PTchain)[["chain"]])
    
    ## Take n1 samples from the chain
    samps <- sample(1:nrow(chain),size = n1,replace=FALSE)
    chain1 <- chain[samps,]
    
    ## Find the left out data point
    cv_data_point <- ferret_titres[ferret_titres$group == left_group &
              ferret_titres$strain == left_strain &
              ferret_titres$indiv == left_indiv,
            left_time + 3]
    
    ## Find the index of the left out data point in the main data matrix
    cv_row <- which(ferret_titres$group == left_group &
      ferret_titres$strain == left_strain &
      ferret_titres$indiv == left_indiv)
    cv_col <- left_time
    
    ## And the index in the model output matrix (not enumerated by individual)
    cv_row_model <- which(ferret_titres_no_indiv$group == left_group &
                            ferret_titres_no_indiv$strain == left_strain)
    
    ## Create the model solving function
    f <- create_model_group_func_cpp(parTab,exposureTab,dat=dat,PRIOR_FUNC=NULL,version="model",
                                       form=options$form,typing = TRUE,cross_reactivity = options$cr,
                                       individuals=c(3,3,3,3,3))

    ## For each sample from the MCMC chain, calculate the posterior probability of observing
    ## the left out point
    tmp_liks <- numeric(nrow(chain1))
    ## for n1 samples
    for(ii in 1:nrow(chain1)){
      pars <- get_index_pars(chain1,ii)
      y <- f(pars, times)
      tmp_liks[ii] <- log(norm_error(y[cv_row_model,cv_col], cv_data_point,pars["S"],pars["MAX_TITRE"]))
    }
    tmp_liks[!is.finite(tmp_liks)] <- -1000000
    
    ## Replace the ELPD for this data point in the stored ELPD estimates
    index_in_loo <- which(tmp_labels$strains == left_strain & 
                            tmp_labels$groups == left_group & 
                            tmp_labels$indiv == left_indiv & 
                            tmp_labels$all_times == times[left_time])
    print(paste0("Pareto k being replaced: ", tmp_loo_estimates$diagnostics$pareto_k[index_in_loo]))
    if(tmp_loo_estimates$diagnostics$pareto_k[index_in_loo] < 0.7){
      print("ERROR PARETO K ESTIMATE FINE")
      stop()
    }
    print(paste0("Updating ELPD estimate of ", tmp_loo_estimates$pointwise[index_in_loo,1], " to ", mean(tmp_liks)))
    tmp_loo_estimates$pointwise[index_in_loo,1] <- mean(tmp_liks)
    old_loo_estimates[[loo_indices[i]]] <- tmp_loo_estimates
    print("...model comparison analysis done")
}

beepr::beep(4)

loo_comparison <- as.data.frame(loo::compare(x=old_loo_estimates))
lpd1 <- lapply(old_loo_estimates,function(x) x$pointwise[,1])
lpd1 <- do.call("cbind",lpd1)
model_weights <- pseudobma_weights(lpd1)
loo_comparison$model <- row.names(loo_comparison)
model_names <- paste0("model",1:length(old_loo_estimates))
model_names <- data.frame("model"=model_names,runName=all_runs$runName)
loo_comparison <- merge(loo_comparison, model_names)
loo_comparison <- loo_comparison[order(loo_comparison$elpd_diff),]

setwd(res_wd)

save(old_loo_estimates,file="all_loo_estimates_cv.RData")

elpds <- unlist(lapply(old_loo_estimates, function(x) sum(x$pointwise[,1])))
elpd_se <- unlist(lapply(old_loo_estimates, function(x) sqrt(var(x$pointwise[,"elpd_loo"])*375)))
results <- data.frame(runID=all_runs$runID,runName=all_runs$runName,cv_elpd=elpds,cv_se=elpd_se)

write.table(results,"convergence_check_cv.csv",sep=",",row.names=FALSE)
write.table(loo_comparison,"loo_comparison.csv",sep=",",row.names=FALSE)
