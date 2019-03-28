######################
## JAMES HAY 11.02.2019 - jameshay218@gmail.com
## Script to assess convergence and calculate parameter summaries
## for simulated single exposure analysis


## Load required packages and auxiliary functions
library(coda)
library(loo)
#library(antibodyKinetics)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")

top_wd <- getwd()

#############################
## Run model comparison?
do_model_comparison <- TRUE
## Check for converence?
check_convergence <- TRUE
## Extract parameter estimates?
do_parameter_estimates <- TRUE


## Input area
## Full file path to giving location to save all results
res_wd <- "~/Drive/Influenza/Ferret/PLOS Comp Bio/sim_single_exposure_analysis"
## Full file path to folder containing the MCMC chains
#chain_wd_base <- "~/net/home/ferret/outputs/"
chain_wd_base <- "/media/james/Storage 2/ferret_results/sim_single_exposure"
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_single_sim.csv",stringsAsFactors=FALSE)

## Location of the data file used
if(!dir.exists(res_wd)) dir.create(res_wd,recursive=TRUE)

## How many initial iterations of the MCMC chains should be disgarded
## as burn in/adaptive period?
adaptive <- adaptive_period <- 1000000

## Where are the input parameter and exposure tables saved?
#runs$parTab_file <- paste0("~/net/home/ferret/",runs$parTab_file)
runs$exposureTab_file <- paste0("~/net/home/ferret/inputs_Jan2019/single_exposure/",runs$exposureTab_file,".csv")
runs$parTab_file <- paste0("~/net/home/ferret/inputs_Jan2019/single_exposure/",runs$parTab_file,".csv")
runs$dat_file <- paste0("~/net/home/ferret/inputs_Jan2019/",runs$dat_file)

times_protocol <- c(0,21,37,49,70)
times_weekly <- seq(0,70,by=7)
times_days <- seq(0,72,by=3)
n <- 1000

#############################
ess_vector <- NULL
ess_names <- NULL
gelman.psrf <- NULL
psrf_names <- NULL
gelman.mpsrf <- NULL

ess_check <- NULL
gelman_check <- NULL
chain_fine <- NULL
rerun <- NULL

bics <- rep(NA, nrow(runs))
names <- NULL
waics <- rep(NA, nrow(runs))
elpd_loos <- rep(NA, nrow(runs))
elpd_loos_se<- rep(NA, nrow(runs))
p_loos <- rep(NA, nrow(runs))
p_loos_se <- rep(NA, nrow(runs))
looics <- rep(NA, nrow(runs))
looics_se <- rep(NA, nrow(runs))

pwaics <- rep(NA, nrow(runs))
residuals <- NULL
waic_names <- NULL

all_estimates <- NULL
all_ess <- NULL
all_gelman <- NULL
all_loo_estimates <- NULL
all_loo_labels <- NULL
all_k_percentages <- NULL
## For each run/model
for(i in 1:nrow(runs)){
    runName <- runs$runName[i]
    runID <- runs$runID[i]
    print(paste0(runID, "_", runName))
    parTab_file <- runs$parTab_file[i]
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
    parTab[parTab$names == "m","upper_bound"] <- 12
    dat_file <- runs$dat_file[i]
    run_times <- runs$run_times[i]
    nindiv <- 3
    
    if(runID == "sim_frequent_10") nindiv <- 10
    if(run_times == "3days") times <- times_days
    if(run_times == "protocol") times <- times_protocol
    if(run_times == "weekly") times <- times_weekly
    
    ## What options were set for this model?
    options <- antibodyKinetics::convert_runName_to_options(runName)
    parTab <- antibodyKinetics::parTab_modification(parTab,options,FALSE)
    setwd(top_wd)
    chain_wd <- paste0(chain_wd_base,"/",runID,"_",runName)
    chains_multi <- load_mcmc_chains(chain_wd,parTab,TRUE,1,adaptive,TRUE,FALSE,TRUE)
    
    gelman_multi <- tryCatch({
        tmp <- gelman.diag(chains_multi[["list"]])
    }, warning = function(w){
        print(w)
    }, error = function(e){
        print(e)
        tmp <- NA
    }, finally = {
        tmp
    })
    ess_multi <- effectiveSize(chains_multi[["list"]])
    
    print("Multivariate chain convergence:")
    print(paste0("ESS results, name: ", names(which.min(ess_multi)),"; value: ",min(ess_multi)))
    
    ## Plot the combined MCMC chains
    pdf(paste0(res_wd,"/",runID,"_",runName,"_multivariate_chains.pdf"))
    plot(chains_multi[["list"]])
    dev.off()
    
    ## Save the minimum ESS, worst gelman diag and multivariate gelman. Record
    ## which parameter this relates to
    ess_vector[i] <- min(ess_multi)
    ess_names[i] <- names(which.min(ess_multi))
    gelman_fine <- FALSE
    if(ess_vector[i] > 200) ess_fine <- TRUE
    if(!any(is.na(gelman_multi))){
        print(paste0("Gelman PSRF results, name:",names(which.max(gelman_multi$psrf[,2])),
                     ";  max: ", max(gelman_multi$psrf[,2])))
        print(paste0("Gelman MPSRF: ", gelman_multi$mpsrf))
        
        
        gelman.psrf[i] <- max(gelman_multi$psrf[,2])
        psrf_names[i] <- names(which.max(gelman_multi$psrf[,2]))
        gelman.mpsrf[i] <- gelman_multi$mpsrf
        ## Threshold for convergence between multiple chains
        if(gelman.psrf[i] < 1.1 && gelman.mpsrf < 1.1) gelman_fine <- TRUE
    } else {
        gelman.psrf_multi[i] <- NA
        psrf_names_multi[i] <- NA
        gelman.mpsrf_multi[i] <- NA
    }
    ess_check <- ess_fine
    gelman_check[i] <- gelman_fine
    chain_fine <- ess_fine & gelman_fine
    rerun[i] <- !chain_fine

    ## Parameter extraction
    if(do_parameter_estimates){
        chain <- chains_multi[["chain"]]
        names <- c(parTab[which(parTab$fixed == 0), "names"],"lnlike")
        names[names == "mod"] <- paste0("mod",1:4)
        types <- c(parTab[which(parTab$fixed == 0), "type"],"all")
        names1 <- names
        types1 <- types
        true_vals <- parTab[parTab$fixed == 0,c("names","type","values")]
        if("mod" %in% true_vals$names) true_vals[true_vals$names == "mod","names"] <- paste0("mod",1:4)
        
        tmpChain <- chain
      
        means <- apply(tmpChain, 2, mean)
        medians <- apply(tmpChain, 2, median)
        modes <- apply(tmpChain, 2, estimate_mode)
        quantile_lower <- apply(tmpChain, 2, quantile, probs=0.025)
        quantile_upper <- apply(tmpChain, 2, quantile, probs=0.975)

        par_estimates <- data.frame(runID=runID,runName=runName,
                                    par_name=names,type=types,
                                    mean=means,median=medians,mode=modes,
                                    lower_quantile=quantile_lower,upper_quantile=quantile_upper,
                                    row.names=NULL)
        colnames(true_vals)[1] <- "par_name"
        par_estimates <- merge(par_estimates, true_vals, id.vars=c("par_name","type"),all=TRUE)
        par_estimates$estimated <- ifelse((par_estimates$values > par_estimates$lower_quantile) & 
                                            (par_estimates$values < par_estimates$upper_quantile),
                                          TRUE, FALSE)
        all_estimates <- rbind(all_estimates,par_estimates)
        tmp_ess_dat <- data.frame("runID"=runID,"runName"=runName,"par_name"=names, "type"=types,"ess"=ess_multi)
        tmp_gelman_dat <- data.frame("runID"=runID,"runName"=runName,"par_name"=names,"type"=types,
                                     "rhat_point"=gelman_multi$psrf[,1],"rhat_upper_CI"=gelman_multi$psrf[,2] )
        all_ess <- rbind(all_ess, tmp_ess_dat)
        all_gelman <- rbind(all_gelman, tmp_gelman_dat)
    }
    
    if(do_model_comparison){
        print("Running model comparison analysis...")
    
        exposureTab <- read.csv(runs$exposureTab_file[i],stringsAsFactors=FALSE)
        model_comparison_res <- model_comparison_analyses(chain_wd, FALSE,
                                                          adaptive_period,
                                                          parTab, exposureTab,
                                                          dat_file, options, TRUE,
                                                          times, n, TRUE,
                                                          individuals=c(nindiv),
                                                          n1=10000)        
        bics[i] <- model_comparison_res[["BIC"]]
        waics[i] <- model_comparison_res[["WAIC"]]
        pwaics[i] <- model_comparison_res[["pwaic"]]

        elpd_loos[i] <- model_comparison_res[["elpd_loo"]]
        elpd_loos_se[i] <- model_comparison_res[["elpd_loo_se"]]
        p_loos[i] <- model_comparison_res[["p_loo"]]
        p_loos_se[i] <- model_comparison_res[["p_loo_se"]]
        looics[i] <- model_comparison_res[["looic"]]
        looics_se[i] <- model_comparison_res[["looic_se"]]

        all_loo_estimates[[i]] <- model_comparison_res[["loo_estimate"]]
        x <- model_comparison_res[["loo_estimate"]]$diagnostics$pareto_k
        all_k_percentages[i] <- length(x[x > 0.7])/length(x) * 100
        all_loo_labels[[i]] <- model_comparison_res[["pareto_k"]]
        
        print(waics[i])
        waic_names[i] <- paste0(runID,"_",runName)
        tmp_residuals <- as.data.frame(model_comparison_res[["mle_res"]])
        tmp_residuals$runName <- runName
        residuals <- rbind(residuals,tmp_residuals)
        print("...model comparison analysis done")
    }
}

#runs <- runs[runs$form == "C",]
results <- data.frame(runID=runs$runID, runName=runs$runName,
                      ess_names,ess_vector,psrf_names,gelman.psrf,gelman.mpsrf,
                      ess_check,gelman_check,
                      chain_fine, rerun
                      )
beepr::beep(4)
model_comparison <- data.frame(bics, waics, elpd_loos, elpd_loos_se, p_loos, p_loos_se, looics, looics_se,all_k_percentages)
if(do_model_comparison) results <- cbind(results, model_comparison)

#loo_comparison <- as.data.frame(loo::compare(x=all_loo_estimates))
#loo_comparison$model <- row.names(loo_comparison)
model_names <- paste0("model",1:length(all_loo_estimates))
model_names <- data.frame("model"=model_names,runName=runs$runName)

setwd(res_wd)
save(all_loo_estimates,file="all_loo_estimates.RData")
save(all_loo_labels, file="all_loo_ks.RData")

write.table(residuals,"mle_residuals.csv",row.names=FALSE,sep=",")
write.table(results,"convergence_check.csv",sep=",",row.names=FALSE)


all_estimates <- merge(all_estimates, all_ess, id.vars=c("runID","runName","par_name","type"))
all_estimates <- merge(all_estimates, all_gelman, id.vars=c("runID","runName","par_name","type"))

if(do_parameter_estimates) write.table(all_estimates,"parameter_estimates.csv",row.names=FALSE,sep=",")

