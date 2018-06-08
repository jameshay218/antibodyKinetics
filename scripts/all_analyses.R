library(coda)
library(antibodyKinetics)
source("~/Documents/Ferret_Model/analysis/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/analysis/model_comparison_functions.R")
top_wd <- getwd()

#############################
## Run model comparison?
do_model_comparison <- TRUE

skips <- read.csv("~/net/home/ferret/rerun_ids.csv")
skips <- skips[skips$chain != "rerun",]
ids <- skips$id
use_chain <- skips$chain

ids <- c(999)
use_chain <- c("wow")
#ids <- c(3,7,17,18,21,22,39,49,50,51,52,53,54,62,74,82,85,86,89,113,114,117,118,119,120,121,122,127,128)
ids <- c(54, 62, 61,117, 53, 119, 126, 22, 30, 125, 118, 113,94)
#use_runs <- read.csv("~/Documents/test_convergence/convergence_check_manual_correction.csv")
#ids <- use_runs[use_runs$manual == "rerun","runID"]
#use_chain <- use_runs[use_runs$manual == "rerun","manual"]

## Check for converence?
check_convergence <- TRUE

## Extract parameter estimates?
do_parameter_estimates <- TRUE

## Input area
res_wd <- "~/Documents/Ferret_Model/effective_pars"
chain_wd_base <- "~/Documents/Ferret_Model/results_112017/outputs/"

dat_file <- "~/net/home/ferret/inputs/real_data_simple.csv"
if(!dir.exists(res_wd)) dir.create(res_wd,recursive=TRUE)

multi_adaptive <- 1000000
univ_adaptive <- 1000000

parTab_loc <- "~/net/home/ferret/inputs/parTabs/"
exposureTab_loc <- "~/net/home/ferret/inputs/exposureTabs/"

runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
#runs <- runs[runs$antigenic_seniority == "Y",]
###
## IF SKIPPING
###
runs <- runs[runs$runID %in% ids,]
times <- c(0,21,36,49,70)
n <- 1000

#############################

ess_vector_multi <- NULL
ess_names_multi <- NULL
gelman.psrf_multi <- NULL
psrf_names_multi <- NULL
gelman.mpsrf_multi <- NULL

ess_vector_univ <- NULL
ess_names_univ <- NULL
gelman.psrf_univ <- NULL
psrf_names_univ <- NULL
gelman.mpsrf_univ <- NULL

ess_check_multi <- NULL
ess_check_univ <- NULL
gelman_check_multi <- NULL
gelman_check_univ <- NULL
univ_chain_fine <- NULL
multi_chain_fine <- NULL
rerun <- NULL

bics <- rep(NA, nrow(runs))
names <- NULL
waics <- rep(NA, nrow(runs))
pwaics <- rep(NA, nrow(runs))
residuals <- NULL
waic_names <- NULL

all_estimates <- NULL


for(i in 1:nrow(runs)){
    runName <- runs$runName[i]
    runID <- runs$runID[i]
    print(paste0(runID, "_", runName))
    parTab_file <- paste0(parTab_loc,runs$parTab_file[i],".csv")
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)

    ## What options were set for this model?
    options <- antibodyKinetics::convert_runName_to_options(runName)
    parTab <- antibodyKinetics::parTab_modification(parTab,options,FALSE)
    parTab[parTab$names == "mod","fixed"][1] <- 1
    setwd(top_wd)
    chain_wd <- paste0(chain_wd_base,"/",runID,"_",runName)
    
    setwd(chain_wd)
    print("Reading in multivariate chains...")
    chains_multi <- load_mcmc_chains(chain_wd,parTab,TRUE,1,multi_adaptive,TRUE)
    
    print("Reading in univariate chains...")
    chains_univ <- load_mcmc_chains(chain_wd,parTab,TRUE,1,univ_adaptive,FALSE)
    print("Done")

    ess_multi_fine <- FALSE
    gelman_multi_fine <- FALSE
    
    if(!is.null(chains_multi)){
        print("Analysing multivariate chains...")
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
        ess_vector_multi[i] <- min(ess_multi)
        ess_names_multi[i] <- names(which.min(ess_multi))
        if(ess_vector_multi[i] > 200) ess_multi_fine <- TRUE
        if(!is.na(gelman_multi)){
            print(paste0("Gelman PSRF results, name:",names(which.max(gelman_multi$psrf[,2])),
                         ";  max: ", max(gelman_multi$psrf[,2])))
            print(paste0("Gelman MPSRF: ", gelman_multi$mpsrf))
            
            
            gelman.psrf_multi[i] <- max(gelman_multi$psrf[,2])
            psrf_names_multi[i] <- names(which.max(gelman_multi$psrf[,2]))
            gelman.mpsrf_multi[i] <- gelman_multi$mpsrf

            if(gelman.psrf_multi[i] < 1.15 && gelman.mpsrf_multi < 1.15) gelman_multi_fine <- TRUE

            ### IF FORCING ACCEPTANCE
            if(runID %in% ids & use_chain == "multi") gelman_multi_fine <- TRUE
            
        } else {
            gelman.psrf_multi[i] <- NA
            psrf_names_multi[i] <- NA
            gelman.mpsrf_multi[i] <- NA
        }
    } else {
        print("Multivariate chain problem")
    }

    ess_univ_fine <- FALSE
    gelman_univ_fine <- FALSE
    
    if(!is.null(chains_univ)){
        print("Analysing univariate chains...")
        gelman_univ <- tryCatch({
            tmp <- gelman.diag(chains_univ[["list"]])
        }, warning = function(w){
            print(w)
        }, error = function(e){
            print(e)
            tmp <- NA
        }, finally = {
            tmp
        })
        
        ess_univ <- effectiveSize(chains_univ[["list"]])

        print("Univariate chain convergence:")
        print(paste0("ESS results, name: ", names(which.min(ess_univ)),"; value: ",min(ess_univ)))
        ess_vector_univ[i] <- min(ess_univ)
        ess_names_univ[i] <- names(which.min(ess_univ))
        if(ess_vector_univ[i] > 200) ess_univ_fine <- TRUE
        
        ## Plot the combined MCMC chains
        pdf(paste0(res_wd,"/",runID,"_",runName,"_univariate_chains.pdf"))
        plot(chains_univ[["list"]])
        dev.off()
        
        if(!is.na(gelman_univ)){
            print(paste0("Gelman PSRF results, name:",names(which.max(gelman_univ$psrf[,2])),
                         ";  max: ", max(gelman_univ$psrf[,2])))
            print(paste0("Gelman MPSRF: ", gelman_univ$mpsrf))
            gelman.psrf_univ[i] <- max(gelman_univ$psrf[,2])
            psrf_names_univ[i] <- names(which.max(gelman_univ$psrf[,2]))
            gelman.mpsrf_univ[i] <- gelman_univ$mpsrf
            if(gelman.psrf_univ[i] < 1.15 && gelman.mpsrf_univ < 1.15) gelman_univ_fine <- TRUE

            ## IF FORCING ACCEPTANCE
            if(runID %in% ids & use_chain == "univ") gelman_univ_fine <- TRUE
        } else {
            gelman.psrf_univ[i] <- NA
            psrf_names_univ[i] <- NA
            gelman.mpsrf_univ[i] <- NA
        }
    } else {
        print("Univariate chain problem")
    }
    print("")
    ess_check_multi[i] <- ess_multi_fine
    ess_check_univ[i] <- ess_univ_fine
    gelman_check_multi[i] <- gelman_multi_fine
    gelman_check_univ[i] <- gelman_univ_fine

    univ_chain_fine[i] <- ess_univ_fine & gelman_univ_fine
    multi_chain_fine[i] <- ess_multi_fine & gelman_multi_fine

    
    rerun[i] <- TRUE
    if(univ_chain_fine[i] | multi_chain_fine[i]) rerun[i] <- FALSE
    if(runID %in% ids) rerun[i] <- TRUE

##################
    ## Extract parameter values
##################
    if(do_parameter_estimates){
        chain <- chains_multi[["chain"]]
        if(univ_chain_fine) chain <- chains_univ[["chain"]]

        names <- c(parTab[which(parTab$fixed == 0), "names"],"lnlike")
        names[names == "mod"] <- paste0("mod",1:4)
        types <- c(parTab[which(parTab$fixed == 0), "type"],"all")
        tmpChain <- chain
        ## Getting post initial waning mu value
        if(!options$monophasic_waning){
            mu_chain <- chain[,colnames(chain) %in% c("mu",paste0("mu.",1:5))]
            dp_chain <- chain[,colnames(chain) %in% c("dp",paste0("dp.",1:5))]
            for(j in 1:ncol(mu_chain)) mu_chain[,j] <- mu_chain[,j] *(1- dp_chain[,j])
            colnames(mu_chain) <- paste0("adj",colnames(mu_chain))
            tmpChain <- cbind(chain, mu_chain)
            names <- c(names, colnames(mu_chain))
            types <- c(types, parTab[which(parTab$names == "mu"),"type"])
        }


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
        all_estimates <- rbind(all_estimates,par_estimates)
    }
    
    if(do_model_comparison){
        print("Running model comparison analysis...")
        if(multi_chain_fine[i]){
            adaptive_period <- multi_adaptive
            use_multi_chain <- TRUE
        } else if(univ_chain_fine[i]){
            adaptive_period <- univ_adaptive
            use_multi_chain <- FALSE
        } else {
            adaptive_period <- univ_adaptive
            chain <- NULL
            use_multi_chain <- FALSE
        }
        exposureTab <- read.csv(paste0(exposureTab_loc,runs$exposureTab_file[i],".csv"),stringsAsFactors=FALSE)
        model_comparison_res <- model_comparison_analyses(chain_wd, use_multi_chain,
                                                          adaptive_period,
                                                          parTab, exposureTab,
                                                          dat_file, options, TRUE,
                                                          times, n)        
        bics[i] <- model_comparison_res[["BIC"]]
        waics[i] <- model_comparison_res[["WAIC"]]
        pwaics[i] <- model_comparison_res[["pwaic"]]
        waic_names[i] <- paste0(runID,"_",runName)
        tmp_residuals <- as.data.frame(model_comparison_res[["mle_res"]])
        tmp_residuals$runName <- runName
        residuals <- rbind(residuals,tmp_residuals)
        print("...model comparison analysis done")
    }
}

results <- data.frame(runID=runs$runID, runName=runs$runName,
                      ess_names_multi,ess_vector_multi,psrf_names_multi,gelman.psrf_multi,gelman.mpsrf_multi,
                      ess_names_univ,ess_vector_univ,psrf_names_univ,gelman.psrf_univ,gelman.mpsrf_univ,
                      ess_check_multi,gelman_check_multi,ess_check_univ,gelman_check_univ,
                      univ_chain_fine,multi_chain_fine, rerun
                      )

model_comparison <- data.frame(bics, waics)
if(do_model_comparison) results <- cbind(results, model_comparison)

setwd(res_wd)
write.table(residuals,"mle_residuals.csv",row.names=FALSE,sep=",")
write.table(results,"convergence_check.csv",sep=",",row.names=FALSE)

if(do_parameter_estimates) write.table(all_estimates,"parameter_estimates.csv",row.names=FALSE,sep=",")

library(ggplot2)
pdf("residuals_final.pdf")
p <- ggplot(residuals) + geom_point(aes(x=y,y=residual)) + facet_wrap(~runName)
print(p)
dev.off()


ggplot(all_estimates[all_estimates$par_name == "mu" & all_estimates$runName %in% tmp$runName,]) + geom_pointrange(aes(x=runName,y=median,ymax=upper_quantile,ymin=lower_quantile)) + facet_wrap(~type)
