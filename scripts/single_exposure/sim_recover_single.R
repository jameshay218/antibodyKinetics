setwd("~/Documents/Ferret_Model/antibodyKinetics/")
#setwd("~/Drive/Influenza/Ferret/PLOS Comp Bio/")
#source("scripts/single_exposure/create_sim_data_funcs.R")
#devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics/",recompile=TRUE)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics/", recompile=FALSE)
setwd("~/Drive/Influenza/Ferret/single_exposure_analysis/")
library(ggplot2)
#paramViewer()

parTab_file <- "~/Drive/Influenza/Ferret/single_exposure_analysis/parTab_basic.csv"
exposureTab_file <- "~/Drive/Influenza/Ferret/single_exposure_analysis/exposureTab_single.csv"
nindiv <- 3
nstrain <- 1
ngroup <- 1
form <- "C"
typing <- FALSE
priming <- FALSE
cr <- FALSE
group <- strain <- NA

times <- c(0,21,37,49,70)
individuals <- rep(nindiv,ngroup)

posterior_func <- create_model_group_func_cpp(parTab,exposureTab, dat=ferret_titres1, version="posterior",
                                              form="C",typing=FALSE,cross_reactivity = FALSE,
                                              individuals=individuals)
posterior_func(exampleParTab$values)

mcmcPars <- c("adaptive_period"=20000,"iterations"=100000,"opt_freq"=2000,"thin"=10,"save_block"=10000,"popt"=0.44)

individuals <- rep(nindiv,1)

## Run the MCMC code

## Attach titre data and clean for submission to MCMC code
data(ferret_titres)
ferret_titres <- ferret_titres[ferret_titres$group == 5 & ferret_titres$strain == "A/Panama/2007/1999 (H3N2)",]
ferret_titres <- as.matrix(ferret_titres[,4:ncol(ferret_titres)])
ferret_titres <- rbind(c(0,21,37,49,70),ferret_titres)
rownames(ferret_titres) <- NULL
filename <- "test_biphasic_full"

posterior_func_mat <- create_model_group_func_cpp(parTab,exposureTab, dat=ferret_titres, version="loo",
                                                  form="C",typing=FALSE,cross_reactivity = FALSE,
                                                  individuals=individuals)
posterior_func_mat(exampleParTab$values)
library(loo)

calculate_loo <- function(tmp_run, mcmcPars, parTab_tmp, N=1000){
  chain <- read.csv(tmp_run$file)
  chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
  f <- create_model_group_func_cpp(parTab_tmp,exposureTab,version="model",
                                   form="C",typing = FALSE,cross_reactivity = FALSE)
  all_liks <- log(calculate_lik_matrix(chain, parTab_tmp, ferret_titres,f, N,nindiv=nindiv))
  rel_n_eff <- relative_eff(exp(all_liks),chain_id=rep(1,each=N))
  loo1 <- loo(all_liks,r_eff=rel_n_eff,cores=1)
  return(loo1)
}

ferret_titres_less <- ferret_titres
ferret_titres_less[2,3] <- NA
filename <- "test_biphasic_less"
run_without_point <- antibodyKinetics::run_MCMC(parTab=parTab, data=ferret_titres_less, mcmcPars=mcmcPars, 
                                  filename=filename,CREATE_POSTERIOR_FUNC=create_model_group_func_cpp,
                                  mvrPars=NULL,PRIOR_FUNC=NULL, 
                                  version="posterior",form="C",
                                  individuals=individuals,exposureTab=exposureTab,
                                  cross_reactivity=FALSE,typing=FALSE)
loo_biphasic_less <- calculate_loo(run_without_point, mcmcPars, parTab,N=1000,ferret_titres_less)

filename <- "test_biphasic_full"
run <- antibodyKinetics::run_MCMC(parTab=parTab, data=ferret_titres, mcmcPars=mcmcPars, 
                                    filename=filename,CREATE_POSTERIOR_FUNC=create_model_group_func_cpp,
                                    mvrPars=NULL,PRIOR_FUNC=NULL, 
                                    version="posterior",form="C",
                                    individuals=individuals,exposureTab=exposureTab,
                                    cross_reactivity=FALSE,typing=FALSE)
loo_biphasic_full <- calculate_loo(run, mcmcPars, parTab, N =5000)
loo_biphasic_full1 <- calculate_loo(run, mcmcPars, parTab, N =100)
loo_biphasic_full2 <- calculate_loo(run, mcmcPars, parTab, N =10000)
plot(loo_biphasic_full)
plot(loo_biphasic_full1)

filename1 <- "test_biphasic_fixed"
parTab1 <- parTab
parTab1[parTab1$names %in% c("m"),"values"] <- 0
parTab1[parTab1$names %in% c("m"),"fixed"] <- 1
run_1 <- antibodyKinetics::run_MCMC(parTab=parTab1, data=ferret_titres, mcmcPars=mcmcPars, 
                                    filename=filename1,CREATE_POSTERIOR_FUNC=create_model_group_func_cpp,
                                    mvrPars=NULL,PRIOR_FUNC=NULL, 
                                    version="posterior",form="C",
                                    individuals=individuals,exposureTab=exposureTab,
                                    cross_reactivity=FALSE,typing=FALSE)
loo_biphasic_fixed <- calculate_loo(run_1, mcmcPars, parTab1)


filename2 <- "test_monophasic"
parTab2 <- parTab
parTab2[parTab2$names %in% c("ts","dp"),"values"] <- 0
parTab2[parTab2$names %in% c("ts","dp"),"fixed"] <- 1
run_2 <- antibodyKinetics::run_MCMC(parTab=parTab2, data=ferret_titres, mcmcPars=mcmcPars, 
                                    filename=filename2,CREATE_POSTERIOR_FUNC=create_model_group_func_cpp,
                                    mvrPars=NULL,PRIOR_FUNC=NULL, 
                                    version="posterior",form="C",
                                    individuals=individuals,exposureTab=exposureTab,
                                    cross_reactivity=FALSE,typing=FALSE)
loo_monophasic <- calculate_loo(run_2, mcmcPars, parTab2)



filename3 <- "test_monophasic_fixed"
parTab3 <- parTab
parTab3[parTab3$names %in% c("ts","dp","m"),"values"] <- 0
parTab3[parTab3$names %in% c("ts","dp","m"),"fixed"] <- 1
run_3 <- antibodyKinetics::run_MCMC(parTab=parTab3, data=ferret_titres, mcmcPars=mcmcPars, 
                                    filename=filename3,CREATE_POSTERIOR_FUNC=create_model_group_func_cpp,
                                    mvrPars=NULL,PRIOR_FUNC=NULL, 
                                    version="posterior",form="C",
                                    individuals=individuals,exposureTab=exposureTab,
                                    cross_reactivity=FALSE,typing=FALSE)
loo_monophasic_fixed <- calculate_loo(run_3, mcmcPars, parTab3)


chain <- read.csv(run$file)
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
chain1 <- read.csv(run_1$file)
chain1 <- chain1[chain1$sampno > mcmcPars["adaptive_period"],]
chain2 <- read.csv(run_2$file)
chain2 <- chain2[chain2$sampno > mcmcPars["adaptive_period"],]
chain3 <- read.csv(run_3$file)
chain3 <- chain3[chain3$sampno > mcmcPars["adaptive_period"],]

WAIC_biphasic_full <- calculate_WAIC(chain, parTab, ferret_titres, f, 1000)
WAIC_biphasic_fixed <- calculate_WAIC(chain1, parTab1, ferret_titres, f, 1000)
WAIC_monophasic <- calculate_WAIC(chain2, parTab2, ferret_titres, f, 1000)
WAIC_monophasic_fixed <- calculate_WAIC(chain3, parTab3, ferret_titres, f, 1000)

#plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))
options <- convert_runName_to_options("CNAN3BN")
antibodyKinetics::plot_single_fit(parTab, exposureTab, ferret_titres, times, chain, options) + theme(legend.position="none")

p_biphasic_full <- plot_fit(parTab, chain) + theme(legend.position="none") + 
  ggtitle(paste0("LOOIC: ", signif(loo_biphasic_full$estimates[3,1],3), " (P_LOO: ", signif(loo_biphasic_full$estimates[2,1],3),")"))
p_biphasic_fixed <- plot_fit(parTab1, chain1) + theme(legend.position="none")+ 
  ggtitle(paste0("LOOIC: ", signif(loo_biphasic_fixed$estimates[3,1],3), " (P_LOO: ", signif(loo_biphasic_fixed$estimates[2,1],3),")"))
p_monophasic <- plot_fit(parTab2, chain2) + theme(legend.position="none")+ 
  ggtitle(paste0("LOOIC: ", signif(loo_monophasic$estimates[3,1],3), " (P_LOO: ", signif(loo_monophasic$estimates[2,1],3),")"))
p_monophasic_fixed <- plot_fit(parTab3, chain3) + theme(legend.position="none")+ 
  ggtitle(paste0("LOOIC: ", signif(loo_monophasic_fixed$estimates[3,1],3), " (P_LOO: ", signif(loo_monophasic_fixed$estimates[2,1],3),")"))
library(cowplot)
plot_grid(p_biphasic_full,p_biphasic_fixed,p_monophasic,p_monophasic_fixed,ncol=2)
compare(loo_biphasic_fixed, loo_biphasic_full, loo_monophasic, loo_monophasic_fixed)
