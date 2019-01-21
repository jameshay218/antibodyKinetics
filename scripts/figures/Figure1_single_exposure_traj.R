library(ggplot2)
library(reshape2)
library(cowplot)
library(coda)
library(extrafont)
library(RColorBrewer)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/figures/plotting_help.R")

times <- c(0,21,37,49,70)
n <- 1000 ## Samples to take from chain

setwd("/media/james/Storage 2/ferret_results/plos_path_review/single_exposure/")
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker.csv",stringsAsFactors=FALSE)
#runs$runName <- "CNAN3BN"
runs$run_loc <- paste0(runs$runID,"_",runs$runName)
runs$parTab_file <- paste0("~/net/home/ferret/inputs_Jan2019/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("~/net/home/ferret/inputs_Jan2019/",runs$exposureTab_file,".csv")

parTab_full <- read.csv(runs$parTab_file[1],stringsAsFactors=FALSE)
parTab_mono_full <- read.csv(runs$parTab_file[2],stringsAsFactors=FALSE)
parTab_biphasic_fixed <- read.csv(runs$parTab_file[3],stringsAsFactors=FALSE)
parTab_monophasic_fixed <- read.csv(runs$parTab_file[4],stringsAsFactors=FALSE)

exposureTab <- read.csv(runs$exposureTab_file[1],stringsAsFactors=FALSE)

chain_full <- as.data.frame(load_mcmc_chains(runs$run_loc[1], parTab_full, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])
chain_mono_full <- as.data.frame(load_mcmc_chains(runs$run_loc[2], parTab_mono_full, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])
chain_biphasic_fixed <- as.data.frame(load_mcmc_chains(runs$run_loc[3], parTab_biphasic_fixed, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])
chain_monophasic_fixed <- as.data.frame(load_mcmc_chains(runs$run_loc[4], parTab_monophasic_fixed, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])

data(ferret_titres)
ferret_titres <- ferret_titres[ferret_titres$group == 5 & ferret_titres$strain == "A/Panama/2007/1999 (H3N2)",]
ferret_titres <- as.matrix(ferret_titres[,4:ncol(ferret_titres)])
ferret_titres <- rbind(c(0,21,37,49,70),ferret_titres)
rownames(ferret_titres) <- NULL
ferret_titres1 <- ferret_titres
options1 <- convert_runName_to_options("CNAN3BN")
options2 <- convert_runName_to_options("CNAN3MN")


f_full <- create_model_group_func_cpp(parTab_full,exposureTab,version="model",
                                 form=options1$form,typing = TRUE,cross_reactivity = options1$cr,
                                 rep(3, 1))
loo_full <- calculate_loo(chain_full, parTab_full, dat=ferret_titres1, f=f_full, 20000, 3)[[1]]

f_mono_full <- create_model_group_func_cpp(parTab_mono_full,exposureTab,version="model",
                                           form=options1$form,typing = TRUE,cross_reactivity = options1$cr,
                                           rep(3, 1))
loo_mono_full <- calculate_loo(chain_mono_full, parTab_mono_full, dat=ferret_titres1, f=f_mono_full, 20000, 3)[[1]]

f_biphasic_fixed <- create_model_group_func_cpp(parTab_biphasic_fixed,exposureTab,version="model",
                                                form=options1$form,typing = TRUE,cross_reactivity = options1$cr,
                                                rep(3, 1))
loo_biphasic_fixed <- antibodyKinetics::calculate_loo(chain_biphasic_fixed, parTab_biphasic_fixed, ferret_titres1, f_biphasic_fixed, 20000, 3)[[1]]

f_monophasic_fixed <- create_model_group_func_cpp(parTab_monophasic_fixed,exposureTab,version="model",
                                                  form=options1$form,typing = TRUE,cross_reactivity = options1$cr,
                                                  rep(3, 1))
loo_monophasic_fixed <- antibodyKinetics::calculate_loo(chain_monophasic_fixed, parTab_monophasic_fixed, ferret_titres1, f_monophasic_fixed, 20000, 3)[[1]]

p1 <- plot_single_fit(parTab_full, exposureTab, ferret_titres1, times, chain_full, options1,TRUE,10000) + xlab("") + ylab("") +
  ggtitle(paste0("ELPD: ", signif(loo_full$estimates[1,1],3), " (SE: ", signif(loo_full$estimates[1,2],3),")")) +
  theme(plot.margin=margin(2,10,1,1),
        plot.title=element_text(family="Arial",size=10))
p2 <- plot_single_fit(parTab_mono_full, exposureTab, ferret_titres1, times, chain_mono_full, options2,TRUE,10000) + xlab("") + ylab("") +
  ggtitle(paste0("ELPD: ", signif(loo_mono_full$estimates[1,1],3), " (SE: ", signif(loo_mono_full$estimates[1,2],3),")")) +
  theme(plot.margin=margin(2,10,1,1),
        plot.title=element_text(family="Arial",size=10))
p3 <- plot_single_fit(parTab_biphasic_fixed, exposureTab, ferret_titres1, times, chain_biphasic_fixed, options1,TRUE,10000) + xlab("") + ylab("") +
  ggtitle(paste0("ELPD: ", signif(loo_biphasic_fixed$estimates[1,1],3), " (SE: ", signif(loo_biphasic_fixed$estimates[1,2],3),")")) +
  theme(plot.margin=margin(2,10,1,1),
        plot.title=element_text(family="Arial",size=10))
p4 <- plot_single_fit(parTab_monophasic_fixed, exposureTab, ferret_titres1, times, chain_monophasic_fixed, options2,TRUE,10000) + xlab("") + ylab("") +
  ggtitle(paste0("ELPD: ", signif(loo_monophasic_fixed$estimates[1,1],3), " (SE: ", signif(loo_monophasic_fixed$estimates[1,2],3),")")) +
  theme(plot.margin=margin(2,10,1,1),
        plot.title=element_text(family="Arial",size=10))
library(cowplot)
all_p <- plot_grid(p1,p3,p2,p4,ncol=2)
svg(paste0("single_model_traj.svg"),width=5.2,height=4,family="Arial")
print(all_p)
dev.off()

all_loos <- list(loo_full, loo_biphasic_fixed, loo_mono_full, loo_monophasic_fixed)
save(all_loos, file="all_loo_estimates_single.RData")
#  ggtitle(paste0("LOOIC: ", signif(loo_biphasic_full$estimates[3,1],3), " (P_LOO: ", signif(loo_biphasic_full$estimates[1,3],3),")"))
