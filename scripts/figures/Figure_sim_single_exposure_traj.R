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

times_protocol <- c(0,21,37,49,70)
times_weekly <- seq(0,70,by=7)
times_frequent <- seq(0,72,by=3)
n <- 1000 ## Samples to take from chain

setwd("/media/james/Storage 2/ferret_results/plos_path_review/sim_single_exposure/")
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_sim.csv",stringsAsFactors=FALSE)
#runs$runName <- "CNAN3BN"
runs$run_loc <- paste0(runs$runID,"_",runs$runName)
runs$parTab_file <- paste0("~/net/home/ferret/inputs_Jan2019/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("~/net/home/ferret/inputs_Jan2019/",runs$exposureTab_file,".csv")
runs$dat_file <- paste0("~/net/home/ferret/inputs_Jan2019/",runs$dat_file)

parTab <- read.csv(runs$parTab_file[1],stringsAsFactors=FALSE)
exposureTab <- read.csv(runs$exposureTab_file[1],stringsAsFactors=FALSE)

chain_protocol_3 <- as.data.frame(load_mcmc_chains(runs$run_loc[1], parTab, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])
chain_frequent_3 <- as.data.frame(load_mcmc_chains(runs$run_loc[2], parTab, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])
chain_frequent_10 <- as.data.frame(load_mcmc_chains(runs$run_loc[3], parTab, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])
chain_weekly_3 <- as.data.frame(load_mcmc_chains(runs$run_loc[4], parTab, FALSE, 1, 100000, FALSE, TRUE,TRUE)[["chain"]])

dat_protocol_3 <- read.csv(runs$dat_file[1])
dat_protocol_3 <- rbind(times_protocol, dat_protocol_3)
dat_frequent_3 <- read.csv(runs$dat_file[2])
dat_frequent_3 <- rbind(times_frequent, dat_frequent_3)
dat_frequent_10 <- read.csv(runs$dat_file[3])
dat_frequent_10 <- rbind(times_frequent, dat_frequent_10)
dat_weekly_3 <- read.csv(runs$dat_file[4])
dat_weekly_3 <- rbind(times_weekly, dat_weekly_3)

options1 <- convert_runName_to_options("CNAN3BN")
f_3 <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                 form=options1$form,typing = TRUE,cross_reactivity = options1$cr,
                                 rep(3, 1))

f_10 <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                           form=options1$form,typing = TRUE,cross_reactivity = options1$cr,
                                           rep(10, 1))
p1 <- plot_single_fit(parTab, exposureTab, dat_protocol_3, times_protocol, chain_protocol_3, options1,TRUE,5000) + xlab("") + ylab("") +
   theme(plot.margin=margin(10,10,1,1),
        plot.title=element_text(family="Arial",size=10))
p2 <- plot_single_fit(parTab, exposureTab, dat_frequent_3, times_frequent, chain_frequent_3, options1,TRUE,5000) + xlab("") + ylab("") +
  theme(plot.margin=margin(10,10,1,1),
        plot.title=element_text(family="Arial",size=10))
p3 <- plot_single_fit(parTab, exposureTab, dat_frequent_10, times_frequent, chain_frequent_10, options1,TRUE,5000,nindiv=10) + xlab("") + ylab("") +
  theme(plot.margin=margin(10,10,1,1),
        plot.title=element_text(family="Arial",size=10))
p4 <- plot_single_fit(parTab, exposureTab, dat_weekly_3, times_weekly, chain_weekly_3, options1,TRUE,5000) + xlab("") + ylab("") +
  theme(plot.margin=margin(10,10,1,1),
        plot.title=element_text(family="Arial",size=10))
library(cowplot)
all_p <- plot_grid(p1,p4,p2,p3,ncol=2)
svg(paste0("sim_single_model_traj.svg"),width=5.2,height=4,family="Arial")
print(all_p)
dev.off()
#  ggtitle(paste0("LOOIC: ", signif(loo_biphasic_full$estimates[3,1],3), " (P_LOO: ", signif(loo_biphasic_full$estimates[2,1],3),")"))
