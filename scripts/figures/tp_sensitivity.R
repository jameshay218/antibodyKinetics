####################################
## This script generates the point-range plots in the supplementary material
## of the antibodyKinetics manuscript for the time to peak sensitivity analysis
## Author: James Hay
## Date: 11/06/2018
## NOTE: PLEASE check all file paths included in these scripts, as they are specific to my machine!
####################################
par_names <- c("mod1"="ρ1", "mod2"="ρ2","mod3"="ρ3","mod4"="ρ4",
               "mod"="ρ1", "mod.1"="ρ2","mod.2"="ρ3","mod.3"="ρ4",
               "adjmu"="μ(1-dp)","adjmu.1"="μ(1-dp)", "adjmu.2"="μ(1-dp)", 
               "adjmu.3"="μ(1-dp)", "adjmu.4"="μ(1-dp)", "adjmu.5"="μ(1-dp)",
               "mu"="μ","sigma"="σ","S"="sd","beta"="β",
               "c"="c","lnlike"="Log likelihood","dp"="dp","ts"="ts","m"="m","y0_mod"="γ", "boost_limit"="y_limit",
               "tau"="τ", "tp"="tp")

library(ggplot2)
library(extrafont)
library(plyr)
library(reshape2)
library(cowplot)
library(coda)
library(mefa)
library(RColorBrewer)
#library(antibodyKinetics)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics/")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/figures/plotting_help.R")
setwd("~/Drive/Influenza/Ferret/PLOS Comp Bio/figures/Supplementary/")

free_tp_res <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/free_tp/parameter_estimates.csv",stringsAsFactors = FALSE)
free_tp_res$model <- "Estimate tp"

lower_tp_res <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/lower_tp/parameter_estimates.csv",stringsAsFactors = FALSE)
lower_tp_res$model <- "Fixed lower tp"

all_res <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/main_results/parameter_estimates.csv",stringsAsFactors = FALSE)
all_res$model <- "Main fits"

all_res <- rbind(free_tp_res, lower_tp_res, all_res)
all_res <- all_res[all_res$runName %in% unique(lower_tp_res$runName),]
all_res <- all_res[all_res$par_name != "lnlike",]
all_res$runName <- sapply(all_res$runName, function(x) strsplit(x,split="C")[[1]][[2]])
all_res$lower_bound <- 0
all_res$upper_bound <- 1000
all_res[all_res$par_name == "mu","upper_bound"] <- 15
all_res[all_res$par_name == "m","upper_bound"] <- 12
all_res[all_res$par_name == "dp","upper_bound"] <- 1
all_res[all_res$par_name == "ts","upper_bound"] <- 30
all_res[all_res$par_name == "sigma","upper_bound"] <- 100
all_res[all_res$par_name == "tp","upper_bound"] <- 14
all_res[all_res$par_name == "c","upper_bound"] <- 20
all_res[all_res$par_name == "beta","upper_bound"] <- 100
all_res[all_res$par_name == "boost_limit","upper_bound"] <- 12
all_res[all_res$par_name == "y0_mod","upper_bound"] <- 1
all_res[all_res$par_name == "tau","upper_bound"] <- 1
all_res[all_res$par_name == "S","upper_bound"] <- 10
all_res$hline <- 0
all_res[all_res$par_name == "mu","hline"] <- 15
all_res[all_res$par_name == "m","hline"] <- 12
all_res[all_res$par_name == "dp","hline"] <- 1
all_res[all_res$par_name == "ts","hline"] <- 30
all_res[all_res$par_name == "sigma","hline"] <- 100
all_res[all_res$par_name == "tp","hline"] <- 14
all_res$hline2 <- 0
all_res[all_res$par_name == "tp","hline2"] <- 5
all_res[all_res$par_name == "tp" & all_res$type == "infection1","hline2"] <- 10

# Relabel some parameters and types for the plot
all_res[all_res$par_name %in% c("beta","c"),"type"] <- "All"
all_res$type <- type_names[all_res$type]
all_res$par_name <- par_names[all_res$par_name]

all_res_all <- all_res[all_res$type == "All",]
all_res_all <- all_res_all[complete.cases(all_res_all),]
all_res_rest <- all_res[all_res$type != "All",]
all_res_rest <- all_res_rest[complete.cases(all_res_rest),]


p1 <- ggplot(all_res_rest[all_res_rest$runName == "YTY6BN",]) + 
  geom_pointrange(aes(x=type, col=model,ymin=lower_quantile,ymax=upper_quantile,y=median), size=0.3,
                  position=position_dodge(width=0.5)) + 
  facet_wrap(~par_name,scales="free_y",ncol=2) +
  geom_blank(aes(y=upper_bound)) +
  geom_blank(aes(y=lower_bound)) +
  geom_hline(aes(yintercept=hline),col="gray40",linetype="dashed") +
  geom_hline(aes(yintercept=hline2),col="red",linetype="dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8,family="Arial",colour="black", angle=45,vjust=0.5, hjust=0.5),
        axis.text.y = element_text(size=8,family="Arial",colour="black"),
        axis.title.y=element_text(size=10,family="Arial",colour="black"),
        legend.position="bottom") +
  xlab("Exposure type") +
  ylab("Posterior estimate")
p2 <- ggplot(all_res_rest[all_res_rest$runName == "NAY6BY",]) + 
  geom_pointrange(aes(x=type, col=model,ymin=lower_quantile,ymax=upper_quantile,y=median), size=0.3,
                  position=position_dodge(width=0.5)) + 
  facet_wrap(~par_name,scales="free_y",ncol=3)+ 
  facet_wrap(~par_name,scales="free_y",ncol=2) +
  geom_blank(aes(y=upper_bound)) +
  geom_blank(aes(y=lower_bound)) +
  geom_hline(aes(yintercept=hline),col="gray40",linetype="dashed") +
  geom_hline(aes(yintercept=hline2),col="red",linetype="dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8,family="Arial",colour="black", angle=45,vjust=0.5, hjust=0.5),
        axis.text.y = element_text(size=8,family="Arial",colour="black"),
        axis.title.y=element_text(size=10,family="Arial",colour="black"),
        legend.position="bottom")  +
  xlab("Exposure type") +
  ylab("Posterior estimate")
p3 <- ggplot(all_res_rest[all_res_rest$runName == "YAY6BY",]) + 
  geom_pointrange(aes(x=type, col=model,ymin=lower_quantile,ymax=upper_quantile,y=median), size=0.3,
                  position=position_dodge(width=0.5)) + 
  facet_wrap(~par_name,scales="free_y",ncol=3)+ 
  facet_wrap(~par_name,scales="free_y",ncol=2) +
  geom_blank(aes(y=upper_bound)) +
  geom_blank(aes(y=lower_bound)) +
  geom_hline(aes(yintercept=hline),col="gray40",linetype="dashed") +
  geom_hline(aes(yintercept=hline2),col="red",linetype="dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8,family="Arial",colour="black", angle=45,vjust=0.5, hjust=0.5),
        axis.text.y = element_text(size=8,family="Arial",colour="black"),
        axis.title.y=element_text(size=10,family="Arial",colour="black"),
        legend.position="bottom")  +
  xlab("Exposure type") +
  ylab("Posterior estimate")

p4 <- ggplot(all_res_all) + 
  geom_pointrange(aes(x=runName,y=median,ymin=lower_quantile,ymax=upper_quantile,col=model), size=0.3,
                  position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept=lower_bound),col="red",linetype="dashed") + 
  geom_hline(aes(yintercept=upper_bound),col="gray40",linetype="dashed") + 
  geom_blank(aes(y=lower_bound)) +
  geom_blank(aes(y=upper_bound)) +
  facet_wrap(~par_name, scales="free_y")+
  theme_classic() +
  theme(axis.text.x = element_text(size=8,family="Arial",colour="black", angle=45,vjust=0.5, hjust=0.5),
        axis.text.y = element_text(size=8,family="Arial",colour="black"),
        axis.title.y=element_text(size=10,family="Arial",colour="black"),
        legend.position="bottom")  +
  xlab("Model variant") +
  ylab("Posterior estimate")

svg("YTY6BN_tp_analysis.svg",width=5.2,height=5.2, family="Arial")
print(p1)
dev.off()

png("YTY6BN_tp_analysis.png",width=5.2,height=5.2,units="in",res=300)
print(p1)
dev.off()

svg("NAY6BY_tp_analysis.svg",width=5.2,height=5.2, family="Arial")
print(p2)
dev.off()

png("NAY6BY_tp_analysis.png",width=5.2,height=5.2,units="in",res=300)
print(p2)
dev.off()

svg("YAY6BY_tp_analysis.svg",width=5.2,height=5.2, family="Arial")
print(p3)
dev.off()

png("YAY6BY_tp_analysis.png",width=5.2,height=5.2,units="in",res=300)
print(p3)
dev.off()

svg("rest_tp_analysis.svg",width=5.2,height=4, family="Arial")
print(p4)
dev.off()

png("rest_tp_analysis.png",width=5.2,height=4,units="in",res=300)
print(p4)
dev.off()


comparison_res_free_tp <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/free_tp/convergence_check.csv",stringsAsFactors = FALSE)
comparison_res_free_tp$model <- "Estimate tp"

comparison_res_lower_tp <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/lower_tp/convergence_check.csv",stringsAsFactors = FALSE)
comparison_res_lower_tp$model <- "Fixed lower tp"

comparison_res <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/main_results/convergence_check.csv",stringsAsFactors = FALSE)
comparison_res$model <- "Main fits"
comparison_res <- comparison_res[,colnames(comparison_res) != "delta_elpd_loo"]

comparison_res <- rbind(comparison_res_free_tp, comparison_res_lower_tp, comparison_res)
comparison_res <- comparison_res[comparison_res$runName %in% unique(comparison_res_free_tp$runName),]
