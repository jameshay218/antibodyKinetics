####################################
## This script generates Figure 3 in the manuscript: the model trajectory plot.
## It assumes that the user has installed the antibodyKinetics package, and has generated
## the required MCMC chains as outlined in the package vignettes
## Author: James Hay
## Date: 11/06/2018
## NOTE: PLEASE check all file paths included in these scripts, as they are specific to my machine!
####################################

library(ggplot2)
library(reshape2)
library(cowplot)
library(coda)
library(extrafont)
library(RColorBrewer)
library(antibodyKinetics)
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/figures/plotting_help.R")

##############################################################
## USER INPUT AREA
##############################################################

## Which run ID would you like to plot?
i <- 62
use_multi <- FALSE ## use multivariate chain?

## Where to save plot to?
res_wd <- "~/Documents/Ferret_Model/"
if(!dir.exists(paste0(res_wd,"plots/",i))) dir.create(paste0(res_wd,"plots/",i))
setwd(paste0(res_wd,"plots/",i))
#chain_wd_base <- "~/Documents/Ferret_Model/raw_results_test/outputs_real/"
## Where are the MCMC chains saved?
#chain_wd_base <- "~/Documents/Ferret_Model/results_112017/outputs"
chain_wd_base <- "/media/james/Storage 2/ferrets_17Aug2018"

## Data and exposure table for plot
infection_times <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/figures/infection_times.csv",stringsAsFactors=FALSE)
#dat_file <- "~/net/home/ferret/inputs/real_data_simple.csv"
dat_file <- "~/Documents/Ferret_Model/antibodyKinetics/inputs/real_data_simple.csv"

## Number of iterations to disgard
adaptive <- 1000000

## Parameter and exposure table file locations
#parTab_loc <- "~/net/home/ferret/inputs/parTabs/"
#exposureTab_loc <- "~/net/home/ferret/inputs/exposureTabs/"
parTab_loc <- "~/Documents/Ferret_Model/antibodyKinetics/inputs/parTabs/"
exposureTab_loc <- "~/Documents/Ferret_Model/antibodyKinetics/inputs/exposureTabs/"

#runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker.csv",stringsAsFactors=FALSE)
convergence <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/waic_table.csv",stringsAsFactors=FALSE)

## Times to solve model over
times <- c(0,21,37,49,70)
n <- 1000 ## Samples to take from chain
##############################################################

##############################################################
## Code running area - should need no modification
##############################################################
runName <- runs$runName[i]
runID <- runs$runID[i]
print(paste0(runID,"_",runName))

parTab_file <- paste0(parTab_loc,runs$parTab_file[i],".csv")
parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
exposureTab <- read.csv(paste0(exposureTab_loc,runs$exposureTab_file[i],".csv"),stringsAsFactors=FALSE)

## What options were set for this model?
options <- antibodyKinetics::convert_runName_to_options(runName)
parTab <- antibodyKinetics::parTab_modification(parTab,options,FALSE)

chain_wd <- paste0(chain_wd_base,"/",runID,"_",runName)
chain <- as.data.frame(load_mcmc_chains(chain_wd, parTab, FALSE, 1, adaptive, use_multi, TRUE,TRUE)[["chain"]])

dat <- read.csv(dat_file)
dat <- as.matrix(rbind(times, dat))
rownames(dat) <- NULL


f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=options$form,
                                 typing = TRUE,cross_reactivity = TRUE)

bestPars <- get_best_pars(chain)
#f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=runs$form[i],
#                                 typing = runs$typing[i],cross_reactivity = runs$cr[i])

mod <- generate_prediction_intervals(chain, 500,seq(0,80,by=1),f,nstrains=5,ngroups=5)

nindiv <- 3
nstrain <- 5
ngroup <- 5
meltedDat <- as.data.frame(dat[2:nrow(dat),])
colnames(meltedDat) <- times
meltedDat <- cbind(meltedDat,expand.grid("indiv"=1:nindiv,"strain"=1:nstrain,"group"=1:ngroup))
meltedDat <- reshape2::melt(meltedDat,id.vars=c("indiv","strain","group"))
meltedDat$variable <- as.numeric(as.character(meltedDat$variable))
meltedDat$group <- as.factor(meltedDat$group)
meltedDat$strain <- as.factor(meltedDat$strain)
meltedDat$indiv <- as.factor(meltedDat$indiv)

bestTraj <- f(bestPars, seq(0,80,by=1))
colnames(bestTraj) <- seq(0,80,by=1)
bestTraj <- cbind(bestTraj,expand.grid("strain"=1:nstrain,"group"=1:ngroup))
bestTraj <- reshape2::melt(bestTraj,id.vars=c("strain","group"))
bestTraj$variable <- as.numeric(as.character(bestTraj$variable))
bestTraj$group <- as.factor(bestTraj$group)
bestTraj$strain <- as.factor(bestTraj$strain)

mod[mod$upper > 14,"upper"] <- 14
mod[mod$lower > 14,"lower"] <- 14
bestTraj[bestTraj$value > 14,"value"] <- 14

convert_group <- c("Group 1","Group 2", "Group 3","Group 4", "Group 5")
convert_strains <- c("A/Panama/2007/1999 (H3N2)","A/Brisbane/10/2007 (H3N2)","A/Wisconsin/67/2005 (H3N2)",
                     "A/Solomon Islands/3/2006 (H1N1)","A/Fukushima/141/2006 (H1N1)")
meltedDat$group <- convert_group[meltedDat$group]
meltedDat$strain <- convert_strains[meltedDat$strain]
mod$group <- convert_group[mod$group]
mod$strain <- convert_strains[mod$strain]
bestTraj$group <- convert_group[bestTraj$group]
bestTraj$strain <- convert_strains[bestTraj$strain]

strains <- convert_strains[c(1,2,3)]
meltedDat_A <- meltedDat[meltedDat$strain %in% strains,]
mod_A <- mod[mod$strain %in% strains,]
bestTraj_A <- bestTraj[bestTraj$strain %in% strains,]
strains <- convert_strains[c(4,5)]
meltedDat_B <- meltedDat[meltedDat$strain %in% strains,]
mod_B <- mod[mod$strain %in% strains,]
bestTraj_B <- bestTraj[bestTraj$strain %in% strains,]


xscale <- c(0,21,37,49,70, infection_times$infection)
xlabels <- c("0","21","37","49","70",paste("\n\n",infection_times$infection,sep=""))
xlabel_colours <- c(rep("gray20",5),rep("red",nrow(infection_times)))
xlabel_sizes <- c(rep(14,5),rep(10,4))

p1 <- ggplot() + 
  geom_hline(yintercept=0,linetype="dashed",col="gray") +    
  geom_hline(yintercept=12,linetype="dashed",col="gray") +
   geom_ribbon(data = mod_A, aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4)+

  geom_line(data=bestTraj_A,aes(x=variable,y=value,col=strain))+
  geom_point(data = meltedDat_A,aes(x=variable,y=value,fill=strain),size=1,shape=21,stroke=0.4,
             col="gray20",position=position_jitter(w=0.25,h=0.25)) +
    geom_vline(data=infection_times,aes(xintercept=time),col="red",linetype="dashed") +
  facet_wrap(~group,ncol=1) +
  scale_y_continuous(limits=c(-1,14),expand=c(0,0),breaks=seq(0,14,by=2)) +
  scale_x_continuous(limits=c(0,80)) +
  #geom_hline(yintercept=13,linetype="dashed",col="gray") +
  ylab("log titre") +
  xlab("Time (days)") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))+
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") + 
  scale_colour_manual(values=brewer.pal(8,"Set1")[1:3])+
  scale_fill_manual(values=brewer.pal(8,"Set1")[1:3])+
  theme(strip.background = element_blank(),
        strip.text=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.direction = "vertical",
        axis.text=element_text(family="Arial"),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line=element_line(colour="gray20"),
        axis.line.x = element_line(colour = "gray20"),
        axis.line.y=element_line(colour="gray20"),
        plot.margin = unit(c(1, 0, 0, 0), "cm"),
        panel.spacing=unit(1,"lines"),
        panel.background=element_blank())

p2 <- ggplot() + 
  geom_hline(yintercept=0,linetype="dashed",col="gray") +
  geom_hline(yintercept=12,linetype="dashed",col="gray") +
  geom_ribbon(data = mod_B, aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4)+
  
  geom_line(data=bestTraj_B,aes(x=variable,y=value,col=strain))+
  geom_point(data = meltedDat_B,aes(x=variable,y=value,fill=strain),size=1,shape=21,stroke=0.4,
             col="gray20",position=position_jitter(w=0.25,h=0.25)) +
   geom_vline(data=infection_times,aes(xintercept=time),col="red",linetype="dashed") +
  facet_wrap(~group,ncol=1) +
  scale_y_continuous(limits=c(-1,14),expand=c(0,0),breaks=seq(0,14,by=2)) +
  scale_x_continuous(limits=c(0,80)) +
  #geom_hline(yintercept=13,linetype="dashed",col="gray") +
  ylab("log titre") +
  xlab("Time (days)") +  
  guides(colour=guide_legend(nrow=1,byrow=TRUE))+
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") + 
  scale_colour_manual(values=brewer.pal(8,"Set1")[4:5])+
  scale_fill_manual(values=brewer.pal(8,"Set1")[4:5])+
  theme(strip.background = element_blank(),
        strip.text=element_blank(),
        legend.position="bottom",
        legend.direction = "vertical",
        legend.title=element_blank(),
        axis.text=element_text(family="Arial"),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line=element_line(colour="gray20"),
        axis.line.x = element_line(colour = "gray20"),
        axis.line.y=element_line(colour="gray20"),
        #axis.text.x=element_text(colour=xlabel_colours,size=xlabel_sizes),
        panel.spacing=unit(1,"lines"),
        plot.margin = unit(c(1, 0, 0, 0), "cm"),
        panel.background=element_blank())
trajP <- plot_grid(p1,p2,ncol=2,align="hv")
trajP
svg(paste0(runName,"_model_traj.svg"),width=5.2,height=6,family="Arial")
print(trajP)
dev.off()
#tiff(paste0(runName,"_model_traj.tiff"),width=7,height=7,units="in",res=300)
#print(trajP)
#dev.off()

