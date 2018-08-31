####################################
## This script generates the point-range plots in the supplementary material
## of the antibodyKinetics manuscript, showing consistency/differences
## in inferred posterior distributions for different parameters across different
## model variants
## To colour by different groupings (eg. model form; waning type), change the "fillBy" arguments
## at the bottom in the "combined_density_plot" calls
## Saves both .png and .svg version of the figures
## Author: James Hay
## Date: 11/06/2018
## NOTE: PLEASE check all file paths included in these scripts, as they are specific to my machine!
####################################

library(ggplot2)
library(reshape2)
library(cowplot)
library(coda)
library(RColorBrewer)
library(antibodyKinetics)
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/figures/plotting_help.R")


##############################################################
## USER INPUT AREA
##############################################################
## Location to save outputs to
res_wd <- "~/Documents/Ferret_Model/PLOS_Pathogens/Figures_15Aug2018/"

## Where are the parameter estimates saved? These come from running the
## scripts/all_analyses.R scripts in the package
all_estimates <- read.csv("~/Documents/Ferret_Model/PLOS_Pathogens/Figures_15Aug2018/parameter_estimates.csv")
convergence <- read.csv("~/Documents/Ferret_Model/PLOS_Pathogens/Figures_15Aug2018/waic_table_complete.csv")
#runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs$runName <- substring(runs$runName,2)
##############################################################
## Code running area - should need no modification
##############################################################
convergence$deltaWAICS <- convergence$WAIC - min(convergence$WAIC)

convergence <- convergence[order(convergence$deltaWAICS),]

## Only look at best fitting model with delta WAIC < 75
convergence <- convergence[convergence$deltaWAICS < 75,]
all_estimates <- all_estimates[all_estimates$runID %in% convergence$runID,]
all_estimates <- merge(all_estimates,convergence[,c("runID")],id.vars="runID")
#all_estimates <- merge(all_estimates,runs[,c("runID","cr","wane","form")],id.vars="runID")

all_estimates <- all_estimates[order(all_estimates$runID),]
convergence <- convergence[order(convergence$runID),]
all_estimates$Exposure.Type <- as.character(all_estimates$Exposure.Type)
#all_estimates$Exposure.Type <- type_names[all_estimates$Exposure.Type]

## Code to reorder the exposure type factor levels
present_names <- ordered_names2[ordered_names2 %in% as.character(unique(all_estimates$Exposure.Type))]
all_estimates$Exposure.Type = factor(all_estimates$Exposure.Type,ordered_names2)

## The commented code below is deprecated, but it was used to colour
## the plots depending on different groupings
#cross_react_change <- c("T"="Type specific","A"="Universal")
#all_estimates$cr <- cross_react_change[as.character(all_estimates$cr)]
#colnames(all_estimates)[which(colnames(all_estimates) == "cr")] <- "Cross_reactivity"

#wane_change <- c("M"="Monophasic waning","B"="Biphasic waning")
#all_estimates$wane <- wane_change[as.character(all_estimates$wane)]
#colnames(all_estimates)[which(colnames(all_estimates) == "wane")] <- "Waning"

#form_change <- c("C"="Competitive","I"="Isolated")
#all_estimates$form <- form_change[as.character(all_estimates$form)]
#colnames(all_estimates)[which(colnames(all_estimates) == "form")] <- "Form"
all_estimates <- merge(all_estimates, runs[,c("runID","runName")],by="runID")
all_estimates$runID <- factor(all_estimates$runID, levels=unique(all_estimates$runID)[order(convergence$deltaWAICS)])
all_estimates[all_estimates$runID %in% c(65,66),"runName"] <- paste0(all_estimates[all_estimates$runID %in% c(65,66),"runName"],"_b")
all_estimates$runName <- factor(all_estimates$runName, levels=unique(all_estimates$runName)[order(convergence$deltaWAICS)])
all_estimates$runID <- all_estimates$runName
tmp <- all_estimates[all_estimates$Parameter.name %in% c("ρ1","ρ2","ρ3"),]
tmp$Parameter.name <- as.character(tmp$Parameter.name)
tmp[tmp$Parameter.name == "ρ3","Parameter.name"] <- "ρ4"
tmp[tmp$Parameter.name == "ρ2","Parameter.name"] <- "ρ3"
tmp[tmp$Parameter.name == "ρ1","Parameter.name"] <- "ρ2"
#tmp$par_name <- paste0(tmp$par_name,c("",paste0(".",1:3)))

p2 <- ggplot(tmp) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`, col=`Cross.Reactivity`)) + 
  facet_wrap(~Parameter.name,nrow=1) +
  geom_hline(yintercept=c(0,1),linetype="dashed")+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        legend.position="none") +
  ylab(paste0("Estimate for ρ")) +
  coord_cartesian(ylim=c(0,1))


tmp1 <- all_estimates[all_estimates$Parameter.name %in% c("γ"),]
tmp1$Parameter.name <- as.character(tmp1$Parameter.name)
p3 <- ggplot(tmp1) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`, col=`Cross.Reactivity`)) + 
  geom_hline(yintercept=c(-1,1),linetype="dashed")+
  geom_hline(yintercept=0,linetype="dashed",col="gray50",size=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        legend.position="none") +
  ylab(paste0("Estimate for γ")) +
  xlab("")+
  coord_cartesian(ylim=c(-1,1))

tmp2 <- all_estimates[all_estimates$Parameter.name %in% c("y_limit"),]
tmp2$Parameter.name <- as.character(tmp2$Parameter.name)
p4 <- ggplot(tmp2) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`, col=`Cross.Reactivity`)) + 
  geom_hline(yintercept=c(0,12),linetype="dashed")+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        legend.position="none") +
  ylab(paste0("Estimate for y_limit")) +
  xlab("") +
  coord_cartesian(ylim=c(0,12))

tmp3 <- all_estimates[all_estimates$Parameter.name %in% c("β"),]
tmp3$Parameter.name <- as.character(tmp3$Parameter.name)
beta_p <- ggplot(tmp3) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`)) + 
  geom_hline(yintercept=c(0,1,5,10),linetype="dashed",col=c("black","red","blue","black"))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none") +
  ylab(paste0("Estimate for β")) +
  xlab("") +
  coord_cartesian(ylim=c(0,10))


all_estimates[all_estimates$Parameter.name == "y_switch","Parameter.name"] <- "y_limit"
mu_p <- combined_density_plot(all_estimates, "μ","μ",0,15,0,15,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
adjmu_p <- combined_density_plot(all_estimates, "μ(1-dp)","μ(1-d)",0,15,0,15,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
sd_p <- combined_density_plot(all_estimates, "sd","sd",0,2,0,10,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy=NULL)
m_p <- combined_density_plot(all_estimates, "m","m",0,0.25,0,12,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
sigma_p <- combined_density_plot(all_estimates, "σ","σ",0,10,0,100,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, blueLine=5,redLine=1, fillBy="Cross.Reactivity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))
dp_p <- combined_density_plot(all_estimates, "dp","d",0,1,0,1,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
ts_p <- combined_density_plot(all_estimates, "ts","ts",0,20,0,20,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
ymax_p <- combined_density_plot(all_estimates, "γ","γ",-1,1,-1,1,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
c_p <- combined_density_plot(all_estimates, "c","c",0,15,0,15,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy=NULL)

svg(paste0(res_wd,"FigS4.svg"),width=7,height=5,family="Arial")
print(mu_p)
dev.off()

svg(paste0(res_wd,"FigS5.svg"),width=7,height=5,family="Arial")
print(adjmu_p)
dev.off()

svg(paste0(res_wd,"FigS6.svg"),width=7,height=5,family="Arial")
print(dp_p)
dev.off()

svg(paste0(res_wd,"FigS7.svg"),width=7,height=5,family="Arial")
print(ts_p)
dev.off()

svg(paste0(res_wd,"FigS8.svg"),width=7,height=5,family="Arial")
print(m_p)
dev.off()

svg(paste0(res_wd,"FigS9.svg"),width=7,height=5,family="Arial")
print(sigma_p)
dev.off()

svg(paste0(res_wd,"FigS11.svg"),width=6.5,height=3,family="Arial")
print(p2)
dev.off()

svg(paste0(res_wd,"FigS12.svg"),width=6.5,height=3,family="Arial")
cowplot::plot_grid(p3,p4,ncol=2)
dev.off()

svg(paste0(res_wd,"FigS10.svg"),width=6.5,height=3,family="Arial")
print(beta_p)
dev.off()
