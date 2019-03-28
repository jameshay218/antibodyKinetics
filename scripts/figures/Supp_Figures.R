####################################
## This script generates the point-range plots in the supplementary material
## of the antibodyKinetics manuscript, showing consistency/differences
## in inferred posterior distributions for different parameters across different
## model variants
## To colour by different groupings (eg. model form; waning type), change the "fillBy" arguments
## at the bottom in the "combined_density_plot" calls
## Saves both .png and .svg version of the figures
## Author: James Hay
## Date: 25/02/2019
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
res_wd <- "~/Drive/Influenza/Ferret/PLOS Comp Bio/figures/Supplementary/"

## Where are the parameter estimates saved? These come from running the
## scripts/all_analyses.R scripts in the package
all_estimates <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/combined_results/main_par_estimates.csv",stringsAsFactors=FALSE)
convergence <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/combined_results/model_comparison_table.csv",stringsAsFactors=FALSE)
#runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_all.csv",stringsAsFactors=FALSE)
runs$runName <- substring(runs$runName,2)
##############################################################
## Code running area - should need no modification
##############################################################
convergence <- convergence[order(convergence$δELPD.LOO),]
colnames(convergence)[1] <- "runID"
colnames(convergence)[2] <- "runName"
## Only look at best fitting model with delta WAIC < 75
convergence <- convergence[convergence$δELPD.LOO < 20,]
all_estimates <- all_estimates[all_estimates$runID %in% convergence$runID,]
all_estimates <- merge(all_estimates,convergence[,c("runID","δELPD.LOO")],id.vars="runID")
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
#all_estimates <- merge(all_estimates, runs[,c("runID","runName")],by="runID")
all_estimates$runID <- all_estimates$runName
all_estimates$runID <- substring(all_estimates$runID,2)
all_estimates$runID <- factor(all_estimates$runID, levels=unique(all_estimates$runID)[order(convergence$δELPD.LOO)])
all_estimates$runName <- factor(all_estimates$runName, levels=unique(all_estimates$runName)[order(convergence$δELPD.LOO)])

tmp1 <- all_estimates[all_estimates$Parameter.name %in% c("γ"),]
tmp1$Parameter.name <- as.character(tmp1$Parameter.name)
p3 <- ggplot(tmp1) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`),size=0.1)+#, col=`Cross.Reactivity`)) + 
  geom_errorbar(aes(x=runID,ymax=`X97.5..CI`,ymin=`X2.5..CI`),width=0.2)+
  geom_hline(yintercept=c(0,1),linetype="dashed")+
  #geom_hline(yintercept=0,linetype="dashed",col="gray50",size=0.1)+
  theme(axis.text.x=element_text(angle=45,hjust=1),
         legend.position="none") +
  ylab(paste0("Estimate for γ")) +
  xlab("")
  #coord_cartesian(ylim=c(0,1))
p3

tmp2 <- all_estimates[all_estimates$Parameter.name %in% c("y_limit"),]
tmp2$Parameter.name <- as.character(tmp2$Parameter.name)
p4 <- ggplot(tmp2) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`),size=0.1)+#, col=`Cross.Reactivity`)) + 
  geom_errorbar(aes(x=runID,ymax=`X97.5..CI`,ymin=`X2.5..CI`),width=0.2)+
  geom_hline(yintercept=c(0,12),linetype="dashed")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none") +
  ylab(paste0("Estimate for y_limit")) +
  xlab("")+
  scale_y_continuous(breaks=seq(0,12,by=2))+
  coord_cartesian(ylim=c(0,12))
p4

tmp3 <- all_estimates[all_estimates$Parameter.name %in% c("β"),]
tmp3$Parameter.name <- as.character(tmp3$Parameter.name)
beta_p <- ggplot(tmp3) + 
  geom_pointrange(aes(x=runID,y=Median,ymax=`X97.5..CI`,ymin=`X2.5..CI`),size=0.1) + 
  geom_errorbar(aes(x=runID,ymax=`X97.5..CI`,ymin=`X2.5..CI`),width=0.2)+
  geom_hline(yintercept=c(0,1,5,100),linetype="dashed",col=c("black","red","blue","black"))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none") +
  ylab(paste0("Estimate for β")) +
  xlab("")+
  coord_cartesian(ylim=c(0,10))
beta_p

all_estimates[all_estimates$Parameter.name == "y_switch","Parameter.name"] <- "y_limit"
mu_p <- combined_density_plot(all_estimates, "μ","μ",0,15,0,15,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Cross.Reactivity")
mu_p
#adjmu_p <- combined_density_plot(all_estimates, "μ(1-dp)","μ(1-d)",0,15,0,15,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy=NULL)
#adjmu_p
sd_p <- combined_density_plot(all_estimates, "sd","sd",0,2,0,10,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy=NULL)
sd_p
m_p <- combined_density_plot(all_estimates, "m","m",0,12,0,12,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Waning")
m_p
sigma_p <- combined_density_plot(all_estimates, "σ","σ",0,10,0,100,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, blueLine=5,redLine=1, 
                                 fillBy="Cross.Reactivity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))
sigma_p
dp_p <- combined_density_plot(all_estimates, "dp","d",0,1,0,1,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Titre.dependent.boosting")
dp_p
ts_p <- combined_density_plot(all_estimates, "ts","ts",0,30,0,30,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Titre.dependent.boosting")
ts_p
ymax_p <- combined_density_plot(all_estimates, "γ","γ",0,1,0,1,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy=NULL)
ymax_p
c_p <- combined_density_plot(all_estimates, "c","c",0,15,0,15,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy=NULL)
c_p
tau_p <- combined_density_plot(all_estimates, "τ","τ",0,1,0,1,saveDir=res_wd,savePNG=FALSE,saveEPS=FALSE, fillBy="Titre.dependent.boosting")
tau_p

svg(paste0(res_wd,"FigS5.svg"),width=7,height=5,family="Arial")
print(mu_p)
dev.off()

#svg(paste0(res_wd,"FigS6.svg"),width=7,height=5,family="Arial")
#print(adjmu_p)
#dev.off()

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

svg(paste0(res_wd,"FigS10.svg"),width=6.5,height=3,family="Arial")
print(beta_p)
dev.off()

svg(paste0(res_wd,"FigS11.svg"),width=6.5,height=3,family="Arial")
cowplot::plot_grid(p3,p4,ncol=2)
dev.off()


svg(paste0(res_wd,"FigS12.svg"),width=4,height=3,family="Arial")
print(tau_p)
dev.off()






# Saving png --------------------------------------------------------------


png(paste0(res_wd,"FigS5.png"),width=7,height=5,family="Arial",units="in",res=300)
print(mu_p)
dev.off()

#png(paste0(res_wd,"FigS6.png"),width=7,height=5,family="Arial",units="in",res=300)
#print(adjmu_p)
#dev.off()

png(paste0(res_wd,"FigS6.png"),width=7,height=5,family="Arial",units="in",res=300)
print(dp_p)
dev.off()

png(paste0(res_wd,"FigS7.png"),width=7,height=5,family="Arial",units="in",res=300)
print(ts_p)
dev.off()

png(paste0(res_wd,"FigS8.png"),width=7,height=5,family="Arial",units="in",res=300)
print(m_p)
dev.off()

png(paste0(res_wd,"FigS9.png"),width=7,height=5,family="Arial",units="in",res=300)
print(sigma_p)
dev.off()

png(paste0(res_wd,"FigS10.png"),width=6.5,height=3,family="Arial",units="in",res=300)
print(beta_p)
dev.off()

png(paste0(res_wd,"FigS11.png"),width=6.5,height=3,family="Arial",units="in",res=300)
cowplot::plot_grid(p3,p4,ncol=2)
dev.off()


png(paste0(res_wd,"FigS12.png"),width=4,height=3,family="Arial",units="in",res=300)
print(tau_p)
dev.off()


