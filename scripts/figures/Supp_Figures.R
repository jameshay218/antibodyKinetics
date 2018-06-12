####################################
## This script generates the point-range plots in the supplementary material
## of the antibodyKinetics manuscript, showing consistency/differences
## in inferred posterior distributions for different parameters across different
## model variants
## To colour by different groupings (eg. model form; waning type), change the "fillBy" arguments
## at the bottom in the "combined_density_plot" calls
## Saves both .png and .eps version of the figures
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
res_wd <- "~/Documents/Ferret_Model/plots1/"

## Where are the parameter estimates saved? These come from running the
## scripts/all_analyses.R scripts in the package
all_estimates <- read.csv("~/Documents/Ferret_Model/results_112017/parameter_estimates.csv")
convergence <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/waic_table.csv")

#runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker.csv",stringsAsFactors=FALSE)


##############################################################
## Code running area - should need no modification
##############################################################
convergence$deltaWAICS <- convergence$WAIC - min(convergence$WAIC)

convergence <- convergence[order(convergence$deltaWAICS),]

## Only look at best fitting model with delta WAIC < 50
convergence <- convergence[convergence$deltaWAICS < 50,]
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
  facet_wrap(~Parameter.name,nrow=2) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ylab(paste0("Estimate for ρ")) +
  coord_cartesian(ylim=c(0,1))
all_estimates[all_estimates$Parameter.name == "y_switch","Parameter.name"] <- "y_limit"
mu_p <- combined_density_plot(all_estimates, "μ","μ",0,15,0,15,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Cross.Reactivity")
adjmu_p <- combined_density_plot(all_estimates, "μ(1-dp)","μ(1-dp)",0,15,0,15,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Cross.Reactivity")
sd_p <- combined_density_plot(all_estimates, "sd","sd",0,2,0,10,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy=NULL)
m_p <- combined_density_plot(all_estimates, "m","m",0,12,0,12,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Waning")
sigma_p <- combined_density_plot(all_estimates, "σ","σ",0,10,0,100,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, blueLine=5,redLine=1, fillBy="Cross.Reactivity")
dp_p <- combined_density_plot(all_estimates, "dp","dp",0,1,0,1,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Cross.Reactivity")
ts_p <- combined_density_plot(all_estimates, "ts","ts",0,20,0,20,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Model.Form")
ymax_p <- combined_density_plot(all_estimates, "γ","γ",-1,1,-1,1,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Model.Form")
boostlim_p <- combined_density_plot(all_estimates, "y_limit","y_limit",0,12,0,12,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy="Model.Form")
c_p <- combined_density_plot(all_estimates, "c","c",0,15,0,15,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, fillBy=NULL)
beta_p <- combined_density_plot(all_estimates, "β","β",0,10,0,100,saveDir=res_wd,savePNG=TRUE,saveEPS=TRUE, blueLine=5,redLine=1, fillBy=NULL)

cairo_ps(paste0(res_wd,"mod_densities.eps"),width=7,height=5,family="Arial")
print(p2)
dev.off()
png(paste0(res_wd,"mod_densities.png"),width=7,height=5,unit="in",res=300,family="Arial")
print(p2)
dev.off()
