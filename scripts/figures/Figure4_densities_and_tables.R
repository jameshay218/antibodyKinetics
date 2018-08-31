####################################
## This script generates Figure 4 (posterior density plot), Figure 5 (cross reactivity profile plot) and Figure S3-5 (posterior densities for antigenic seniority and titre dependent boosting) in the manuscript. It also generates tables of posterior estimates for this run and a correlation plot from the MCMC chain.
## Note that these plots will not be saved as "Figure 4" etcm but with more descriptive names given the run ID.
## It assumes that the user has installed the antibodyKinetics package, and has generated
## the required MCMC chains as outlined in the package vignettes
## Author: James Hay
## Date: 11/06/2018
## NOTE: PLEASE check all file paths included in these scripts, as they are specific to my machine!
## NOTE 2: the "density plot" area should be modified if manually generating these plots. 
##         The "skip_pars" argument should be set to NULL to plot all parameters
####################################
library(ggplot2)
library(extrafont)
library(plyr)
library(reshape2)
library(cowplot)
library(coda)
library(mefa)
library(RColorBrewer)
library(antibodyKinetics)

source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/convergence_check_funcs.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/analyses/model_comparison_functions.R")
source("~/Documents/Ferret_Model/antibodyKinetics/scripts/figures/plotting_help.R")

##############################################################
## USER INPUT AREA
##############################################################
i <- 62

## Input area
res_wd <- "~/Documents/Ferret_Model/"

if(!dir.exists(paste0(res_wd,"plots/",i))) dir.create(paste0(res_wd,"plots/",i))
setwd(paste0(res_wd,"plots/",i))

## Where are the MCMC chains saved?
chain_wd_base <- "~/Documents/Ferret_Model/rerun_correct_times_28082018/outputs"

dat_file <- "~/Documents/Ferret_Model/antibodyKinetics/inputs/real_data_simple.csv"

## Number of iterations to disgard
adaptive <- 1000000

## Parameter and exposure table file locations
#parTab_loc <- "~/net/home/ferret/inputs/parTabs/"
#exposureTab_loc <- "~/net/home/ferret/inputs/exposureTabs/"
parTab_loc <- "~/Documents/Ferret_Model/antibodyKinetics/inputs/parTabs/"
exposureTab_loc <- "~/Documents/Ferret_Model/antibodyKinetics/inputs/exposureTabs/"

runs <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker.csv",stringsAsFactors=FALSE)
#runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)

## Times to solve model over
times <- c(0,21,37,49,70)
n <- 1000 ## Samples to take from chain

##############################################################
## Code running area - should need no modification
##############################################################
runName <- runs$runName[i]
runID <- runs$runID[i]
print(paste0(runID,"_",runName))
runName <- as.character(runName)

parTab_file <- paste0(parTab_loc,runs$parTab_file[i],".csv")
parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)

## What options were set for this model?
options <- antibodyKinetics::convert_runName_to_options(runName)
parTab <- antibodyKinetics::parTab_modification(parTab,options,FALSE)

chain_wd <- paste0(chain_wd_base,"/",runID,"_",runName)
chain <- as.data.frame(load_mcmc_chains(chain_wd, parTab, FALSE, 1, adaptive, FALSE, TRUE,TRUE)[["chain"]])

##########################
## Density plots
##########################
## USER NOTE - may need to change the "skipped" parameters. Change to NULL
## if you want to include all parameters
extraTheme <- theme(text=element_text(family="Arial"),axis.text.x=element_text(size=8,angle=45,hjust=1),
                    axis.text.y=element_text(size=8),axis.title.y=element_text(size=10),plot.margin=unit(c(0.5,0,0,0),"cm"))
densities_mu <- den_plot(chain, "mu", parTab, options, "Maximum homologous\n boost, μ", 15,add_priming_blank=FALSE,use_pointrange=TRUE)
densities_adjmu <- den_plot(chain, "adjmu", parTab, options, "Homologous boost after\ninitial waning, μ(1-d)", 15,use_pointrange=TRUE)
densities_ts <- den_plot(chain, "ts", parTab, options, "Duration of initial\n waning phase, ts", 20,use_pointrange=TRUE)
densities_dp <- den_plot(chain, "dp", parTab, options, "Initial proportion of\n boost lost, d", 1,use_pointrange=TRUE)
densities_m <- den_plot(chain, "m", parTab,options, "Long term waning\n rate, m", 0.3,FALSE,skip_pars=c("m","m.2","m.3"),yupper=12,use_pointrange=TRUE)
densities_m[[2]]
densities_sigma <- den_plot(chain, "sigma", parTab, options, "Cross reactivity\ngradient,σ", 10, skip_pars=c("sigma.1","sigma.2"),
                            yupper=100,add_priming_blank = FALSE,use_pointrange=TRUE)
densities_mod <- den_plot(chain, "mod", parTab, options, "Antigenic seniority \nmodifiers, ρ", 1,use_pointrange=TRUE)
densities_y0mod <- den_plot(chain, "y0_mod", parTab, options, "Titre dependence gradient, γ", ymax=1,ymin=-1,use_pointrange=TRUE)
densities_boostlim <- den_plot(chain, "boost_limit", parTab, options, "Maximum titre dependence,\ny_switch", ymax=12,use_pointrange=TRUE)

all_dens <- plot_grid(densities_mu[[2]] + extraTheme,
                      densities_adjmu[[2]] + extraTheme,
                      densities_ts[[2]] + extraTheme, 
                      densities_dp[[2]] + extraTheme,
                      densities_m[[2]] + extraTheme, 
                      densities_sigma[[2]] +
                        geom_hline(yintercept=1,col="gray30",linetype="dotted",size=0.5) +extraTheme,
                      ncol=2)

svg(paste0(runName,"_all_densities.svg"),width=5.2,height=6.5,family="Arial")
print(all_dens)
dev.off()

png(paste0(runName,"_all_densities.png"),width=5.2,height=6.5,family="Arial",units="in",res=300)
print(all_dens)
dev.off()


svg(paste0(runName,"_mod_density.svg"),width=4,height=3,family="Arial")
print(densities_mod[[2]])
dev.off()

#png(paste0(runName,"_all_densities.png"),width=9,height=10,units="in",res=300)
#all_dens
#dev.off()

cr_plots <- generate_cr_profiles(chain, options, parTab, 10)
cr_table <- cr_plots[[1]]
cr_plots <- cr_plots[[2]]

svg(paste0(runName,"_cr.svg"),width=5.2,height=5.2,family="Arial")
print(cr_plots)
dev.off()

png(paste0(runName,"_cr.png"),width=5.2,height=5.2,units="in",res=300)
print(cr_plots)
dev.off()


#cr_table$x <- round(cr_table$x,2)
#cr_table <- cr_table[cr_table$x %in% c(0,6.2,6.1,0.7,1.5),]
#cr_table$res <- paste0(signif(cr_table$mid,3), " (",signif(cr_table$lower,3),"-",signif(cr_table$upper,3),")")
#cr_table <- cr_table[,c("type","x","res")]

## Table of cross reactivity against antigenic distance
#wow <- dcast(cr_table,x~type)
write.table(cr_table,paste0(runName,"_antigenicdistances.csv"),sep=",",row.names=FALSE)



######
## y0 plot, if appropriate for this run ID
######
if(options$y0_mod){
  p <- plot_grid(densities_y0mod[[2]]+extraTheme + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) + 
              geom_hline(yintercept=0,linetype="dashed", alpha=0.3),
            densities_boostlim[[2]]+extraTheme + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) + 
              scale_y_continuous(limits=c(0,12),breaks=seq(0,12,by=2)))
  
  svg(paste0(runName,"_y0_mod_density.svg"),width=5.2,height=2.5,family="Arial")
  print(p)
  dev.off()
  
  y0_plot <- y0_mod_plot(chain, 15,8,12,2,1000,12)
  svg(paste0(runName,"_y0mod.svg"),width=5.2,height=4,family="Arial")
  print(y0_plot)
  dev.off()
}

##################################
## Plotting correlations
##################################
parTab[parTab$names %in% c("y0_mod","boost_limit"),"fixed"] <- 0
parTab[parTab$names == "sigma" & parTab$type == "infection1","fixed"] <- 1
parTab[parTab$names == "mod","fixed"][1] <- 1
parTab$order <- ""
parTab[parTab$names == "mod","order"] <- 1:4
#parTab$names <- par_names_2[parTab$names]
parTab[is.na(parTab$names),"names"] <- "tmp"
parTab$names2 <- NA
parTab$names2 <- paste0(parTab$names,".",parTab$type)
parTab[parTab$names == "mod","names2"] <- paste0(parTab[parTab$names == "mod","names"], ".",parTab[parTab$names == "mod","order"])

unfixed_chain <- chain[seq(1,nrow(chain),by=500),which(parTab$fixed == 0)+1]
colnames(unfixed_chain) <- parTab[parTab$fixed == 0,"names2"]
conditions <- which(abs(cor(unfixed_chain)) > 0.3,arr.ind=TRUE)
conditions <- conditions[conditions[,1] != conditions[,2],]
for(i in 1:nrow(conditions)){
  x <- conditions[i,1]
  y <- conditions[i,2]
  if(x < y){
    conditions[i,1] <- x
    conditions[i,2] <- y
  } else {
    conditions[i,1] <- y
    conditions[i,2] <- x
  }
}
conditions <- unique(conditions)

pair_plots <- NULL
for(i in 1:nrow(conditions)){
  tmpDat <- unfixed_chain[,conditions[i,]]
  colnames(tmpDat) <- c("x","y")
  pair_plots[[i]] <- ggplot(tmpDat) + 
    geom_point(aes(x=x,y=y)) + 
    xlab(colnames(unfixed_chain)[conditions[i,1]]) +
    ylab(colnames(unfixed_chain)[conditions[i,2]]) +
    theme_bw() +
      theme(text=element_text(size=8,family="Arial"),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8))
}

svg(paste0(runName,"_correlations.svg"),width=7.5,height=8,family="Arial")
print(do.call("plot_grid",c(pair_plots,ncol=5)))
dev.off()

## Save posterior densities for all inferred parameters


write.table(densities_mu[[1]][,c("type","mean","median","lower","upper")],"mu_densities.csv",sep=",",row.names=FALSE)
write.table(densities_adjmu[[1]][,c("type","mean","median","lower","upper")],"adjmu_densities.csv",sep=",",row.names=FALSE)
write.table(densities_m[[1]][,c("type","mean","median","lower","upper")],"m_densities.csv",sep=",",row.names=FALSE)
write.table(densities_sigma[[1]][,c("type","mean","median","lower","upper")],"sigma_densities.csv",sep=",",row.names=FALSE)
write.table(densities_ts[[1]][,c("type","mean","median","lower","upper")],"ts_densities.csv",sep=",",row.names=FALSE)
write.table(densities_dp[[1]][,c("type","mean","median","lower","upper")],"dp_densities.csv",sep=",",row.names=FALSE)
write.table(densities_mod[[1]][,c("type","mean","median","lower","upper")],"mod_densities.csv",sep=",",row.names=FALSE)

