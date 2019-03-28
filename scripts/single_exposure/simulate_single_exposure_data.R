######################
## JAMES HAY 11.02.2019 - jameshay218@gmail.com
## This script simulates single-exposure antibody trajectories
## with different sampling protocols and model assumptions.
## Simulate data are saved in the set working directory, save_wd.
## Note that file paths are set for my local PC, and will need to be
## changed to run on a new PC. The parameter tables and exposure table
## are included in the R package.

library(ggplot2)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics/", recompile=FALSE)
parTab_file <- "~/net/home/ferret/inputs_Jan2019/single_exposure/parTab_basic.csv"
exposureTab_file <- "~/net/home/ferret/inputs_Jan2019/single_exposure/exposureTab_single.csv"
nindiv <- 3
runName <- "CNAN3BN"
save_wd <- "~/net/home/ferret/inputs_Jan2019/single_exposure/sim/"

times <- c(0,21,37,49,70)
data <- create_data(runName,"sim_protocol_3",parTab_file, exposureTab_file,
                    ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)
times <- seq(0,72,by=3)
data <- create_data(runName,"sim_frequent_3",parTab_file, exposureTab_file,
                    ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)
data <- create_data(runName,"sim_frequent_10",parTab_file, exposureTab_file,
                    ngroup=1,nstrain=1,nindiv=10,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)
times <- seq(0,70,by=7)
data <- create_data(runName,"sim_weekly_3",parTab_file, exposureTab_file,
                    ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)


parTab_biphasic_full <- parTab_monophasic_full <- "~/net/home/ferret/inputs_Jan2019/single_exposure/parTab_basic.csv"
parTab_biphasic_fixed <- "~/net/home/ferret/inputs_Jan2019/single_exposure/parTab_basic_fixed_biphasic.csv"
parTab_monophasic_fixed <- "~/net/home/ferret/inputs_Jan2019/single_exposure/parTab_basic_fixed_monophasic.csv"
exposureTab_file <- "~/net/home/ferret/inputs_Jan2019/single_exposure/exposureTab_single.csv"

times <- c(0,21,37,49,70)
data <- create_data("CNAN3BN","sim_biphasic_full",parTab_biphasic_full, exposureTab_file,
                      ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)
data <- create_data("CNAN3MN","sim_monophasic_fixed",parTab_monophasic_fixed, exposureTab_file,
                     ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)
data <- create_data("CNAN3MN","sim_monophasic_full",parTab_monophasic_full, exposureTab_file,
                      ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)
data <- create_data("CNAN3BN","sim_biphasic_fixed",parTab_biphasic_fixed, exposureTab_file,
                     ngroup=1,nstrain=1,nindiv=nindiv,
                    times=times, wd=save_wd,
                    group=NA,strain=NA,normal=TRUE)

