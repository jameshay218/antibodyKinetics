library(ggplot2)
setwd("~/net/home/ferret")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")

save_wd <- "~/net/home/ferret/inputs_Jan2019/sim_all/isolated"
runs <- read.csv("~/net/home/ferret/inputs_Jan2019/run_tracker_isolated.csv",stringsAsFactors=FALSE)
runs$parTab_file <- paste0("inputs/parTabs/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")
times <- c(0,21,37,49,70)
all_data <- NULL
runs <- runs[runs$form == "I",]
for(i in 1:nrow(runs)){
  print(runs[i,])
  all_data[[i]] <- create_data(runs$runName[i], runs$runID[i],
                                 runs$parTab_file[i], runs$exposureTab_file[i], 
                                 ngroup=5, nstrain=5,nindiv=3,times=times,
                                  wd=save_wd, normal=TRUE)
}
runs$dat_file <- paste0("inputs_Jan2019/isolated/",runs$runID, "_",runs$runName,"_data.csv")
write.table(runs, "~/net/home/ferret/inputs_Jan2019/run_tracker_isolated.csv",sep=",",row.names=FALSE)
