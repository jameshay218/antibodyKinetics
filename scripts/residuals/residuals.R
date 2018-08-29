source("~/Documents/Ferret_Model/analysis/get_residuals.R")
n <- 1000
library(coda)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
source("~/Documents/Ferret_Model/analysis/model_comparison_functions.R")

models <- read.csv("~/Dropbox/Wellcome Trust/ferret_model/results/model_comparisons_predictive_models.csv",stringsAsFactors=FALSE)

times <- c(0,21,37,49,70)
runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs$typing <- as.logical(runs$typing)
runs$cr <- as.logical(runs$cr)
runs$priming <- as.logical(runs$priming)
runs$parTab_file <- paste0("inputs/parTab_new/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")
runs$dat_file <- "inputs/real_data_simple.csv"
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$form <- as.character(runs$form)
runs$runName <- as.character(runs$runName)

all_res <- NULL
names <- models[models$bi.phasic.waning=="Yes" & models$titre_dependent_boosting=="No","names"]
tmpRuns <- runs[runs$runName %in% names,]

## Biphasic waning
for(i in 1:nrow(tmpRuns)){
  print(i)
  tmp <- get_residuals(tmpRuns$runName[i],tmpRuns$parTab_file[i],tmpRuns$exposureTab_file[i],
                                                                     tmpRuns$form[i],tmpRuns$typing[i],tmpRuns$cr[i],
                                                                     tmpRuns$dat_file[i],100,c("adaptive_period"=1000000),
                                                                     univariate_files, biphasic_waning=TRUE,y0_mod=FALSE,
                  wd="~/Documents/Ferret_Model/results/outputs_real_2/")
  all_res <- rbind(all_res,tmp[[1]])
}

## Monophasic waning
names <- models[models$bi.phasic.waning=="No" & models$titre_dependent_boosting=="No","names"]
tmpRuns <- runs[runs$runName %in% names,]

for(i in 1:nrow(tmpRuns)){
  print(i)
  tmp <- get_residuals(tmpRuns$runName[i],tmpRuns$parTab_file[i],tmpRuns$exposureTab_file[i],
                  tmpRuns$form[i],tmpRuns$typing[i],tmpRuns$cr[i],
                  tmpRuns$dat_file[i],100,c("adaptive_period"=1000000),
                  univariate_files, biphasic_waning=FALSE,y0_mod=FALSE,
                  wd="~/Documents/Ferret_Model/results/outputs_simple/")
  all_res <- rbind(all_res,tmp[[1]])
}
## y0 mod
names <- models[models$bi.phasic.waning=="Yes" & models$titre_dependent_boosting=="Yes","names"]
tmpRuns <- runs[runs$runName %in% names,]

for(i in 1:nrow(tmpRuns)){
  print(i)
  tmp <- get_residuals(tmpRuns$runName[i],tmpRuns$parTab_file[i],tmpRuns$exposureTab_file[i],
                  tmpRuns$form[i],tmpRuns$typing[i],tmpRuns$cr[i],
                  tmpRuns$dat_file[i],100,c("adaptive_period"=1000000),
                  univariate_files, biphasic_waning=TRUE,y0_mod=TRUE,
                  wd="~/Documents/Ferret_Model/results/outputs_y0mod/")
  all_res <- rbind(all_res,tmp[[1]])
}


dat <- tmp[[2]]

write.csv(all_res,"~/Dropbox/Wellcome Trust/ferret_model/results/all_residuals.csv")
write.csv(dat,"~/Dropbox/Wellcome Trust/ferret_model/results/data_residual_format.csv")
