library(antibodyKinetics)
library(coda)
source("~/Documents/Ferret_Model/analysis/model_comparison_functions.R")
n <- 1000

#############################
## Input area
VER <- "outputs"
chain_wd <- "~/Documents/Ferret_Model/results_112017/outputs/"
burnin <- 1000000
monophasic_waning <- FALSE
y0_mod <- FALSE
antigenic_seniority <- FALSE

output_file <- "~/Documents/Ferret_Model/results_112017/model_comparison_tmp.csv"
output_file_residuals <- "~/Documents/Ferret_Model/results_112017/model_residuals.csv"

times <- c(0,21,36,49,70)
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

use_runs <- c("G","I","J","L","M","N","O","W","X","Z","AC","AF","AJ")
runs <- runs[runs$runName %in% use_runs,]

## For y0mod  run
univariate_files <- c("A","D","S","V","B","H","I","K","T","U","Y","AA","AB","AD")

bics <- NULL
names <- NULL
waics <- NULL
pwaics <- NULL
residuals <- NULL
index <- 1

for(i in 1:nrow(runs)){
    tmp <- analyses(runs$runName[i],runs$parTab_file[i],runs$exposureTab_file[i],
                    runs$form[i],runs$typing[i],runs$cr[i], y0_mod,
                    antigenic_seniority, runs$priming[i],
                    runs$dat_file[i],n,burnin,
                    univariate_files,
                    monophasic_waning, y0_mod,
                    wd="~/Documents/Ferret_Model/results/outputs_y0mod/")
    bics[index] <- tmp[["BIC"]]
    waics[index] <- tmp[["WAIC"]]
    pwaic[index] <- tmp[["pwaic"]]
    names[index] <- runs$runName[i]
    tmp_residuals <- as.data.frame(tmp[["mle_res"]])
    tmp_residuals$runName <- runs$runName[i]
    residuals <- rbind(residuals,tmp_residuals)
    index <- index + 1
}

res <- data.frame(names,bics, waics, pwaics)
write.table(res,output_file,row.names=FALSE,sep=",")
#write.table(residuals,output_file_residuals,row.names=FALSE,sep=",")
