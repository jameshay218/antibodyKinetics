setwd("~/net/home/ferret")
source("scripts/create_data.R")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")

runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
runs$typing <- as.logical(runs$typing)
runs$cr <- as.logical(runs$cr)
runs$priming <- as.logical(runs$priming)
runs$parTab_file <- paste0("inputs/parTab_new/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$form <- as.character(runs$form)
runs$runName <- as.character(runs$runName)
times <- seq(0,100,by=15)
times <- c(0,21,36,49,70)
wows <- NULL
for(i in 1:nrow(runs)){
    wows[[i]] <- create_data(runs$runName[i], runs$parTab_file[i], runs$exposureTab_file[i], runs$form[i],
                runs$typing[i],runs$cr[i],runs$priming[i],runs$ngroup[i],runs$nstrain[i],
                runs$nindiv[i],times)
}



runs <- read.csv("~/net/home/ferret/inputs/base_runs.csv",stringsAsFactors=FALSE)
runs$typing <- as.logical(runs$typing)
runs$cr <- as.logical(runs$cr)
runs$priming <- as.logical(runs$priming)
runs$parTab_file <- paste0("inputs/parTabs/",runs$parTab_file,".csv")
runs$exposureTab_file <- paste0("inputs/exposureTabs/",runs$exposureTab_file,".csv")
runs$parTab_file <- as.character(runs$parTab_file)
runs$exposureTab_file <- as.character(runs$exposureTab_file)
runs$form <- as.character(runs$form)
runs$runName <- as.character(runs$runName)

strains <- NA
groups <- 1:5
tmp <- expand.grid(strain=strains,group=groups,runName=c("A","S"))
runsA <- dplyr::right_join(runs,tmp,by="runName")
runsA$nstrain <- 5
runsA$ngroup <- 1
chainNos <- rep(1:maxChain,each=nrow(runsA))
runsA$parTab_file <- as.character(runsA$parTab_file)
runsA$exposureTab_file <- as.character(runsA$exposureTab_file)
runsA$form <- as.character(runsA$form)
runsA$runName <- as.character(runsA$runName)

tmp2 <- expand.grid(strain=NA,group=groups,runName=c("D","V"),stringsAsFactors=FALSE)
runsB <- dplyr::right_join(runs,tmp2,by="runName")
runsB$nstrain <- 5
runsB$ngroup <- 1
runsB$parTab_file <- as.character(runsB$parTab_file)
runsB$exposureTab_file <- as.character(runsB$exposureTab_file)
runsB$form <- as.character(runsB$form)
runsB$runName <- as.character(runsB$runName)

runsA$strain <- as.character(runsA$strain)
runsB$strain <- as.character(runsB$strain)


for(i in 1:nrow(runsA)){
    create_data(runsA$runName[i], runsA$parTab_file[i], runsA$exposureTab_file[i], runsA$form[i],
                runsA$typing[i],runsA$cr[i],runsA$priming[i],runsA$ngroup[i],runsA$nstrain[i],
                runsA$nindiv[i],times,group=runsA$group[i],strain=runsA$strain[i])
}

for(i in 1:nrow(runsB)){
    create_data(runsB$runName[i], runsB$parTab_file[i], runsB$exposureTab_file[i], runsB$form[i],
                runsB$typing[i],runsB$cr[i],runsB$priming[i],runsB$ngroup[i],runsB$nstrain[i],
                runsB$nindiv[i],times,group=runsB$group[i],strain=runsB$strain[i])
}
