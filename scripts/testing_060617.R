setwd("~/Documents/Ferret_Model/antibodyKinetics")
library(devtools)
load_all()
parTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inst/shiny-examples/paramViewer/Parameters/parTab_base.csv",stringsAsFactors = FALSE)
exposures <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inst/shiny-examples/paramViewer/Parameters/exposureTab_base.csv",stringsAsFactors=FALSE)

group <- 1
cr <- FALSE
typing <- FALSE
form <- "competitive"

exposureTab <- exposures[exposures$group == group,]
parTab <- parTab[parTab$id %in% c(NA,"all",unique(exposureTab$id)) | parTab$names == "x",]

parTab[parTab$names %in% c("mu","dp","ts","m"),"fixed"] <- 0
parTab[parTab$names %in% c("mu","dp","ts","m"),"upper_bound"] <- c(15,1,25,1)

## Create data
times <- seq(0,100,by=0.11)
f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=form,typing = typing,cross_reactivity = cr)
dat <- f(parTab$values,times)

plot(dat[1,],ylim=c(0,15),type='l')
lines(dat[2,])
lines(dat[3,])
lines(dat[4,])
lines(dat[5,])


startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bound"],startTab[i,"upper_bound"])
}

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=50000,"thin"=1,
               "adaptive_period"=10000,"save_block"=1000)
run_1 <- run_MCMC(startTab,data, mcmcPars1, "testing",create_model_group_func_cpp,
                  NULL,NULL, version="posterior",form=form,
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=cr,typing=typing)
chain <- read.csv(run_1$file)
bestPars <- get_best_pars(chain)