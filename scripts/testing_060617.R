setwd("~/Documents/Ferret_Model/antibodyKinetics")
library(devtools)
load_all()
<<<<<<< HEAD
parTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inst/shiny-examples/paramViewer/Parameters/parTab_base.csv",stringsAsFactors = FALSE)
exposures <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inst/shiny-examples/paramViewer/Parameters/exposureTab_base.csv",stringsAsFactors=FALSE)
=======
parTab <- read.csv("~/Documents/antibodyKinetics/inst/shiny-examples/paramViewer/Parameters/parTab_base.csv",stringsAsFactors = FALSE)
exposures <- exposureTab <- read.csv("~/Documents/antibodyKinetics/inst/shiny-examples/paramViewer/Parameters/exposureTab_base.csv",stringsAsFactors=FALSE)
>>>>>>> a55e9fb2ab962d974093274362fdc0de00497f12

group <- 1
cr <- FALSE
typing <- FALSE
form <- "competitive"

exposureTab <- exposures[exposures$group == group,]
parTab <- parTab[parTab$id %in% c(NA,"all",unique(exposureTab$id)) | parTab$names == "x",]

parTab[parTab$names %in% c("mu","dp","ts","m"),"fixed"] <- 0
parTab[parTab$names %in% c("mu","dp","ts","m"),"upper_bound"] <- c(15,1,25,1)

## Create data
<<<<<<< HEAD
times <- seq(0,100,by=0.11)
f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=form,typing = typing,cross_reactivity = cr)
dat <- f(parTab$values,times)

plot(dat[1,],ylim=c(0,15),type='l')
lines(dat[2,])
lines(dat[3,])
lines(dat[4,])
lines(dat[5,])
=======
times <- seq(0,100,by=10)
f <- create_model_group_func_cpp(parTab,exposureTab,form=form,version="model",typing = typing,cross_reactivity = cr)
data <- f(parTab$values,times)
plot(data[1,])
data <- rbind(times,data)
data <- floor(data)

data <- matrix(c(0,8,8,8,6),nrow=1)
times <- c(0,21,36,49,70)
data <- rbind(times,data)
>>>>>>> a55e9fb2ab962d974093274362fdc0de00497f12

post <- create_model_group_func_cpp(parTab,exposureTab,dat=data,form=form,version="posterior",typing = typing,cross_reactivity = cr,individuals =1)
post(parTab$values)

startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bound"],startTab[i,"upper_bound"])
}

mcmcPars1 <- c("iterations"=500000,"popt"=0.44,"opt_freq"=500,"thin"=10,
               "adaptive_period"=100000,"save_block"=1000)
run_1 <- run_MCMC(startTab,data, mcmcPars1, "testing",create_model_group_func_cpp,
                  NULL,NULL, version="posterior",form=form,
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=cr,typing=typing)
chain <- read.csv(run_1$file)
bestPars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

parTab1 <- parTab
parTab1$values <- bestPars

mcmcPars1 <- c("iterations"=500000,"popt"=0.234,"opt_freq"=2000,"thin"=10,"adaptive_period"=50000,
               "save_block"=500)
run_2 <- run_MCMC(parTab1,data, mcmcPars1, "test2_comp",create_model_group_func_cpp,
                  mvrPars,NULL, version="posterior",form=form,
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=cr,typing=typing)

chain1 <- read.csv(run_2$file)
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period"],]
plot(coda::as.mcmc(chain1[,which(parTab$fixed == 0)+1]))
