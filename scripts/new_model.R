library(antibodyKinetics)
library(microbenchmark)
setwd("~/Documents/Ferret Model/antibodyKinetics/scripts/")

parTab <- read.csv("parTab_new.csv",stringsAsFactors=FALSE)
parTab[which(parTab$names %in% c("sigma","beta")),"values"] <- log(parTab[which(parTab$names %in% c("sigma","beta")),"values"])
parTab[which(parTab$values < parTab$lower_bounds),"values"] <- parTab[which(parTab$values < parTab$lower_bounds),"lower_bounds"]

f_new_cpp <- create_model_group_func_cpp(parTab)

times <- seq(0,100,by=10)
dat <- y1 <- f_new_cpp(parTab$values, times)

startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bounds"],startTab[i,"upper_bounds"])
}

mcmcPars1 <- c("iterations"=10000000,"popt"=0.44,"opt_freq"=5000,"thin"=100,"adaptive_period"=1000000,"save_block"=5000)
run_1 <- run_MCMC(startTab,dat, mcmcPars1, "log_scale",create_posterior_group_func_cpp,NULL,NULL)
chain <- read.csv(run_1$file)
chain <- read.csv("log_scale_chain.csv")
bestPars <- get_best_pars(chain)

chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

parTab1 <- parTab
parTab1$values <- bestPars

mcmcPars1 <- c("iterations"=10000000,"popt"=0.1,"opt_freq"=5000,"thin"=100,"adaptive_period"=2000000,"save_block"=5000)
devtools::load_all()
run_2 <- run_MCMC(parTab1,dat, mcmcPars1, "test2",create_posterior_group_func_cpp,mvrPars,NULL)
chain1 <- read.csv("test2_chain.csv")
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period"],]
plot(coda::as.mcmc(chain1))


mod <- generate_prediction_intervals(chain1, 1000,seq(0,100,by=1),f_new_cpp,5,5)]
meltedDat <- as.data.frame(dat)
dat <- floor(dat)
colnames(meltedDat) <- times
meltedDat$strain <- rep(seq(1,5,by=1),5)
meltedDat$group <- rep(seq(1,5,by=1),each=5)
meltedDat <- reshape2::melt(meltedDat,id.vars=c("group","strain"))
meltedDat$variable <- as.numeric(as.character(meltedDat$variable))
meltedDat$group <- as.factor(meltedDat$group)
meltedDat$strain <- as.factor(meltedDat$strain)
ggplot(mod) + 
  geom_ribbon(aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4) +
  geom_point(data = meltedDat,aes(x=variable,y=value,col=strain)) +
  facet_wrap(~group) +
  coord_cartesian(ylim=c(0,20)) + theme_bw()

