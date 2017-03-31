
source("funcs.R")
library(antibodyKinetics)

parTab <- read.csv("parTab.csv",stringsAsFactors=FALSE)
tmp_func <- create_model_func(parTab)

times <- seq(0,100,by=1)
test <- tmp_func(parTab$values, times)

tmp <- reshape2::melt(test)
colnames(tmp) <- c("time","strain","value")
tmp$time <- rep(times,3)
tmp$strain <- as.factor(tmp$strain)

times1 <- seq(0,100,by=15)
test1 <- tmp_func(parTab$values, times1)


tmp1 <- reshape2::melt(test1)
colnames(tmp1) <- c("time","strain","value")
tmp1$time <- rep(times1,3)
tmp1$strain <- as.factor(tmp1$strain)


p <- ggplot() + 
  geom_point(data=tmp1,aes(x=time,col=strain,y=value)) + 
  geom_line(data=tmp,aes(x=time,col=strain,y=value)) + 
  scale_y_continuous(limits=c(0,15))

dat <- test1
dat <- cbind(time=times1,dat)

lik <- create_model_post(parTab,dat,NULL)
lik(parTab$values)

lik <- create_model_groups_post(parTab,dat,NULL)

startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bounds"],startTab[i,"upper_bounds"])
}

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,"thin"=1,"adaptive_period"=50000,"save_block"=100)
run_1 <- run_MCMC(startTab,dat, mcmcPars1, "test",create_model_groups_post,NULL,NULL)
chain <- read.csv(run_1$file)
chain <- chain[chain$sampno > mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.8)

run_2 <- run_MCMC(parTab,dat, mcmcPars1, "test2",create_model_post,mvrPars,NULL)
#plot(coda::as.mcmc(chain))
chain1 <- read.csv(run_2$file)
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period",]]
mod <- generate_prediction_intervals_new(chain1, 100,seq(0,100,by=1),tmp_func,3)
wow <- melt(mod,id.vars=c("time","strain"))
tmp <- melt(as.data.frame(dat),id.vars=c("time"))
mod$strain <- as.factor(mod$strain)
ggplot(mod) + 
  geom_ribbon(aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.2) +
  geom_point(data = tmp,aes(x=time,y=value,col=variable)) +
  coord_cartesian(ylim=c(0,20))

