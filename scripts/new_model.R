library(antibodyKinetics)
setwd("~/Documents/Ferret Model/testing/")
devtools::load_all()
parTab <- read.csv("parTab_new.csv",stringsAsFactors=FALSE)
parTab[which(parTab$names %in% c("sigma","beta")),"values"] <- log(parTab[which(parTab$names %in% c("sigma","beta")),"values"])
parTab[which(parTab$values < parTab$lower_bounds),"values"] <- parTab[which(parTab$values < parTab$lower_bounds),"lower_bounds"]

f_new_cpp <- create_model_group_func_cpp(parTab)

times <- c(0,21,36,49,70)
dat <- y1 <- f_new_cpp(parTab$values, times)
individuals <- 3
dat <- floor(dat)
dat <- apply(dat, 2,function(x) matrix(rep(x, individuals),nrow=individuals,byrow=TRUE))

for(i in 1:nrow(dat)){
  for(j in 1:ncol(dat)){
    dat[i,j] <- add_noise(c("MAX_TITRE"=15,"S"=0.79,"EA"=0.2),dat[i,j])
  }
}
dat <- rbind(times, dat)
individuals <- rep(3,5)
posterior <- create_model_group_func_cpp(parTab,dat=dat,version="posterior",individuals=individuals)
posterior(parTab$values)

real_trajectory <- f_new_cpp(parTab$values, seq(0,100,by=1))
real_trajectory <- rbind(seq(0,100,by=1), real_trajectory)




startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bounds"],startTab[i,"upper_bounds"])
}

mcmcPars1 <- c("iterations"=10000000,"popt"=0.44,"opt_freq"=5000,"thin"=100,"adaptive_period"=1000000,"save_block"=5000)
run_1 <- run_MCMC(startTab,dat, mcmcPars1, "log_scale",create_model_group_func_cpp,NULL,NULL, version="posterior",individuals=individuals)
chain <- read.csv(run_1$file)
chain <- read.csv("log_scale_chain.csv")
bestPars <- get_best_pars(chain)

chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

parTab1 <- parTab
parTab1$values <- bestPars

mcmcPars1 <- c("iterations"=10000000,"popt"=0.234,"opt_freq"=10000,"thin"=100,"adaptive_period"=2000000,"save_block"=5000)
run_2 <- run_MCMC(parTab1,dat, mcmcPars1, "test2",create_model_group_func_cpp,mvrPars,NULL, version="posterior",individuals=individuals)
chain1 <- read.csv("test2_chain.csv")
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period"],]
plot(coda::as.mcmc(chain1))


mod <- generate_prediction_intervals(chain, 1000,seq(0,100,by=1),f_new_cpp,5,5)

meltedDat <- as.data.frame(dat[2:nrow(dat),])
dat <- floor(dat)
colnames(meltedDat) <- times
meltedDat$strain <- rep(rep(seq(1,5,by=1),5),each=3)
meltedDat$group <- rep(seq(1,5,by=1),each=5*3)
meltedDat <- reshape2::melt(meltedDat,id.vars=c("group","strain"))
meltedDat$variable <- as.numeric(as.character(meltedDat$variable))
meltedDat$group <- as.factor(meltedDat$group)
meltedDat$strain <- as.factor(meltedDat$strain)

realTraj <- as.data.frame(real_trajectory[2:nrow(real_trajectory),])
colnames(realTraj) <- seq(0,100,by=1)
realTraj$strain <- rep(seq(1,5,by=1),5)
realTraj$group <- rep(seq(1,5,by=1),each=5)
realTraj <- reshape2::melt(realTraj,id.vars=c("group","strain"))
realTraj$variable <- as.numeric(as.character(realTraj$variable))
realTraj$group <- as.factor(realTraj$group)
realTraj$strain <- as.factor(realTraj$strain)


ggplot(mod) + 
  geom_ribbon(aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4) +
  geom_line(data=realTraj,aes(x=variable,y=value,col=strain)) +
  geom_point(data = meltedDat,aes(x=variable,y=value,col=strain),position=position_jitter(w=1,h=0.5)) +
  facet_wrap(~group) +
  coord_cartesian(ylim=c(0,20)) + theme_bw()

