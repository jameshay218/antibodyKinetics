source("funcs.R")
library(antibodyKinetics)

parTab1 <- read.csv("parTab_new.csv",stringsAsFactors=FALSE)
allA <- NULL
allB <- NULL
for(group in 1:5){
  parTab <- parTab1[parTab1$group %in% c("all",group),]
  tmp_func <- create_model_func(parTab)
  
  times <- seq(0,100,by=1)
  test <- tmp_func(parTab$values, times)
  
  tmp <- reshape2::melt(test)
  colnames(tmp) <- c("time","strain","value")
  tmp$time <- rep(times,5)
  tmp$strain <- as.factor(tmp$strain)
  
  times1 <- c(0,21,37,48,70)
  test1 <- tmp_func(parTab$values, times1)
  
  
  tmp1 <- reshape2::melt(test1)
  colnames(tmp1) <- c("time","strain","value")
  tmp1$time <- rep(times1,5)
  tmp1$strain <- as.factor(tmp1$strain)
  
  tmp$group <- group
  tmp1$group <- group
  
  allA <- rbind(allA,tmp)
  allB <- rbind(allB,tmp1)
}


crs <- NULL
for(type in c("infection","vacc","adj")){
  sigma <- parTab1[parTab1$names =="sigma" & parTab1$type == type,"values"]
  mu <- parTab1[parTab1$names == "mu" & parTab1$type == type,"values"]
  tmp <- exp(-sigma*0:500)*mu
  tmp <- data.frame(x=0:500,value=tmp,type=type)
  crs <- rbind(crs, tmp)
  }
tmp <- exp(-parTab1[parTab1$names == "beta","values"]*0:500)*parTab1[parTab1$names == "c","values"]
tmp <- data.frame(x=0:500,value=tmp,type="priming")
crs <- rbind(crs,tmp)
colnames(crs) <- c("time","value","strain")
crs$group <- "cross reactive boosting"
crs <- crs[,c("time","strain","value","group")]

allA <- rbind(allA,crs)
allB$value <- floor(allB$value)
p <- ggplot() + 
  geom_point(data=allB,aes(x=time,col=strain,y=value)) + 
  geom_line(data=allA,aes(x=time,col=strain,y=value)) + 
  facet_wrap(~group,scales="free_x")+
  coord_cartesian(ylim=c(-1,15)) + 
  ylab("Antibody titre/mu") +
  xlab("Time (days)/Antigenic distance") +
  theme_bw()


#dat <- test1
#dat <- cbind(time=times1,dat)
dat <- allB
lik <- create_model_post(parTab,dat,NULL)
lik(parTab$values)

lik <- create_model_groups_post(parTab1,dat,NULL)
lik(parTab$values)
startTab <- parTab1
for(i in which(parTab1$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bounds"],startTab[i,"upper_bounds"])
}

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=100000,"save_block"=1000)
run_1 <- run_MCMC(startTab,dat, mcmcPars1, "test",create_model_groups_post,NULL,NULL)
chain <- read.csv(run_1$file)
chain <- chain[chain$sampno > mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.8)

mcmcPars1 <- c("iterations"=250000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=200000,"save_block"=1000)
run_2 <- run_MCMC(parTab,dat, mcmcPars1, "test2",create_model_post,mvrPars,NULL)
#plot(coda::as.mcmc(chain))
chain1 <- read.csv(run_2$file)
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period",]]

tmp_func <- create_model_group_func(parTab1)

mod <- generate_prediction_intervals_new(chain, 100,seq(0,100,by=1),tmp_func,5)

#tmp <- melt(as.data.frame(dat),id.vars=c("time"))
mod$strain <- as.factor(mod$strain)
mod$group <- as.factor(mod$group)
ggplot(mod) + 
  geom_ribbon(aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.2) +
  geom_point(data = dat,aes(x=time,y=value,col=strain)) +
  facet_wrap(~group) +
  coord_cartesian(ylim=c(0,20))

