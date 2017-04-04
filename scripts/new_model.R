library(antibodyKinetics)
library(microbenchmark)
setwd("~/Documents/Ferret Model/antibodyKinetics/scripts/")
source("funcs.R")
Rcpp::sourceCpp("../src/test.cpp")
parTab <- read.csv("parTab_new.csv",stringsAsFactors=FALSE)
#parTab <- parTab[parTab$group %in% c("all",5),]
f_old <- create_model_group_func(parTab)
f_new <- create_model_group_func_fast(parTab)
f_new_cpp <- create_model_group_func_cpp(parTab)
posterior <- create_posterior_group_func_cpp(parTab,y1)

posterior(parTab$values)
times <- seq(0,100,by=10)
dat <- y1 <- f_new_cpp(parTab$values, times)
y2 <- f_new(parTab$values,times)
y3 <- f_old(parTab$values, times)

t1 <- microbenchmark(f_old(parTab$values, seq(0,100,by=10)))
t2 <- microbenchmark(f_new(parTab$values, seq(0,100,by=10)))
t3 <- microbenchmark(f_new_cpp(parTab$values, seq(0,100,by=10)))
t4 <- microbenchmark(posterior(parTab$values))

y1 <- rbind(t(y1[1:5,]),t(y1[6:10,]),t(y1[11:15,]),t(y1[16:20,]),t(y1[21:25,]))

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

i <- 21
group <- 5

plot(y1[i,],type='l',ylim=c(0,15))
for(i in (i+1):(i+4)) lines(y1[i,])

tmpY <- y3[y3$group == group,]
for(i in 2:6) lines(tmpY[,i],col=i)

startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bounds"],startTab[i,"upper_bounds"])
}

mcmcPars1 <- c("iterations"=1000000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=400000,"save_block"=100)
run_1 <- run_MCMC(startTab,dat, mcmcPars1, "test",create_posterior_group_func_cpp,NULL,NULL)
chain <- read.csv(run_1$file)
bestPars <- get_best_pars(chain)
#plot(coda::as.mcmc(chain))
chain <- chain[chain$sampno > mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.8)

parTab1 <- parTab
parTab1$values <- bestPars

mcmcPars1 <- c("iterations"=1000000,"popt"=0.234,"opt_freq"=1000,"thin"=10,"adaptive_period"=400000,"save_block"=100)
run_2 <- run_MCMC(parTab1,dat, mcmcPars1, "test2",create_posterior_group_func_cpp,mvrPars,NULL)
#plot(coda::as.mcmc(chain1))
chain1 <- read.csv(run_2$file)
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period"],]
plot(coda::as.mcmc(chain1))








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
#lik <- create_model_post(parTab,dat,NULL)
#lik(parTab$values)

lik <- create_model_groups_post(parTab1,dat,NULL)
lik(parTab1$values)
startTab <- parTab1
for(i in which(parTab1$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bounds"],startTab[i,"upper_bounds"])
}

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=100000,"save_block"=100)


run_1 <- run_MCMC(startTab,dat, mcmcPars1, "test",create_posterior_group_func_cpp,NULL,NULL)

chain <- read.csv(run_1$file)
chain <- chain[chain$sampno > mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.8)

mcmcPars1 <- c("iterations"=250000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=200000,"save_block"=100)
run_2 <- run_MCMC(parTab1,dat, mcmcPars1, "test2",create_model_groups_post,mvrPars,NULL)
plot(coda::as.mcmc(chain1))
chain1 <- read.csv(run_2$file)
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period",]]

tmp_func <- create_model_group_func(parTab1)

mod <- generate_prediction_intervals_new(chain1, 100,seq(0,100,by=1),tmp_func,5)

#tmp <- melt(as.data.frame(dat),id.vars=c("time"))
mod$strain <- as.factor(mod$strain)
mod$group <- as.factor(mod$group)
ggplot(mod) + 
  geom_ribbon(aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.2) +
  geom_point(data = dat,aes(x=time,y=value,col=strain)) +
  facet_wrap(~group) +
  coord_cartesian(ylim=c(0,20))

