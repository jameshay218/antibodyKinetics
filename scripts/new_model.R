#model_func(parTab1,cr_table,order_tab,exposures,strains,times)

parTab <- read.csv("~/GitHub/antibodyKinetics/scripts/parTab.csv",stringsAsFactors=FALSE)
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
  scale_y_continuous(limits=c(0,10))

dat <- test1
dat <- cbind(time=times1,dat)


mcmcPars1 <- c("iterations"=5000,"popt"=0.44,"opt_freq"=100,"thin"=1,"adaptive_period"=1000,"save_block"=100)
run_1 <- run_MCMC(parTab,dat, mcmcPars1, "test",create_model_post,NULL,NULL)
chain <- read.csv(run_1$file)
plot(coda::as.mcmc(chain))
