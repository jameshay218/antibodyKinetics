setwd("~/Documents/antibodyKinetics")
library(devtools)
load_all()
parTab <- read.csv("scripts/parTab_grp1_base.csv",stringsAsFactors = FALSE)
exposures <- read.csv("scripts/infections_cr.csv",stringsAsFactors=FALSE)
f <- create_model_group_func_cpp(parTab,exposureTab,form="isolated",
                                 cross_reactivity=FALSE,typing=TRUE)
times <- seq(0,100,by=1)
y <- f(parTab$values, times)
dat <- reshape2::melt(y)
ggplot(dat) + geom_line(aes(x=Var2,col=Var1,group=Var1,y=value))


f1 <- create_model_group_func(parTab,exposureTab,form="competitive")
y1 <- f1(parTab$values,times)

#devtools::install_github("jameshay218/antibodyKinetics")
#library(antibodyKinetics)

devtools::load_all("~/Documents/Ferret Model/antibodyKinetics")
setwd("~/Documents/Ferret Model/test/")
dat <- read.csv("ferret_data.csv")
dat <- as.matrix(dat[,3:ncol(dat)])
form <- "isolated"
times <- c(0,21,36,49,70)

#strains <- c("A/Panama/2007/1999 (H3N2)" = "A", "A/Brisbane/10/2007 (H3N2)" = "B",
              #"A/Wisconsin/67/2005 (H3N2)" = "C", "A/Solomon Islands/3/2006 (H1N1)" = "D",
              #"A/Fukushima/141/2006 (H1N1)" = "E")
#dat$strain <- strains[dat$strain]


parTab <- read.csv("parTab_new3.csv",stringsAsFactors=FALSE)
parTab[which(parTab$names %in% c("sigma","beta")),"values"] <- log(parTab[which(parTab$names %in% c("sigma","beta")),"values"])
parTab[which(parTab$values < parTab$lower_bounds),"values"] <- parTab[which(parTab$values < parTab$lower_bounds),"lower_bounds"]

f_new_cpp <- create_model_group_func_cpp(parTab,form = form)
posterior <- create_model_group_func_cpp(parTab,dat=dat,form=form,version="posterior")

#dat <- y1 <- f_new_cpp(parTab$values, times)
#n_indiv <- 3
#dat <- floor(dat)
#dat <- apply(dat, 2,function(x) matrix(rep(x, n_indiv),nrow=n_indiv,byrow=TRUE))

#for(i in 1:nrow(dat)){
#  for(j in 1:ncol(dat)){
#    dat[i,j] <- add_noise(c("MAX_TITRE"=15,"S"=0.79,"EA"=0.2),dat[i,j])
#  }
#}
dat <- rbind(times, dat)
individuals <- rep(3,5)

real_trajectory <- f_new_cpp(parTab$values, seq(0,100,by=1))
real_trajectory <- rbind(seq(0,100,by=1), real_trajectory)


startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bound"],startTab[i,"upper_bound"])
}

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=50000,"thin"=1,
               "adaptive_period"=10000,"save_block"=1000)
run_1 <- run_MCMC(startTab,data, mcmcPars1, "log_scale_comp",create_model_group_func_cpp,
                  NULL,NULL, version="posterior",form=form,
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=TRUE,typing=TRUE)
chain <- read.csv(run_1$file)
bestPars <- get_best_pars(chain)

chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
#plot(coda::as.mcmc(chain[chain$sampno >= mcmcPars1["adaptive_period"],]))
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

parTab1 <- parTab
parTab1$values <- bestPars

mcmcPars1 <- c("iterations"=500000,"popt"=0.234,"opt_freq"=2000,"thin"=10,"adaptive_period"=50000,
               "save_block"=500)
run_2 <- run_MCMC(parTab1,data, mcmcPars1, "test2_comp",create_model_group_func_cpp,
                  mvrPars,NULL, version="posterior",form=form,
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=TRUE,typing=TRUE)
  
  #run_MCMC(parTab1,dat, mcmcPars1, "test2_comp",create_model_group_func_cpp,mvrPars,NULL, form = "form",version="posterior",individuals=individuals)
chain1 <- read.csv("test2_comp_chain.csv")
chain1 <- chain1[chain1$sampno > mcmcPars1["adaptive_period"],]


nindiv <- 3
nstrain <- 5
ngroup <- 5
mod <- generate_prediction_intervals(chain1, 1000,seq(0,100,by=1),f_isolated,nstrains=nstrain,ngroups=ngroup)


meltedDat <- as.data.frame(data_isolated[2:nrow(data_isolated),])
data <- floor(data)
colnames(meltedDat) <- times
meltedDat <- cbind(meltedDat,expand.grid("indiv"=1:nindiv,"strain"=1:nstrain,"group"=1:ngroup))
meltedDat <- reshape2::melt(meltedDat,id.vars=c("indiv","strain","group"))
meltedDat$variable <- as.numeric(as.character(meltedDat$variable))
meltedDat$group <- as.factor(meltedDat$group)
meltedDat$strain <- as.factor(meltedDat$strain)
meltedDat$indiv <- as.factor(meltedDat$indiv)

ggplot() + 
  geom_ribbon(data = mod, aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4)+
  geom_point(data = meltedDat,aes(x=variable,y=value,col=strain),position=position_jitter(w=0.1,h=0.)) +
  facet_wrap(~group) +
  #coord_cartesian(ylim=c(0,20)) + 
  theme_bw()



realTraj <- as.data.frame(real_trajectory[2:nrow(real_trajectory),])
colnames(realTraj) <- seq(0,100,by=1)
realTraj$strain <- rep(seq(1,5,by=1),3)
realTraj$group <- rep(seq(1,3,by=1),each=5)
realTraj <- reshape2::melt(realTraj,id.vars=c("group","strain"))
realTraj$variable <- as.numeric(as.character(realTraj$variable))
realTraj$group <- as.factor(realTraj$group)
realTraj$strain <- as.factor(realTraj$strain)


