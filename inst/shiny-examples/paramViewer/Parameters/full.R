setwd("~/Documents/Ferret_Model/antibodyKinetics")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")

parTab <- read.csv("~/Documents/Ferret_Model/Parameters/parTab_full.csv",stringsAsFactors = FALSE)
exposureTab <- read.csv("~/Documents/Ferret_Model/Parameters/exposureTab_full.csv",stringsAsFactors=FALSE)

times <- seq(0,100,by=10)
#times <- c(0,21,37,49,56)
times1 <- seq(0,100,by=10)
nstrain <- 5
ngroup <- 5
nindiv <- 3
individuals <- rep(nindiv,ngroup)
cr <- TRUE
typing <- TRUE

pars <- c("S"=0.79,"EA"=0.2,"MAX_TITRE"=15)


f_isolated <- create_model_group_func_cpp(parTab,exposureTab,form="isolated",
                                 cross_reactivity=cr,typing=typing)


data_isolated <- f_isolated(parTab$values, times)
data_isolated <- floor(data_isolated)
data_isolated <- apply(data_isolated,2,function(x) rep(x, each=nindiv))
for(i in 1:nrow(data_isolated)){
  for(j in 1:ncol(data_isolated)){
    data_isolated[i,j] <- add_noise(pars,data_isolated[i,j])
  }
}
data_isolated <- rbind(times, data_isolated)
rownames(data_isolated) <- NULL
y_isolated <- f_isolated(parTab$values, times1)

p_isolated <- create_model_group_func_cpp(parTab,exposureTab,dat=data_isolated,version="posterior",form="isolated",individuals = rep(nindiv,nstrain),
                                          cross_reactivity=cr,typing=typing)
print(p_isolated(parTab$values))


f_competitive <- create_model_group_func_cpp(parTab,exposureTab,form="competitive",
                                          cross_reactivity=cr,typing=typing)
data_competitive <- f_competitive(parTab$values, times)
data_competitive <- floor(data_competitive)
data_competitive <- apply(data_competitive,2,function(x) rep(x, each=nindiv))
for(i in 1:nrow(data_competitive)){
  for(j in 1:ncol(data_competitive)){
    data_competitive[i,j] <- add_noise(pars,data_competitive[i,j])
  }
}
data_competitive <- rbind(times, data_competitive)
rownames(data_competitive) <- NULL
y_competitive <- f_competitive(parTab$values, times1)

p_competitive <- create_model_group_func_cpp(parTab,exposureTab,dat=data_competitive,version="posterior",form="competitive",individuals = individuals,
                                          cross_reactivity=cr,typing=typing)
print(p_competitive(parTab$values))


parTab[parTab$names %in% c("mu","dp","m"),"fixed"] <- 0
parTab[parTab$names %in% c("mu"),"upper_bound"] <- 25

startTab <- parTab
for(i in which(parTab$fixed == 0)){
  startTab[i,"values"] <- runif(1,startTab[i,"lower_bound"],startTab[i,"upper_bound"])
}


mcmcPars1 <- c("iterations"=10000000,"popt"=0.44,"opt_freq"=2000,"thin"=1000,
               "adaptive_period"=5000000,"save_block"=5000)
run_1 <- run_MCMC(startTab,data_isolated, mcmcPars1, "testing_new",create_model_group_func_cpp,
                  NULL,NULL, version="posterior",form="isolated",
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=cr,typing=typing)

chain <- data.table::fread(run_1$file,data.table=FALSE)
bestPars <- get_best_pars(chain)

chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]

covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

parTab1 <- parTab
parTab1$values <- bestPars


mcmcPars1 <- c("iterations"=20000000,"popt"=0.234,"opt_freq"=2000,"thin"=1000,"adaptive_period"=5000000,
               "save_block"=5000)
run_2 <- run_MCMC(parTab1,data_isolated, mcmcPars1, "test2_comp",create_model_group_func_cpp,
                  mvrPars,NULL, version="posterior",form="isolated",
                  individuals=individuals,exposureTab=exposureTab,
                  cross_reactivity=cr,typing=typing)
chain1 <- read.csv(run_2$file)
plot(coda::as.mcmc(chain1[chain1$sampno > mcmcPars1["adaptive_period"],which(parTab1$fixed == 0) + 1]))


effectiveSize(coda::as.mcmc(chain1[chain1$sampno > mcmcPars1["adaptive_period"],which(parTab1$fixed == 0) + 1]))










## Plot isolated
y_isolated <- as.data.frame(y_isolated)
y_isolated$group <- rep(1:ngroup,each=nstrain)
y_isolated$strain <- rep(1:nstrain,ngroup)
colnames(y_isolated) <- c(times1,"group","strain")
dat <- reshape2::melt(y_isolated,id.vars=c("group","strain"))
colnames(dat) <- c("group","strain","times","value")
convert_strain <- c("A","B","C","D","E")
dat$times <- as.numeric(as.character(dat$times)) -1
dat$strain <- convert_strain[dat$strain]
dat$strain <- as.factor(dat$strain)
dat$group <- as.factor(dat$group)

data_isolated <- as.data.frame(data_isolated[2:nrow(data_isolated),])
data_isolated <- cbind(data_isolated,expand.grid(1:nindiv,1:nstrain,1:ngroup))
colnames(data_isolated) <- c(times,"indiv","strain","group")
dat1 <- reshape2::melt(data_isolated,id.vars=c("indiv","strain","group"))
colnames(dat1) <- c("indiv","strain","group","times","value")
#convert_strain <- c("A","B","C","D","E")
dat1$times <- as.numeric(as.character(dat1$times)) -1
dat1$strain <- convert_strain[dat1$strain]
dat1$strain <- as.factor(dat1$strain)
dat1$group <- as.factor(dat1$group)
dat1$indiv <- as.factor(dat1$indiv)



p1 <- ggplot() + 
  geom_line(data=dat,aes(x=times,col=strain,y=value)) + 
  geom_point(data=dat1,aes(x=times,col=strain,y=value),position = "jitter",size=0.3) +
  facet_wrap(~group) + 
  theme_bw()








## Plot competitive
y_competitive <- as.data.frame(y_competitive)
y_competitive$group <- rep(1:ngroup,each=nstrain)
y_competitive$strain <- rep(1:nstrain,ngroup)
colnames(y_competitive) <- c(times1,"group","strain")
dat_c <- reshape2::melt(y_competitive,id.vars=c("group","strain"))
colnames(dat_c) <- c("group","strain","times","value")
convert_strain <- c("A","B","C","D","E")
dat_c$times <- as.numeric(as.character(dat_c$times)) -1
dat_c$strain <- convert_strain[dat_c$strain]
dat_c$strain <- as.factor(dat_c$strain)
dat_c$group <- as.factor(dat_c$group)

data_competitive <- as.data.frame(data_competitive[2:nrow(data_competitive),])
data_competitive <- cbind(data_competitive,expand.grid(1:nindiv,1:nstrain,1:ngroup))
colnames(data_competitive) <- c(times,"indiv","strain","group")
dat1_c <- reshape2::melt(data_competitive,id.vars=c("indiv","strain","group"))
colnames(dat1_c) <- c("indiv","strain","group","times","value")
dat1_c$times <- as.numeric(as.character(dat1_c$times)) -1
dat1_c$strain <- convert_strain[dat1_c$strain]
dat1_c$strain <- as.factor(dat1_c$strain)
dat1_c$invid <- as.factor(dat1_c$indiv)
dat1_c$group <- as.factor(dat1_c$group)



p2 <- ggplot() + 
  geom_line(data=dat_c,aes(x=times,col=strain,y=value)) + 
  geom_point(data=dat1_c,aes(x=times,col=strain,y=value),size=0.3,position="jitter") +
  facet_wrap(~group) + 
  theme_bw()


