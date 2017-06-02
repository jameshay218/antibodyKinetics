setwd("~/Documents/Ferret_Model/antibodyKinetics")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
parTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/parTab_grp1_base.csv",stringsAsFactors = FALSE)
exposureTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/infections.csv",stringsAsFactors=FALSE)

parTab <- read.csv("~/Documents/Ferret_Model/Parameters/parTab_base.csv",stringsAsFactors = FALSE)
exposureTab <- read.csv("~/Documents/Ferret_Model/Parameters/exposureTab_base.csv",stringsAsFactors=FALSE)



f <- create_model_group_func_cpp(parTab,exposureTab,form="isolated",
                                 cross_reactivity=FALSE,typing=FALSE)





times1 <- seq(0,100,by=10)
times <- seq(0,100,by=0.1)
data <- f(parTab$values, times1)
y <- f(parTab$values, times)
data <- rbind(times1, data)
rownames(data) <- NULL
post <- create_model_group_func_cpp(parTab,exposureTab,dat=data,version="posterior",form="isolated",
                                    cross_reactivity=FALSE,typing=FALSE)

nstrain <- 5
ngroup <- 5

post(parTab$values)
post(pars)

y <- as.data.frame(y)
y$group <- rep(1:ngroup,each=nstrain)
y$strain <- rep(1:nstrain,ngroup)
colnames(y) <- c(times,"group","strain")
dat <- reshape2::melt(y,id.vars=c("group","strain"))
colnames(dat) <- c("group","strain","times","value")
convert_strain <- c("A","B","C","D","E")
dat$times <- as.numeric(as.character(dat$times)) -1
dat$strain <- convert_strain[dat$strain]
dat$strain <- as.factor(dat$strain)
dat$group <- as.factor(dat$group)

data <- as.data.frame(data[2:nrow(data),])
data$group <- rep(1:ngroup,each=nstrain)
data$strain <- rep(1:nstrain,ngroup)
colnames(data) <- c(times1,"group","strain")
dat1 <- reshape2::melt(data,id.vars=c("group","strain"))
colnames(dat1) <- c("group","strain","times","value")
convert_strain <- c("A","B","C","D","E")
dat1$times <- as.numeric(as.character(dat1$times)) -1
dat1$strain <- convert_strain[dat1$strain]
dat1$strain <- as.factor(dat1$strain)
dat1$group <- as.factor(dat1$group)



ggplot() + 
  geom_line(data=dat,aes(x=times,col=strain,y=value)) + 
  geom_point(data=dat1,aes(x=times,col=strain,y=value)) +
  facet_wrap(~group) + 
  theme_bw()





sum <- 0
for(group in unique(dat$group)){
  #print(group)
  for(strain in unique(dat$strain)){
    #print(strain)
    for(time in unique(dat$times)){ 
      #print(time)
      tmp <- (dat[dat$strain==strain & dat$times == time & dat$group == group,"value"] - 
                dat1[dat1$group == group & dat1$strain ==strain & dat1$times == time,"value"])^2
      if(tmp != 0){
        print(tmp)
        print(group)
        print(strain)
        print(time)
      }
      sum <- sum + tmp
    }
  }
}
sum

