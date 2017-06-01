setwd("~/Documents/Ferret_Model/antibodyKinetics")
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
parTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/parTab_grp1_base.csv",stringsAsFactors = FALSE)
exposureTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/infections.csv",stringsAsFactors=FALSE)

parTab1 <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/downloaded_partab.csv",stringsAsFactors = FALSE)
exposureTab1 <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/scripts/downloaded_exposures.csv",stringsAsFactors=FALSE)

parTab1 <- read.csv("~/Downloads/parTab(1).csv",stringsAsFactors = FALSE)
exposureTab1 <- read.csv("~/Downloads/exposureTab(3).csv",stringsAsFactors=FALSE)

parTab1 <- read.csv("~/Downloads/pTab.csv",stringsAsFactors = FALSE)
exposureTab1 <- read.csv("~/Downloads/expTab.csv",stringsAsFactors=FALSE)



f <- create_model_group_func_cpp(parTab1,exposureTab1,form="isolated",
                                 cross_reactivity=FALSE,typing=FALSE)





times <- seq(0,100,by=1)
y <- f(parTab1$values, times)
data <- y
post <- create_model_group_func_cpp(parTab,exposureTab,dat=data,version="posterior",form="isolated",
                                    cross_reactivity=TRUE,typing=TRUE)

post(parTab$values)
y <- as.data.frame(y)
y$group <- rep(c(1,2,3),each=5)
y$strain <- rep(c(1,2,3,4,5),3)
colnames(y) <- c(times,"group","strain")
dat <- reshape2::melt(y,id.vars=c("group","strain"))
colnames(dat) <- c("group","strain","times","value")
convert_strain <- c("A","B","C","D","E")
dat$times <- as.numeric(dat$times) -1
dat$strain <- convert_strain[dat$strain]
dat$strain <- as.factor(dat$strain)
dat$group <- as.factor(dat$group)
#dat$group <- 1
ggplot(dat) + geom_line(aes(x=times,col=strain,y=value)) + facet_wrap(~group)




f1 <- create_model_group_func(parTab,exposureTab,form="isolated",cross_reactivity=TRUE,typing=TRUE)
y1 <- f1(parTab$values,times)
dat1 <- reshape2::melt(y1,id.vars=c("times","group"))
colnames(dat1) <- c("times","group","strain","value")
ggplot(dat1) + geom_line(aes(x=times,col=strain,y=value)) + facet_wrap(~group)


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

