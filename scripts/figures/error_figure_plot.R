####################################
## This script generates a figure of the observation error matrix used
## in the antibodyKinetics model fitting
## Author: James Hay
## Date: 11/06/2018
## NOTE: PLEASE check all file paths included in these scripts, as they are specific to my machine!
####################################

library(ggplot2)
library(reshape2)
library(cowplot)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
MAX_TITRE <- 12
norms <- matrix(nrow=2000,ncol=13)
for(x in 1:2000){
  for(i in 0:12){
    use <- (x-500)/100
    if(use >= MAX_TITRE) use <- MAX_TITRE
    if(use < 0) use <- 0
    norms[x,i+1] <- norm_error(use,i,1.1,12)
  }
}

unmelted <- norms
norms <- reshape2::melt(norms)
colnames(norms) = c("Var1","Var2","P(obs=k)")
p1 <- ggplot(norms) + 
  geom_raster(aes(x=Var1-1,y=Var2-1,fill=`P(obs=k)`)) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,12,by=1),labels=seq(0,12,by=1)) +
  scale_x_continuous(expand=c(0,1),breaks=seq(0,2001,by=100),labels=seq(-5,15,by=1)) +
  xlab("True log titre (k)") +
  ylab("Observed log titre (obs)") +
  theme_bw() +
  theme(axis.text=element_text(size=12,colour="black"),
        axis.title=element_text(size=16,colour="black"),
        legend.position="bottom",
        plot.margin=unit(units = "cm",x = c(0.5,0.5,0.5,0.5))) +
 scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.6,limits=c(0,1))
  

dat <- unmelted[c(110,560,1180),]
dat <- reshape2::melt(dat)
convert_names <- c("True log titre = 1.10","True log titre = 5.60","True log titre = 12.8")
dat$Var1 <- convert_names[dat$Var1]

p2 <- ggplot(dat) + geom_bar(aes(x=Var2-1,y=value),col="black",fill="red",alpha=0.8,stat="identity") + 
  facet_wrap(~Var1, ncol=1) +
  theme_bw() +
  theme(axis.title=element_text(size=16,colour="black"),
        axis.text=element_text(size=12,colour="black"),
        strip.text=element_text(size=16,colour="black")) +
  ylab("Probability of observation") +
  xlab("Observed log titre") +
  scale_y_continuous(limits=c(0,0.6),expand=c(0,0), breaks=seq(0,0.6,by=0.1)) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,12,by=1),labels=seq(0,12,by=1)) +
  theme(panel.grid.minor = element_blank())

cairo_ps("~/Documents/Ferret_Model/plots/error_plot.eps",width=10,height=6,family="Arial")
print(plot_grid(p2,p1,align="v",rel_widths=c(1,1.5)))
dev.off()

cairo_ps("~/Documents/Ferret_Model/plots/error_matrix.eps",width=7,height=6)
print(p1)
dev.off()

png("~/Documents/Ferret_Model/plots/error_matrix.png",width=7,height=6,units="in",res=300)
print(p1)
dev.off()
