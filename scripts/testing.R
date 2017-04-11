## Normal trajectory
test1_normaltrajectory<- function(){
    pars <- c("lower_bound"=0,
              "S"=0.79,
              "EA"=0.2,
              "MAX_TITRE"=13,    
              "mu"=8,
              "tp"=12,
              "dp"=0.5,
              "ts"=10,
              "m"=0.003,
              "sigma"=-6,
              "beta"=log(0.005),
              "c"=4,
              "y0_mod"=-20,
              "primed"=0,
              "mod"=1,
              "x"=0,
              "t_i"=0,
              "y0"=0,
              "eff_y0"=0)

    times <- seq(0,100,by=1)

    y <- model_trajectory_cpp(pars, times)
    y1 <- model_trajectory(pars, times, TRUE)
    return(y - y1)
}

## Primed
test2_primed <- function(){
    pars <- c("lower_bound"=0,
              "S"=0.79,
              "EA"=0.2,
              "MAX_TITRE"=13,    
              "mu"=8,
              "tp"=12,
              "dp"=0.5,
              "ts"=10,
              "m"=0.003,
              "sigma"=-6,
              "beta"=log(0.005),
              "c"=4,
              "y0_mod"=-20,
              "primed"=1,
              "mod"=1,
              "x"=0,
              "t_i"=0,
              "y0"=0,
              "eff_y0"=0)

    times <- seq(0,100,by=1)

    y <- model_trajectory_cpp(pars, times)
    y1 <- model_trajectory(pars, times, TRUE)
    return(y - y1)
}


## Cross reactivity
test3_cr<- function(){
    pars <- c("lower_bound"=0,
              "S"=0.79,
              "EA"=0.2,
              "MAX_TITRE"=13,    
              "mu"=8,
              "tp"=12,
              "dp"=0.5,
              "ts"=10,
              "m"=0.003,
              "sigma"=-6,
              "beta"=log(0.005),
              "c"=4,
              "y0_mod"=-20,
              "primed"=0,
              "mod"=1,
              "x"=100,
              "t_i"=10,
              "y0"=0,
              "eff_y0"=5)

    times <- seq(0,100,by=1)

    y <- model_trajectory_cpp(pars, times)
    y1 <- model_trajectory(pars, times, TRUE)
    return(y - y1)
}

devtools::load_all()

parTab <- read.csv("~/Documents/Ferret Model/test/parTab_new3.csv",stringsAsFactors=FALSE)
parTab <- parTab[parTab$group %in% c("all","1"),]
parTab1 <- parTab[!(parTab$names %in% c("t_i","x")),]
exposures <- parTab[parTab$names == "t_i" & parTab$group == 1,]
cr_table <- parTab[parTab$names == "x",]

order_tab <- parTab[parTab$names == "mod",]
y <- model_func(parTab1,cr_table,order_tab,exposures,strains,times,1)


colnames(y) <- strains
dat <- reshape2::melt(y)
library(ggplot2)
p1 <- ggplot() + geom_line(data=dat,aes(y=value,col=Var2,x=Var1))






devtools::load_all()
times <- seq(0,100,by=1)
strains <- unique(parTab$strain)
strains <- strains[!is.na(strains)]
parTab <- read.csv("~/Documents/Ferret Model/test/parTab_new2.csv",stringsAsFactors=FALSE)
parTab <- parTab[parTab$group %in% c("all","1"),]
f <- create_model_func(parTab, "competitive")
y <- f(parTab$values, times)
colnames(y) <- strains
dat <- reshape2::melt(y)
p2 <- ggplot() + geom_line(data=dat,aes(y=value,col=Var2,x=Var1))


f2 <- create_model_group_func_cpp(parTab,version="model",form="competitive")
y1 <- f2(parTab$values, times)
y1 <- t(y1)
colnames(y1) <- strains
dat1 <- reshape2::melt(y1)
library(ggplot2)
p1 <- ggplot() + geom_line(data=dat1,aes(y=value,col=Var2,x=Var1))

