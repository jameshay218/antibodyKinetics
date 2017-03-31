model_trajectory <- function(pars, times){
  mu <- pars["mu"]*pars["cr"]*pars["mod"]
  dp <- pars["dp"]
  tp <- pars["tp"]
  ts <- pars["ts"]
  m <- pars["m"]
  y0 <- pars["y0"]
  t_i <- pars["t_i"]
  
  y <- numeric(length(times))
  i <- 1
  for(t in times){
    y[i] <- (
        (t <= t_i)*0
      ) +
        (
          (t > t_i)*(t <= (tp + t_i))*((mu/tp)*t-(mu/tp)*t_i)
        ) +
        (
          (t > (tp + t_i))*(t <= (ts + tp + t_i))*((-(dp*mu)/ts)*t+((mu*dp)/ts)*(t_i+tp) + mu)
        ) +
        (
          (t > (ts + tp + t_i))*(-m*t+m*(t_i+tp+ts)+(1-dp)*mu)
        ) + y0
      i <- i + 1
    }
  return(y)
}


## For one strain
model_func <- function(parTab, cr_table, order_tab, exposures, strains, times){
  trajectories <- matrix(0, ncol=length(strains),nrow=length(times))
  index <- 1
  
  new_strains <- character(length(strains))

  for(strain in strains){
    for(i in 1:nrow(exposures)){
      t_i <- exposures[i,"values"]
      type <- exposures[i,"type"]
      exposure <- exposures[i,"exposure"]
      order <- exposures[i,"order"]
     
      tmpParTab <- parTab[parTab$type==type,]
      pars <- tmpParTab$values
      names(pars) <- tmpParTab$names

      mod <- order_tab[order_tab$order == order,"values"]
      
      ind <- sort(c(exposure,strain))
      
      cr <- cr_table[cr_table$exposure==ind[1] & cr_table$strain == ind[2] & cr_table$type == type,"values"]

      pars <- c("t_i"=t_i,pars,"mod"=mod,"cr"=cr,y0=0)

      y <- model_trajectory(pars, times)

      trajectories[,index] <- trajectories[,index] + y
      
    }
    index <- index + 1
    new_strains[index] <- strain
  }
  colnames(trajectories) <- strains
  return(trajectories)
}


create_model_func <- function(parTab){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names =="t_i",]
  exposure_indices <- which(parTab$names =="t_i")

  cr_table <- parTab[parTab$names == "cr",]
  cr_indices <- which(parTab$names == "cr")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","cr","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","cr","mod")))

  test <- function(pars,times){
    parTab1$values <- pars[parTab_indices]
    cr_table$values <- pars[cr_indices]
    order_tab$values <- pars[order_indices]
    exposures$values <- pars[exposure_indices]
    y <- model_func(parTab1,cr_table,order_tab,exposures,strains,times)
    return(y)
  }
  return(test)  
}

#model_func(parTab1,cr_table,order_tab,exposures,strains,times)

parTab <- read.csv("~/testing_mcmc/parTab.csv",stringsAsFactors=FALSE)
tmp_func <- create_model_func(parTab)

times <- seq(0,100,by=1)
test <- tmp_func(parTab$values, times)

tmp <- reshape2::melt(test)
colnames(tmp) <- c("time","strain","value")
tmp$time <- rep(times,3)
tmp$strain <- as.factor(tmp$strain)

times1 <- seq(0,100,by=15)
test1 <- tmp_func(parTab$values, times1)


plot(wow[[1]][,2],type='l',col="red",ylim=c(0,10))
lines(wow[[1]][,3],col="red")
lines(wow[[1]][,4],col="red")

points(dat[,c("time","A")],col="red")
points(dat[,c("time","B")],col="blue")
points(dat[,c("time","C")],col="green")

lines(wow[[2]][,2],col="blue")
lines(wow[[2]][,3],col="blue")
lines(wow[[2]][,4],col="blue")

lines(wow[[3]][,2],col="green")
lines(wow[[3]][,3],col="green")
lines(wow[[3]][,4],col="green")





tmp1 <- reshape2::melt(test1)
colnames(tmp1) <- c("time","strain","value")
tmp1$time <- rep(times1,3)
tmp1$strain <- as.factor(tmp1$strain)


p <- ggplot() + 
  geom_point(data=tmp1,aes(x=time,col=strain,y=value)) + 
  geom_line(data=tmp,aes(x=time,col=strain,y=value)) + 
  scale_y_continuous(limits=c(0,10))

create_model_post <- function(parTab, data, PRIOR_FUNC=NULL){
  times <- data[,"time"]
  
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names=="t_i",]
  exposure_indices <- which(parTab$names=="t_i")
  
  cr_table <- parTab[parTab$names == "cr",]
  cr_indices <- which(parTab$names == "cr")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","cr","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","cr","mod")))
  
  test <- function(pars){
    parTab1$values <- pars[parTab_indices]
    cr_table$values <- pars[cr_indices]
    order_tab$values <- pars[order_indices]
    exposures$values <- pars[exposure_indices]
    y <- model_func(parTab1,cr_table,order_tab,exposures,strains,times)
    ln <- 0
    for(strain in strains){
      ln <- ln + sum(dnorm(x=data[,strain],mean=y[,strain],sd=1,log=TRUE))
    }
    
    return(ln)
  }
  return(test)  
}

dat <- test1
dat <- cbind(time=times1,dat)


library(antibodyKinetics)
mcmcPars1 <- c("iterations"=5000,"popt"=0.44,"opt_freq"=100,"thin"=1,"adaptive_period"=1000,"save_block"=100)
run_1 <- run_MCMC(parTab,dat, mcmcPars1, "test",create_model_post,NULL,NULL)

chain <- read.csv(run_1$file)
plot(coda::as.mcmc(chain))
