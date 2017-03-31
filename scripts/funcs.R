library(ggplot2)
library(plyr)
#'library(antibodyKinetics)

model_trajectory <- function(pars, times){
  mu <- pars["mu"]*pars["cr"]*pars["mod"]*(pars["priming"]*pars["primed"] + 1*!pars["primed"])
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
      isPrimed <- exposures[i,"primed"]
    
      
      tmpParTab <- parTab[parTab$type %in% c(type,"all"),]
      pars <- tmpParTab$values
      names(pars) <- tmpParTab$names
      
      
      
      mod <- order_tab[order_tab$order == order,"values"]
      
      ind <- sort(c(exposure,strain))
      
      cr <- cr_table[cr_table$exposure==ind[1] & cr_table$strain == ind[2] & cr_table$type == type,"values"]
      
      pars <- c("t_i"=t_i,pars,"mod"=mod,"cr"=cr,y0=0,"primed"=isPrimed)
      
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
