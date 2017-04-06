library(ggplot2)
library(plyr)
#'library(antibodyKinetics)

calculate_x <- function(y, m){
  return(-log(y)/m)
}

model_trajectory_fast <- function(pars, times){
  mu <- pars[4]
  tp <- pars[5]
  dp <- pars[6]
  ts <- pars[7]
  m <- pars[8]
  sigma <- pars[9]
  beta <- pars[10]
  c <- pars[11]
  primed <- pars[12]
  mod <- pars[13]
  x <- pars[14]
  t_i <- pars[15]
  y0 <- 0
  
  cr <- exp(-sigma*x)
  prime_cr <- c*exp(-beta*x)*primed
  
  mu <- mu*cr*mod + prime_cr
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


model_trajectory1 <- function(pars, times){

  x <- pars["x"]
  sigma <- pars["sigma"]
  
  cr <- exp(-sigma*x)
  prime_cr <- pars["c"]*exp(-pars["beta"]*x)*pars["primed"]

  mu <- pars["mu"]*cr*pars["mod"] + prime_cr
 
  dp <- pars["dp"]
  tp <- pars["tp"]
  ts <- pars["ts"]
  m <- pars["m"]
  y0 <- pars["y0"]
  t_i <- pars["t_i"]
  
  #print(paste0("mu: ",pars["mu"]))
  #print(paste0("Combined mu: ", mu))
  
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

## For multiple groups
model_func_groups <- function(parTab, cr_table, order_tab, exposures, strains, times){
  groups <- unique(exposures$group)

  y <- matrix(nrow=length(times)*length(groups),ncol=length(strains))
  colnames(y) <- strains
  y <- data.frame(times=rep(times, length(groups)),y, group=rep(groups,each=length(times)))
  for(group in groups){
    #print(paste0("Group: ", group))
   tmpExposures <- exposures[exposures$group == group,] 
   y[y$group==group,strains] <- model_func(parTab, cr_table, order_tab, tmpExposures, strains, times)
   
  }
  return(y)
}

## For one strain
model_func <- function(parTab, cr_table, order_tab, exposures, strains, times){
  trajectories <- matrix(0, ncol=length(strains),nrow=length(times))
  index <- 1
  new_strains <- character(length(strains))
  for(strain in strains){
   # print(paste0("Strain: ", strain))
    for(i in 1:nrow(exposures)){
      
      t_i <- exposures[i,"values"]
     # print(paste0("Exposure: ",t_i))
      type <- exposures[i,"type"]
      exposure <- exposures[i,"exposure"]
      order <- exposures[i,"order"]
      isPrimed <- exposures[i,"primed"]
      tmpParTab <- parTab[parTab$type %in% c(type,"all","priming"),]
      pars <- tmpParTab$values
      names(pars) <- tmpParTab$names
    
      mod <- order_tab[order_tab$order == order,"values"]
      ind <- sort(c(exposure,strain))
      x <- cr_table[cr_table$exposure==ind[1] & cr_table$strain == ind[2],"values"]
      
      #print(paste0("Exposure strain: ",exposure))
     # print(paste0("cr: ", x))
      
      pars <- c("t_i"=t_i,pars,"mod"=mod,"x"=x,y0=0,"primed"=isPrimed)
      #print(pars)
      y <- model_trajectory1(pars, times)
      #y1 <- model_trajectory_fast(pars[c("S","EA","MAX_TITRE","mu","tp","dp","ts","m","sigma","beta","c","primed","mod","x","t_i")], times)
      #print(y)
      #print(y1)
      #print("")
           #y <- simple_model(pars, times)
      
      trajectories[,index] <-  trajectories[,index] +y
    
    }
    index <- index + 1
    new_strains[index] <- strain
  }
  colnames(trajectories) <- strains
  return(trajectories)
}



create_model_group_func <- function(parTab){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names =="t_i",]
  exposure_indices <- which(parTab$names =="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))
  
  test <- function(pars,times){
    parTab1$values <- pars[parTab_indices]
    cr_table$values <- pars[cr_indices]
    order_tab$values <- pars[order_indices]
    exposures$values <- pars[exposure_indices]
    y <- model_func_groups(parTab1,cr_table,order_tab,exposures,strains,times)
    return(y)
  }
  return(test)  
}


create_model_func <- function(parTab){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names =="t_i",]
  exposure_indices <- which(parTab$names =="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))
  
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
  times <- unique(data[,"time"])
  
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names=="t_i",]
  exposure_indices <- which(parTab$names=="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))

  test <- function(pars){
    parTab1$values <- pars[parTab_indices]
    cr_table$values <- pars[cr_indices]
    order_tab$values <- pars[order_indices]
    exposures$values <- pars[exposure_indices]
    y <- model_func(parTab1,cr_table,order_tab,exposures,strains,times)
    ln <- 0
    for(strain in strains){
      ln <- ln + antibodyKinetics::obs_likelihood(y[,strain],data[,strain], pars)
    }
    
    return(ln)
  }
  return(test)  
}



create_model_groups_post <- function(parTab, data, PRIOR_FUNC=NULL){
  times <- unique(data[,"time"])
  
  groups <- unique(parTab$group)
  groups <- groups[groups != "all"]

  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names=="t_i",]
  exposure_indices <- which(parTab$names=="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))

  test <- function(pars){
    parTab1$values <- pars[parTab_indices]
    cr_table$values <- pars[cr_indices]
    order_tab$values <- pars[order_indices]
    exposures$values <- pars[exposure_indices]
    y <- model_func_groups(parTab1,cr_table,order_tab,exposures,strains,times)
    ln <- 0

    for(group in groups){
      for(strain in strains){
        ln <- ln + antibodyKinetics::obs_likelihood(y[y$group==group,strain],data[data$group == group & data$strain==strain,"value"], pars)
      }
    }
    
    return(ln)
  }
  return(test)  
}






create_model_group_func_fast <- function(parTab){
##########################################################
    ## Firstly, we need to isolate group specific exposures
##########################################################
    ## Get unique groups
    groups <- unique(parTab$group)
    groups <- groups[groups != "all"]
    exposure_indices <- NULL
    exposure_i_lengths <- NULL
    
    ## For each group, isolate the parameter table indices that correspond to exposures for
    ## that group. Save these indices in a contiguous vector, and also store the indices
    ## of THIS vector that each group corresponds to
    for(group in groups){
        tmp <- which(parTab$names == "t_i" & parTab$group == group)
        exposure_indices <- c(exposure_indices, tmp)
        exposure_i_lengths <- c(exposure_i_lengths, length(tmp))
    }
    
    ## The length of this vector is the number of groups plus 1
    exposure_i_lengths <- c(1,cumsum(exposure_i_lengths))
    
#########################################################
    ## Order modifier parameters
#########################################################
    order_indices <- which(parTab$names == "mod")
    
#########################################################
    ## Cross reactivity parameters
#########################################################
    strains <- unique(parTab$strain)
    strains <- strains[!is.na(strains)]
    cr_inds <- NULL
    ## In blocks of length(strains),
    cr_lengths <- rep(length(strains),length(strains))
    cr_lengths <- c(1,cumsum(cr_lengths))
    for(strain1 in strains){
        for(strain2 in strains){
            tmpStrains <- sort(c(strain1,strain2))
            cr_inds <- c(cr_inds,which(parTab$exposure == tmpStrains[1] & parTab$strain == tmpStrains[2] & parTab$names == "x"))
        }
    }
    
#########################################################
    ## Model parameters
#########################################################
    ## Indices of model parameters. 1 = infection, 2 = vaccine, 3 = adjuvant
    param_indices <- which(parTab$index == "parameter")
    par_lengths <- c(length(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")]),
                     length(param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")]),
                     length(param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")]))
    par_type_ind <- c(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")])
    par_lengths <- c(1,cumsum(par_lengths))
    
    
#########################################################
    ## Model parameters
#########################################################
    ## Now we can loop through each group, each strain and
    ## each exposure and get the correct model parameters
    pars <- parTab$values
    exposure_types <- parTab$type
    exposure_strains <- parTab$exposure
    measured_strains <- parTab$strain
    exposure_orders <- parTab$order
    exposure_primes <- parTab$primed
    
    convert_types <- c("all"=0,"infection"=1,"vacc"=2,"adj"=3,"mod"=4,"NA"=5)
    convert_strains <- c("A"=1,"B"=2,"C"=3,"D"=4,"E"=5)
    
    convert_types_back <- c("infection","vacc","adj")
    convert_strains_back <- c("A","B","C","D","E")
    
    exposure_types <- convert_types[exposure_types]
    exposure_strains <- convert_strains[exposure_strains]
    measured_strains <- convert_strains[measured_strains]
    
    strains <- convert_strains[strains]
    
    f <- function(pars, times){
        model_func_groups_fast(pars, times, 
                               groups, strains,
                               exposure_types, exposure_strains, measured_strains, exposure_orders, exposure_primes, 
                               exposure_indices, cr_inds, par_type_ind, order_indices,
                               exposure_i_lengths,  par_lengths, cr_lengths)
    }
    return(f)
    
    
    
}




model_func_groups_fast <- function(pars, times, ## Model solving stuff
                                   groups, strains, ## Admin
                                   exposure_types, exposure_strains, measured_strains, exposure_orders, exposure_primes, ## Vectors from the parameter table
                                   exposure_indices, cr_inds, par_type_ind, order_indices, ## Calculated indices for the various model components
                                   exposure_i_lengths,  par_lengths, cr_lengths ## Length of the various components
                                   ){

    y <- matrix(0,nrow=length(groups)*length(strains),ncol=length(times))
    index <- 1
    for(i in 1:length(groups)){
        group <- groups[i]
                                        #print(paste0("Group: ", group))
        tmp_exposures <- exposure_indices[(exposure_i_lengths[i]:exposure_i_lengths[i+1])]  
        for(strain in strains){
                                        #print(paste0("Strain: ", strain))
            for(j in 1:length(tmp_exposures)){
                t_i <- pars[tmp_exposures[j]]
                                        #print(paste0("Exposure: ",t_i))
                type <- exposure_types[tmp_exposures[j]]
                                        #print(paste0("Exposure type: ", convert_types_back[type]))
                
                order <- exposure_orders[tmp_exposures[j]]
                exposure_strain <- exposure_strains[tmp_exposures[j]]
                                        #print(paste0("Exposure strain: ",convert_strains_back[exposure_strain]))            
                mod <- pars[order_indices[order]]
                isPrimed <- exposure_primes[tmp_exposures[j]]
                cr <- pars[cr_inds[cr_lengths[strain] + exposure_strain]]
                
                fullPars <- pars[par_type_ind[par_lengths[type]:par_lengths[type+1]]]
                fullPars <- c(fullPars, isPrimed, mod, cr,t_i)
                                        #print(fullPars)
                                        #print(model_trajectory_cpp(fullPars,times))
                omg <- model_trajectory_fast(fullPars, times)
                omg1 <- model_trajectory_cpp(fullPars, times)
                                        #print(fullPars)
                                        #print(omg)
                                        #print(omg1)
                                        #print("")
                y[index,] <- y[index,] + model_trajectory_cpp(fullPars,times)
                                        #print("")
            }
            index <- index + 1
                                        #print("")
        }
                                        #print("")
    }
    return(y)
}

create_model_group_func_cpp <- function(parTab){
##########################################################
    ## Firstly, we need to isolate group specific exposures
##########################################################
    ## Get unique groups
    groups <- unique(parTab$group)
    groups <- groups[groups != "all"]
    exposure_indices <- NULL
    exposure_i_lengths <- NULL
    
    ## For each group, isolate the parameter table indices that correspond to exposures for
    ## that group. Save these indices in a contiguous vector, and also store the indices
    ## of THIS vector that each group corresponds to
    for(group in groups){
        tmp <- which(parTab$names == "t_i" & parTab$group == group)
        exposure_indices <- c(exposure_indices, tmp)
        exposure_i_lengths <- c(exposure_i_lengths, length(tmp))
    }
    
    ## The length of this vector is the number of groups plus 1
    exposure_i_lengths <- c(0,cumsum(exposure_i_lengths))
    
#########################################################
    ## Order modifier parameters
#########################################################
    order_indices <- which(parTab$names == "mod")
    
#########################################################
    ## Cross reactivity parameters
#########################################################
    strains <- unique(parTab$strain)
    strains <- strains[!is.na(strains)]
    cr_inds <- NULL
    ## In blocks of length(strains),
    cr_lengths <- rep(length(strains),length(strains))
    cr_lengths <- c(0,cumsum(cr_lengths))
    for(strain1 in strains){
        for(strain2 in strains){
            tmpStrains <- sort(c(strain1,strain2))
            cr_inds <- c(cr_inds,which(parTab$exposure == tmpStrains[1] & parTab$strain == tmpStrains[2] & parTab$names == "x"))
        }
    }
    
#########################################################
    ## Model parameters
#########################################################
    ## Indices of model parameters. 1 = infection, 2 = vaccine, 3 = adjuvant
    param_indices <- which(parTab$index == "parameter")
    par_lengths <- c(length(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")]),
                     length(param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")]),
                     length(param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")]))
    par_type_ind <- c(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")])
    par_lengths <- c(0,cumsum(par_lengths))
    
    
#########################################################
    ## Model parameters
#########################################################
    ## Now we can loop through each group, each strain and
    ## each exposure and get the correct model parameters
    pars <- parTab$values
    exposure_types <- parTab$type
    exposure_strains <- parTab$exposure
    measured_strains <- parTab$strain
    exposure_orders <- parTab$order
    exposure_primes <- parTab$primed
    
    convert_types <- c("all"=0,"infection"=1,"vacc"=2,"adj"=3,"mod"=4,"NA"=5)
    convert_strains <- c("A"=1,"B"=2,"C"=3,"D"=4,"E"=5)
    convert_groups <- c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5)
    
    convert_types_back <- c("infection","vacc","adj")
    convert_strains_back <- c("A","B","C","D","E")
    
    
    exposure_types <- convert_types[exposure_types]
    exposure_strains <- convert_strains[exposure_strains]
    measured_strains <- convert_strains[measured_strains]
    
    strains <- convert_strains[strains]
    groups <- convert_groups[groups]
    
    exposure_indices <- exposure_indices - 1
    cr_inds <- cr_inds - 1
    par_type_ind <- par_type_ind - 1
    order_indices <- order_indices - 1
    exposure_i_lengths <- exposure_i_lengths
    par_lengths <- par_lengths
    cr_lengths <- cr_lengths
    
    f <- function(pars, times){
        model_func_group_cpp(pars, times, groups, strains,
                             exposure_types, exposure_strains, measured_strains, exposure_orders, exposure_primes, 
                             exposure_indices, cr_inds, par_type_ind, order_indices,
                             exposure_i_lengths,  par_lengths, cr_lengths)
    }
    return(f)
    
    
    
}

#' Create posterior calculation function CPP model
#'
#' Creates a function pointer to the cpp version of the ferret model which returns a single likelihood value given a vector of parameters
#' @export
create_posterior_group_func_cpp <- function(parTab, dat, PRIOR_FUNC=NULL){
##########################################################
    ## Firstly, we need to isolate group specific exposures
##########################################################
    ## Get unique groups
    groups <- unique(parTab$group)
    groups <- groups[groups != "all"]
    exposure_indices <- NULL
    exposure_i_lengths <- NULL
    
    ## For each group, isolate the parameter table indices that correspond to exposures for
    ## that group. Save these indices in a contiguous vector, and also store the indices
    ## of THIS vector that each group corresponds to
    for(group in groups){
        tmp <- which(parTab$names == "t_i" & parTab$group == group)
        exposure_indices <- c(exposure_indices, tmp)
        exposure_i_lengths <- c(exposure_i_lengths, length(tmp))
    }
    
    ## The length of this vector is the number of groups plus 1
    exposure_i_lengths <- c(0,cumsum(exposure_i_lengths))
    
#########################################################
    ## Order modifier parameters
#########################################################
    order_indices <- which(parTab$names == "mod")
    
#########################################################
    ## Cross reactivity parameters
#########################################################
    strains <- unique(parTab$strain)
    strains <- strains[!is.na(strains)]
    cr_inds <- NULL
    ## In blocks of length(strains),
    cr_lengths <- rep(length(strains),length(strains))
    cr_lengths <- c(0,cumsum(cr_lengths))
    for(strain1 in strains){
        for(strain2 in strains){
            tmpStrains <- sort(c(strain1,strain2))
            cr_inds <- c(cr_inds,which(parTab$exposure == tmpStrains[1] & parTab$strain == tmpStrains[2] & parTab$names == "x"))
        }
    }
    
#########################################################
    ## Model parameters
#########################################################
    ## Indices of model parameters. 1 = infection, 2 = vaccine, 3 = adjuvant
    param_indices <- which(parTab$index == "parameter")
    par_lengths <- c(length(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")]),
                     length(param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")]),
                     length(param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")]))
    par_type_ind <- c(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")])
    par_lengths <- c(0,cumsum(par_lengths))
    
    
#########################################################
    ## Model parameters
#########################################################
    ## Now we can loop through each group, each strain and
    ## each exposure and get the correct model parameters
    pars <- parTab$values
    exposure_types <- parTab$type
    exposure_strains <- parTab$exposure
    measured_strains <- parTab$strain
    exposure_orders <- parTab$order
    exposure_primes <- parTab$primed
    
    convert_types <- c("all"=0,"infection"=1,"vacc"=2,"adj"=3,"mod"=4,"NA"=5)
    convert_strains <- c("A"=1,"B"=2,"C"=3,"D"=4,"E"=5)
    convert_groups <- c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5)
    
    convert_types_back <- c("infection","vacc","adj")
    convert_strains_back <- c("A","B","C","D","E")
    
    
    exposure_types <- convert_types[exposure_types]
    exposure_strains <- convert_strains[exposure_strains]
    measured_strains <- convert_strains[measured_strains]
    
    strains <- convert_strains[strains]
    groups <- convert_groups[groups]
    
    exposure_indices <- exposure_indices - 1
    cr_inds <- cr_inds - 1
    par_type_ind <- par_type_ind - 1
    order_indices <- order_indices - 1
    exposure_i_lengths <- exposure_i_lengths
    par_lengths <- par_lengths
    cr_lengths <- cr_lengths
    
    times <- dat[1,]
    dat <- dat[2:nrow(dat),]
    
    f <- function(pars){
        posterior_func_group_cpp(pars, times, groups, strains,
                                 exposure_types, exposure_strains, measured_strains, exposure_orders, exposure_primes, 
                                 exposure_indices, cr_inds, par_type_ind, order_indices,
                                 exposure_i_lengths,  par_lengths, cr_lengths, dat)
    }
    return(f)
}



