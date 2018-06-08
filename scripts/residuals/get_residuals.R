get_residuals <- function(runName, parTab_file, exposureTab_file, form, typing, cr, 
                          dat_file, n, mcmcPars, univariate_runs, biphasic_waning=TRUE,
                          y0_mod=TRUE,wd){
  wd <- paste0(wd,runName)
  print(wd)
  setwd(wd)
  
  chains <- NULL
  index <- 1
  for(i in 1:5){
    if(runName %in% univariate_runs){
      filename <- paste0(runName,"_",i,"_univariate_chain.csv")
    } else {
      filename <- paste0(runName,"_",i,"_multivariate_chain.csv")
    }
    if(file.exists(filename)){
      #print(index)
      chain <- read.csv(filename)
      chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
      chains[[index]] <- chain
      index <- index + 1
    }
  }
  
  chain <- do.call("rbind",chains)
  
  exposureTab <- read.csv(paste0("~/net/home/ferret/",exposureTab_file),stringsAsFactors=FALSE)
  parTab <- read.csv(paste0("~/net/home/ferret/",parTab_file),stringsAsFactors=FALSE)
  parTab[parTab$names == "MAX_TITRE","values"] <- 13
  parTab[parTab$names == "tp","values"] <- 12
  if(!biphasic_waning){
    parTab[parTab$names %in% c("ts","dp"),"fixed"] <- 1
    parTab[parTab$names %in% c("ts","dp"),"values"] <- 0
  }
  dat <- read.csv(paste0("~/net/home/ferret/",dat_file))
  dat <- as.data.frame(dat)
  indices <- expand.grid(1:3,1:5,1:5)
  colnames(indices) <- c("indiv","strain","group")
  dat <- cbind(dat,indices)
  
  f <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                   form=form,typing = typing,cross_reactivity = cr)
  samples <- sample(1:nrow(chain),n,replace=FALSE)
  res <- NULL
  all_calcs <- NULL
  index <- 1
  for(j in samples){
    pars <- as.numeric(chain[j,2:(ncol(chain)-1)])
    y <- f(pars,times)
    y <- apply(y,2,function(x) rep(x, each=3))
    y[y > pars[4]] <- pars[4]
    y <- floor(y)
    y <- as.data.frame(y)
    y <- cbind(y,indices)
    y$samp <- index
    y$runName <- runName
    y <- reshape2::melt(y,id.vars=c("runName","samp","indiv","strain","group"))
    res <- rbind(y,res)
    #res[index] <- sum((y - dat[2:nrow(dat),])^2)
    index <- index + 1
  }
  
  bestPars <- get_best_pars(chain)
  y <- f(pars,times)
  y <- apply(y,2,function(x) rep(x, each=3))
  y[y > pars[4]] <- pars[4]
  y <- floor(y)
  y <- as.data.frame(y)
  y <- cbind(y,indices)
  y$samp <- "MLE"
  y$runName <- runName
  y <- reshape2::melt(y,id.vars=c("runName","samp","indiv","strain","group"))
  res <- rbind(y,res)
  res$variable <- times[as.numeric(res$variable)]
  dat <- reshape2::melt(dat,id.vars=c("indiv","strain","group"))
  dat$variable <- times[as.numeric(dat$variable)]
  return(list(res,dat))
}