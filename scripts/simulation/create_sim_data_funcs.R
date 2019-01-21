create_data <- function(runName, parTab_file,
                     exposureTab_file, form,typing,cr,priming,
                     ngroup=5,nstrain=5,nindiv=3,
                     times,wd="~/net/home/ferret/inputs/data_normal",
                     group=NA,strain=NA,
                     normal=TRUE){
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
    exposureTab <- read.csv(exposureTab_file,stringsAsFactors=FALSE)
    if(!is.na(group)){
        exposureTab <- exposureTab[exposureTab$group == group,]
        if(!is.na(strain)){
            exposureTab <- exposureTab[exposureTab$strain == strain,]
        }
        
        parTab <- parTab[parTab$id %in% c(NA,"all",unique(exposureTab$id)) | parTab$names == "x",]
    }
    
    ## Simulate data and save
    individuals <- rep(nindiv,ngroup)
    pars <- parTab[parTab$names %in% c("S","EA","MAX_TITRE"),"values"]
    names(pars) <- c("S","EA","MAX_TITRE")
    print(pars)
                                        #pars <- c("S"=0.79,"EA"=0.2,"MAX_TITRE"=)
    f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=form,typing = typing,cross_reactivity = cr)
    dat <- f(parTab$values, times)
    dat <- floor(dat)
    dat <- apply(dat,2,function(x) rep(x, each=nindiv))
    index <- 1
    difs <- NULL
    for(i in 1:nrow(dat)){
        for(j in 1:ncol(dat)){
            tmp <- dat[i,j]
            dat[i,j] <- add_noise(pars,dat[i,j],normal)
            tmp <- tmp - dat[i,j]
            difs[index] <- tmp
            index <- index + 1
        }
    }
    #dat <- rbind(times, dat)
    rownames(dat) <- NULL

    ## Filenames
    filename <- paste0(wd,"/",runName,"_data.csv")
    if(!is.na(group)){
        filename <- paste0(wd,"/",runName,"_",group,"_data.csv")
    }
    write.csv(dat, filename,row.names=FALSE)
    return(list("data"=dat,"residuals"=difs,"filename"=filename))
}

plot_fit <- function(parTab, chain){
  
  ## Create function to plot titre trajectories, same options as above
  f <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                   form="C",typing = FALSE,cross_reactivity = FALSE)
  ## Times to solve model over
  times1 <- seq(0,max(ferret_titres[1,]),by=0.1)
  
  ## Get MLE parameters
  bestPars <- get_best_pars(chain)
  
  mod <- generate_prediction_intervals(chain, 1000,times1,f,nstrains=1,ngroups=1)
  
  ## Format data for plotting
  ## Just adding labels and changing variable classes
  ## for correct plotting
  meltedDat <- as.data.frame(ferret_titres[2:nrow(ferret_titres),])
  colnames(meltedDat) <- times
  meltedDat <- cbind(meltedDat,expand.grid("indiv"=1:nindiv,"strain"=1:nstrain,"group"=1:ngroup))
  meltedDat <- reshape2::melt(meltedDat,id.vars=c("indiv","strain","group"))
  meltedDat$variable <- as.numeric(as.character(meltedDat$variable))
  meltedDat$group <- as.factor(meltedDat$group)
  meltedDat$strain <- as.factor(meltedDat$strain)
  meltedDat$indiv <- as.factor(meltedDat$indiv)
  
  ## Upper observable bound = 14
  mod[mod$upper > 14,"upper"] <- 14
  mod[mod$lower > 14,"lower"] <- 14
  
  ## Plot trajectory using the MLE parameter set and format this for plotting
  bestTraj <- f(bestPars, times1)
  colnames(bestTraj) <- times1
  bestTraj <- cbind(bestTraj,expand.grid("strain"=1:nstrain,"group"=1:ngroup))
  bestTraj <- reshape2::melt(bestTraj,id.vars=c("strain","group"))
  bestTraj$variable <- as.numeric(as.character(bestTraj$variable))
  bestTraj$group <- as.factor(bestTraj$group)
  bestTraj$strain <- as.factor(bestTraj$strain)
  bestTraj[bestTraj$value > 14,"value"] <- 14
  
  ## Plot model fit over data
  p <- ggplot() + 
    geom_ribbon(data = mod, aes(x=time,ymax=upper,ymin=lower),alpha=0.4)+
    geom_line(data=bestTraj,aes(x=variable,y=value))+
    geom_line(data=meltedDat,aes(x=variable,y=value,group=indiv,linetype=indiv))+
    scale_linetype_manual(values=c("twodash","dotted","dashed"))+
   # scale_linetype_manual(values=1:10)+
    geom_point(data = meltedDat,aes(x=variable,y=value,shape=indiv)) +
    scale_y_continuous(limits=c(0,14),breaks=seq(0,14,by=2),expand=c(0,0))+
    ylab("HI titre") +
    xlab("Time post infection") +
    #facet_wrap(~group,ncol=1) +
    theme(legend.position = "none") +
    theme_bw()
  return(p)
}

