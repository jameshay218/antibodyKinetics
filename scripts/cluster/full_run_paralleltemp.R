full_run_paralleltemp <- function(runName,
                       runID,
                       chainNo,
                       parTab_file,
                       exposureTab_file,
                       dat_file=NA,
                       y0mod_ver=-1,
                       fixed_S=TRUE,
                       typing=TRUE,
                       ngroup=5,nstrain=5,nindiv=3,
                       times,
                       mcmcPars,
                       sim=FALSE,
                       run_times="protocol"
                       ){

    if(runID == "3days" | run_times == "3days") times <- seq(0,72,by=3)
    if(runID == "protocol" | run_times == "protocol") times <- c(0,21,37,49,70)
    if(runID == "weekly" | run_times == "weekly") times <- seq(0,70,by=7)

    
  
    ## How many temperature chains to run?
    n_temperatures <- length(mcmcPars[["temperature"]])

    ## Read in correct parameter and exposure tables for this model run
    ## Generate list of options for this run
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
    options <- convert_runName_to_options(runName)
    parTab <- parTab_modification(parTab, options, fixed_S)
    exposureTab <- read.csv(exposureTab_file,stringsAsFactors=FALSE)
    
    ## Enumerate out number of individuals per group
    individuals <- rep(nindiv,ngroup)
    
    f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=options$form,typing = typing,cross_reactivity = options$cr)
    ## Change save location depending on simulated or real data
    if(sim){
        if(!fixed_S){
            wd <- paste0(getwd(),"/outputs_sim/",runID,"_",runName)
        } else {
            wd <- paste0(getwd(),"/outputs_sim_fixedS/",runID,"_",runName)
        }
    } else {
        if(!fixed_S){
            wd <- paste0(getwd(),"/outputs_real/",runID,"_",runName)
        } else {
            wd <- paste0(getwd(),"/outputs_real_fixedS/",runID,"_",runName)
        }
    }
    
    if(!dir.exists(wd)) dir.create(wd,recursive=TRUE)
    
    ## Filenames
    filename <- paste0(wd,"/",runID,"_",runName,"_",chainNo,"_data.csv")
    chain_file <- paste0(wd,"/",runID,"_",runName,"_",chainNo)
    
    ## If no data file is provided, create some data from the given
    ## parTab file
    pars <- parTab[parTab$names %in% c("S","EA","MAX_TITRE"),"values"]
    
#####
    ## Changing upper bound on m for monophasic waning
#####
    parTab[parTab$names == "m","upper_bound"] <- 12
    if(y0mod_ver != -1){
        if(y0mod_ver == 1){
            parTab[parTab$names == "y0_mod","lower_bound"] <- 0
        } else {
            parTab[parTab$names == "y0_mod","upper_bound"] <- 0
        }        
    }

    
    names(pars) <- c("S","EA","MAX_TITRE")

    ## Process either real data, or simulate data if no data file is
    ## provided
    if(is.na(dat_file)){
        dat <- f(parTab$values, times)
        dat <- floor(dat)
        dat <- apply(dat,2,function(x) rep(x, each=nindiv))
        ## Add some noise to each 
        for(i in 1:nrow(dat)){
            for(j in 1:ncol(dat)){
                dat[i,j] <- add_noise(pars,dat[i,j])
            }
        }
        dat <- rbind(times, dat)
        rownames(dat) <- NULL
        write.csv(dat, filename)
    } else {
        dat <- read.csv(dat_file)
        dat <- as.matrix(rbind(times, dat))
        rownames(dat) <- NULL
    }
    
    startTab <- parTab
      ## Generate random starting conditions for the chain
    ## Note that a list of starting tables is created,
    ## one for each temperature chain
    if(n_temperatures > 1){
        startTab <- rep(list(startTab),n_temperatures)
        for(j in 1:length(startTab)){
            for(i in which(parTab$fixed == 0)){
                startTab[[j]][i,"values"] <- runif(1,startTab[[j]][i,"lower_bound"],startTab[[j]][i,"upper_bound"])
            }
        }
    }

    ## Run the parallel tempering MCMC chain
    run_1 <- lazymcmc::run_MCMC(parTab=startTab,data=dat, mcmcPars=mcmcPars, 
                                filename=chain_file,CREATE_POSTERIOR_FUNC = create_model_group_func_cpp,
                                mvrPars=NULL,PRIOR_FUNC=NULL, 
                                version="posterior",form=options$form,
                                individuals=individuals,exposureTab=exposureTab,
                                cross_reactivity=options$cr,typing=typing)
    chain <- read.csv(run_1$file)
    ## Plot and save MCMC traces
    pdf(paste0(chain_file, ".pdf"))
    plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],c(which(parTab$fixed == 0) + 1,ncol(chain))]))
    dev.off()
    
    times1 <- seq(0,100,by=0.1)
    ## Plot and save inferred model trajectories and 95% prediction intervals over data
    bestPars <- get_best_pars(chain)
    
    mod <- generate_prediction_intervals(chain=chain, samp_no=200,ts=seq(0,100,by=1),MODEL_FUNC=f,nstrains=nstrain,ngroups=ngroup)[[1]]
    
    meltedDat <- as.data.frame(dat[2:nrow(dat),])
    colnames(meltedDat) <- times
    meltedDat <- cbind(meltedDat,expand.grid("indiv"=1:nindiv,"strain"=1:nstrain,"group"=1:ngroup))
    meltedDat <- reshape2::melt(meltedDat,id.vars=c("indiv","strain","group"))
    meltedDat$variable <- as.numeric(as.character(meltedDat$variable))
    meltedDat$group <- as.factor(meltedDat$group)
    meltedDat$strain <- as.factor(meltedDat$strain)
    meltedDat$indiv <- as.factor(meltedDat$indiv)
    mod[mod$upper > 12,"upper"] <- 12
    mod[mod$lower > 12,"lower"] <- 12
    
    bestTraj <- f(bestPars, seq(0,100,by=1))
    colnames(bestTraj) <- seq(0,100,by=1)
    bestTraj <- cbind(bestTraj,expand.grid("strain"=1:nstrain,"group"=1:ngroup))
    bestTraj <- reshape2::melt(bestTraj,id.vars=c("strain","group"))
    bestTraj$variable <- as.numeric(as.character(bestTraj$variable))
    bestTraj$group <- as.factor(bestTraj$group)
    bestTraj$strain <- as.factor(bestTraj$strain)
    bestTraj[bestTraj$value > 12,"value"] <- 12
    
    p <- ggplot() + 
        geom_ribbon(data = mod, aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4)+
        geom_line(data=bestTraj,aes(x=variable,y=value,col=strain))+
        geom_point(data = meltedDat,aes(x=variable,y=value,col=strain),position=position_jitter(w=0.5,h=0.5)) +
        facet_wrap(~group,ncol=1) +
        theme_bw()
    
    pdf(paste0(chain_file, "_plot.pdf"))
    print(p)
    dev.off()
}
