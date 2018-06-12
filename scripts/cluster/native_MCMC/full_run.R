full_run <- function(runName,
                     runID,
                     chainNo,
                     parTab_file,
                     exposureTab_file,
                     dat_file=NA,
                     fixed_S=TRUE,
                     typing=TRUE,
                     ngroup=5,nstrain=5,nindiv=3,
                     times,
                     mcmcPars1,
                     mcmcPars2,
                     sim=FALSE
                     ){
    ## Read in correct parameter and exposure tables for this model run
    ## Generate list of options for this run
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
    options <- convert_runName_to_options(runName)
    parTab <- parTab_modification(parTab, options, fixed_S)
    exposureTab <- read.csv(exposureTab_file,stringsAsFactors=FALSE)

    ## Enumerate out number of individuals per group
    individuals <- rep(nindiv,ngroup)
    f <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                     form=options$form,typing = typing,cross_reactivity = options$cr)

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
    chain1_file <- paste0(wd,"/",runID,"_",runName,"_",chainNo,"_univariate")
    chain2_file <- paste0(wd,"/",runID,"_",runName,"_",chainNo,"_multivariate")

    ## If no data file is provided, create some data from the given
    ## parTab file
    pars <- parTab[parTab$names %in% c("S","EA","MAX_TITRE"),"values"]

    #####
    ## Changing upper bound on m for monophasic waning
    #####
    parTab[parTab$names == "m","upper_bound"] <- 12

    ## Process either real data, or simulate data if no data file is
    ## provided
    names(pars) <- c("S","EA","MAX_TITRE")
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

    ## Generate random starting conditions for the chain
    startTab <- parTab
    for(i in which(parTab$fixed == 0)){
        startTab[i,"values"] <- runif(1,startTab[i,"lower_bound"],startTab[i,"upper_bound"])
    }

    ## Run the first chain using the univariate sampler
    run_1 <- antibodyKinetics::run_MCMC(startTab,dat, mcmcPars1, chain1_file,create_model_group_func_cpp,
                      NULL,NULL, version="posterior",form=options$form,
                      individuals=individuals,exposureTab=exposureTab,
                      cross_reactivity=options$cr,typing=typing)

    chain <- data.table::fread(run_1$file,data.table=FALSE)
    bestPars <- get_best_pars(chain)

    chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]

    ## Use this initial run to generate a covariance matrix for the multivariate proposal
    ## chain
    covMat <- cov(chain)
    mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

    parTab1 <- parTab
    parTab1$values <- bestPars

    run_2 <- antibodyKinetics::run_MCMC(parTab1,dat, mcmcPars2, chain2_file,create_model_group_func_cpp,
                      mvrPars,NULL, version="posterior",form=options$form,
                      individuals=individuals,exposureTab=exposureTab,
                      cross_reactivity=options$cr,typing=typing)
    chain1 <- read.csv(run_2$file)

    ## Plot and save MCMC traces
    pdf(paste0(chain2_file, ".pdf"))
    plot(coda::as.mcmc(chain1[chain1$sampno > mcmcPars1["adaptive_period"],c(which(parTab1$fixed == 0) + 1,ncol(chain1))]))
    dev.off()


    ## Plot and save inferred model trajectories and 95% prediction intervals over data
    times1 <- seq(0,100,by=0.1)
    
    bestPars <- get_best_pars(chain1)
    
    mod <- generate_prediction_intervals(chain1, 200,seq(0,100,by=1),f,nstrains=nstrain,ngroups=ngroup)
    
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

    ##y <- f(parTab$values, seq(0,100,by=0.1))
    ##y <- rbind(times1,y)
    ##realTraj <- as.data.frame(y[2:nrow(y),])
    ##realTraj[realTraj > pars["MAX_TITRE"]] <- pars["MAX_TITRE"]
    ##colnames(realTraj) <- times1
    ##realTraj <- cbind(realTraj,expand.grid("strain"=1:nstrain,"group"=1:ngroup))
    ##realTraj <- reshape2::melt(realTraj,id.vars=c("strain","group"))
    ##realTraj$variable <- as.numeric(as.character(realTraj$variable))
    ##realTraj$group <- as.factor(realTraj$group)
    ##realTraj$strain <- as.factor(realTraj$strain)
    
    p <- ggplot() + 
        geom_ribbon(data = mod, aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.4)+
        geom_line(data=bestTraj,aes(x=variable,y=value,col=strain))+
        geom_point(data = meltedDat,aes(x=variable,y=value,col=strain),position=position_jitter(w=0.5,h=0.5)) +
        facet_wrap(~group,ncol=1) +
        theme_bw()
    
    pdf(paste0(chain2_file, "_plot.pdf"))
    print(p)
    dev.off()
}

