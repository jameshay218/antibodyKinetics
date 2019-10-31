#' Generate prediction intervals from MCMC chain
#'
#' Produces a data frame of lower 2.5%, upper 97.5% and median prediction intervals for a given model function.
#' @param chain the MCMC chain to sample from
#' @param samp_no the number of random samples to take from the MCMC chain
#' @param ts the times over which to solve the model
#' @param ts_obs times at which observations were made
#' @param MODEL_FUNCTION pointer to the model solving function that takes two arguments: pars and ts. This should return a matrix of trajectories with times across columns
#' @param nstrains the number of strains for which measurements are available
#' @param ngroups the number of exposure groups
#' @param calc_obs if TRUE, also calculates predicted observations
#' @param ci_range defaults to 95%, but can be changed to give say 50% CI
#' @export
generate_prediction_intervals <- function(chain, samp_no=1000,ts=seq(0,70,by=1), ts_obs=ts, MODEL_FUNCTION,nstrains=3,ngroups=5, calc_obs=TRUE, ci_range=0.95){
    ## Take random samples from the MCMC chain
    samp_max <- max(nrow(chain), samp_no)
    samps <- sample(nrow(chain),samp_max, replace=FALSE)

    ci_lower <- (1-ci_range)/2
    ci_upper <- ci_range + ci_lower
    
    ## The result of the model solving procedure will be a matrix of titres, with columns as time
    ## and rows as groups/strains. As such, we need a way to store samp_no*(ngroups*nstrains) trajectories
    ## of length length(ts). Best way off the top of my head - have a data frame for each
    ## group and strain combo, with samples as rows.
    dats <- replicate(ngroups,
                      replicate(
                          nstrains,matrix(nrow=samp_no,ncol=length(ts)),
                          simplify=FALSE),
                      simplify=FALSE)
    
    observations <- replicate(ngroups,
                      replicate(
                          nstrains,matrix(nrow=samp_no,ncol=length(ts_obs)),
                          simplify=FALSE),
                      simplify=FALSE)
    
    ## For each sample
    for(i in 1:samp_no){
        samp <- samps[i]
        ## Get parameters
        pars <- as.numeric(chain[samp,!(colnames(chain) %in% c("sampno","lnlike"))])
        names(pars) <- colnames(chain[,!(colnames(chain) %in% c("sampno","lnlike"))])

        ## Solve model
        y <- MODEL_FUNCTION(pars, ts)
        if(calc_obs){
            y_obs <- MODEL_FUNCTION(pars, ts_obs)
            ## Save this trajectory
            index <- 1
            
            for(x in 1:ngroups){
                for(j in 1:nstrains){
                    dats[[x]][[j]][i,] <- y[index,]
                    #obs <- floor(y_obs[index,])
                                        #print(obs)
                    obs <- add_noise_serosolver(y_obs[index,], pars)
                    #obs <- sapply(obs, add_noise, pars=pars, normal=TRUE)
                                        #print(obs)
                    observations[[x]][[j]][i,] <- obs
                    index <- index + 1
                }
            }
        }
    }
    allQuantiles <- NULL
    allQuantiles_observations <- NULL
    for(x in 1:ngroups){
        tmpQuantiles <- NULL
        tmpQuantiles_observations <- NULL
        for(j in 1:nstrains){
            quantiles <- t(apply(dats[[x]][[j]], 2, function(x) quantile(x,c(ci_lower,0.5,ci_upper))))
            quantiles <- data.frame(time=ts,quantiles)
            colnames(quantiles) <- c("time","lower","median","upper")
            quantiles$strain <- j
            tmpQuantiles[[j]] <- quantiles

            quantiles_obs <- t(apply(observations[[x]][[j]], 2, function(x) quantile(x,c(ci_lower,0.5,ci_upper))))
            quantiles_obs <- data.frame(time=ts_obs,quantiles_obs)
            colnames(quantiles_obs) <- c("time","lower","median","upper")
            quantiles_obs$strain <- j
            tmpQuantiles_observations[[j]] <- quantiles_obs

            
        }
        tmpQuantiles <- do.call("rbind",tmpQuantiles)
        tmpQuantiles$group <- x
        allQuantiles[[x]] <- tmpQuantiles
        
        tmpQuantiles_observations<- do.call("rbind",tmpQuantiles_observations)
        tmpQuantiles_observations$group <- x
        allQuantiles_observations[[x]] <- tmpQuantiles_observations
    }
    allQuantiles <- do.call("rbind",allQuantiles)
    allQuantiles$group <- as.factor(allQuantiles$group)
    allQuantiles$strain <- as.factor(allQuantiles$strain)


    allQuantiles_observations <- do.call("rbind",allQuantiles_observations)
    allQuantiles_observations$group <- as.factor(allQuantiles_observations$group)
    allQuantiles_observations$strain <- as.factor(allQuantiles_observations$strain)
    
    return(list(allQuantiles, allQuantiles_observations))
}




#' @export
generate_prediction_intervals_original <- function(chain, samp_no=1000,ts){
  samps <- sample(nrow(chain),samp_no)
  
  dats <- matrix(nrow=samp_no,ncol=length(ts))
  
  for(i in 1:nrow(dats)){
    samp <- samps[i]
    pars <- as.numeric(chain[samp,!(colnames(chain) %in% c("sampno","lnlike"))])
    names(pars) <- colnames(chain[,!(colnames(chain) %in% c("sampno","lnlike"))])
    y <- predict_titres(pars[!(names(pars) %in% c("S","EA","MAX_TITRE"))],ts)
    dats[i,] <- y[,2]
  }
  quantiles <- t(apply(dats, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  quantiles <- data.frame(time=ts,quantiles)
  colnames(quantiles) <- c("time","lower","median","upper")
  
  return(quantiles)
}

#' Find best MCMC chain parameters
#'
#' Given an MCMC chain, finds the row with the highest log likelihood. There must be a column named "lnlike"
#' which is the log likelihood. Will return the vector of all parameters other than "sampno","lnlike" and "strain"
#' @param chain the MCMC chain
#' @return a named vector of parameters
#' @export
get_best_pars <- function(chain){
    pars <- as.numeric(chain[which.max(chain[,"lnlike"]),])
    names(pars) <- colnames(chain)
    pars <- pars[!(names(pars) %in% c("sampno","lnlike","strain"))]
    return(pars)    
}


#' @export
base_plot <- function(chains, data, samps=1000,infection_times=NULL){
    ts <- seq(min(data$time),max(data$time),by=1)
    quants <- plyr::ddply(chains,~strain,
                          function(x) generate_prediction_intervals_original(x[,colnames(x) != "strain"],samps,ts))
    best_pars <-  plyr::ddply(chains,~strain,get_best_pars)
    best_y <- plyr::ddply(best_pars, ~strain,
                          function(x){
                              pars <- as.numeric(x)
                              names(pars) <- colnames(x)
                              as.data.frame(predict_titres(pars[!(names(pars) %in% c("S","EA","MAX_TITRE","strain"))],ts))
                          })
    colnames(best_y) <- c("strain","time","value")


    xscale <- c(0,21,37,49,70, infection_times$time)
    xlabels <- c("0","21","37","49","70",paste("\n\n",infection_times$infection,sep=""))
    xlabel_colours <- c(rep("gray20",5),rep("red",nrow(infection_times)))
    xlabel_sizes <- c(rep(14,5),rep(10,nrow(infection_times)))
    
    p <- ggplot() +
        geom_point(data=data,aes(x=time,y=value,col=strain),size=4,position=position_jitter(w=0.5,h=0.1)) +
        geom_line(data=best_y,aes(x=time,y=value,col=strain),size=0.8)+
        geom_ribbon(data=quants,aes(x=time,ymax=upper,ymin=lower,fill=strain),alpha=0.2) +
        coord_cartesian(ylim=c(0,best_pars[1,"MAX_TITRE"]))+
        scale_fill_brewer(palette="Dark2") +
        scale_colour_brewer(palette="Dark2") +
        xlab("Time (days)") +
        ylab("Log Titre") +
        theme(
            legend.justification=c(0,1),
            legend.position=c(0,1),
            text=element_text(size=16,colour="gray20"),
            plot.title=element_text(size=28),
            legend.text=element_text(size=14,colour="gray20"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line=element_line(colour="gray20"),
            axis.line.x = element_line(colour = "gray20"),
            axis.line.y=element_line(colour="gray20"),
            axis.text.x=element_text(colour=xlabel_colours,size=xlabel_sizes),
            plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),
            panel.background=element_blank(),
            axis.text.y=element_text(colour="gray20",size=14))+
        scale_x_continuous(breaks=xscale,labels=xlabels)
    
    
    if(!is.null(infection_times)){
        for(i in 1:nrow(infection_times)){
            p <- p + geom_vline(xintercept=infection_times[i,"time"],colour="red",linetype="longdash")
        }
    }
    return(p)    
}

#' Plot single trajectory model fit
#'
#' 
#' @export
plot_single_fit <- function(parTab, exposureTab,
                            ferret_titres, times,
                            chain, options,
                            add_sim_obs=FALSE,
                            n=1000,
                            nindiv=3,
                            nstrain=1,
                            ngroup=1){
    ## Create function to plot titre trajectories, same options as above
    f <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                     form=options$form,typing = TRUE,cross_reactivity = TRUE)
    ## Times to solve model over
    times1 <- seq(0,max(ferret_titres[1,]),by=0.1)
    
    ## Get MLE parameters
    bestPars <- get_best_pars(chain)  
    pred_intervals <- generate_prediction_intervals(chain, n,times1,times,f,nstrains=1,ngroups=1)
    mod <- pred_intervals[[1]]
    sim_obs <- pred_intervals[[2]]
    
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
        geom_ribbon(data = mod, aes(x=time,ymax=upper,ymin=lower),fill="grey60",col="grey60")+
        geom_line(data=bestTraj,aes(x=variable,y=value))+
        geom_line(data=meltedDat,aes(x=variable,y=value,group=indiv,linetype=indiv),col="royalblue2")#
    if(nindiv <= 3){
        p <- p + scale_linetype_manual(values=c("twodash","dotted","dashed"))
    }
  
    p <- p +
        geom_point(data = meltedDat,aes(x=variable,y=value), size=2, shape=4, col="blue")
    if(add_sim_obs){
        p <- p + geom_point(data=sim_obs, aes(x=time,y=median),stat="identity",
                            size=2,col="black",shape=18) +
            geom_errorbar(data=sim_obs,aes(x=time,ymin=lower,ymax=upper),
                          col="black",stat="identity",
                          width=1.5)
    }
    p <- p + scale_y_continuous(limits=c(0,14),breaks=seq(0,14,by=2),expand=c(0,0))+
        scale_x_continuous(expand=c(0.02,0.02),limits=c(-2,72),breaks=seq(0,70,by=10))+
        ylab("log titre") +
        xlab("Time post infection (days)") +
        theme_bw() +
        theme(legend.position = "none",
              strip.background = element_blank(),
              strip.text=element_blank(),
              legend.title=element_blank(),
              legend.direction = "vertical",
              axis.text=element_text(family="Arial",colour="gray20"),
              axis.text.x=element_text(size=8),
              axis.text.y=element_text(size=8),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10),
              legend.text=element_text(size=8),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              axis.line=element_line(colour="gray20"),
              axis.line.x = element_line(colour = "gray20"),
              axis.line.y=element_line(colour="gray20"),
              plot.margin = unit(c(1, 0, 0, 0), "cm"),
              panel.spacing=unit(1,"lines"),
              panel.background=element_blank())
    return(p)
}
