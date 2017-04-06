#' Generate prediction intervals from MCMC chain
#'
#' Produces a data frame of lower 2.5%, upper 97.5% and median prediction intervals for a given model function.
#' @param chain the MCMC chain to sample from
#' @param samp_no the number of random samples to take from the MCMC chain
#' @param ts the times over which to solve the model
#' @param MODEL_FUNCTION pointer to the model solving function that takes two arguments: pars and ts. This should return a matrix of trajectories with times across columns
#' @param nstrains the number of strains for which measurements are available
#' @param ngroups the number of exposure groups
#' @export
generate_prediction_intervals <- function(chain, samp_no=1000,ts, MODEL_FUNCTION,nstrains=3,ngroups=5){
    ## Take random samples from the MCMC chain 
    samps <- sample(nrow(chain),samp_no)

    ## The result of the model solving procedure will be a matrix of titres, with columns as time
    ## and rows as groups/strains. As such, we need a way to store samp_no*(ngroups*nstrains) trajectories
    ## of length length(ts). Best way off the top of my head - have a data frame for each
    ## group and strain combo, with samples as rows.
    dats <- replicate(ngroups,
                      replicate(
                          nstrains,matrix(nrow=samp_no,ncol=length(ts)),
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

        ## Save this trajectory
        index <- 1
        for(x in 1:ngroups){
            for(j in 1:nstrains){
                dats[[x]][[j]][i,] <- y[index,]
                index <- index + 1
            }
        }
    }
    
    allQuantiles <- NULL
    for(x in 1:ngroups){
        tmpQuantiles <- NULL
        for(j in 1:nstrains){
            quantiles <- t(apply(dats[[x]][[j]], 2, function(x) quantile(x,c(0.025,0.5,0.975))))
            quantiles <- data.frame(time=ts,quantiles)
            colnames(quantiles) <- c("time","lower","median","upper")
            quantiles$strain <- j
            tmpQuantiles[[j]] <- quantiles
        }
        tmpQuantiles <- do.call("rbind",tmpQuantiles)
        tmpQuantiles$group <- x
        allQuantiles[[x]] <- tmpQuantiles
    }
    allQuantiles <- do.call("rbind",allQuantiles)
    allQuantiles$group <- as.factor(allQuantiles$group)
    allQuantiles$strain <- as.factor(allQuantiles$strain)
    return(allQuantiles)
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
