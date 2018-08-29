type_names <- c("all"="All","vacc"="TIV","adj"="TIV + adjuvant","infection"="Infection",
                "vacc1"="TIV 1","vacc2"="TIV 2","adj1"="TIV 1 + \nadjuvant","adj2"="TIV 2 + \nadjuvant",
                "infection1"="Infection 1","infection2"="Infection 2","Priming"="Priming","priming"="Priming")
ordered_names <- c("All","TIV","TIV + adjuvant","Infection",
                   "TIV 1","TIV 2","TIV 1 + \nadjuvant","TIV 2 + \nadjuvant",
                   "Infection 1","Infection 2","Priming")
ordered_names2 <- c("All","TIV","TIV + adjuvant","Infection",
                   "TIV 1","TIV 2","TIV 1 + adjuvant","TIV 2 + adjuvant",
                   "Infection 1","Infection 2","Priming")
den_plot <- function(chain, parName, parTab, options, ylab, ymax, plot_intervals=TRUE, skip_pars=NULL,yupper=NULL,ymin=0, add_priming_blank=TRUE){
    print(parName)
    if(is.null(yupper)) yupper <- ymax
    print(yupper)
    if(parName %in% c("ts","dp") & options$monophasic_waning ==TRUE) return(NULL)
    if(parName == "mod" & !options$antigenic_seniority) return(NULL)
    if(parName == "y0_mod" & !options$y0_mod) return(NULL)
    
    if(parName == "mu" & options$priming == TRUE) parName <- c(parName, "c")
    if(parName == "sigma" & options$priming == TRUE) parName <- c(parName, "beta")
    
    tmp_chain <- chain[,which(colnames(chain) %in% parName), drop=FALSE]
    if(parName == "adjmu"){
        tmp_chain1 <- chain[,which(colnames(chain) == "mu"), drop=FALSE]
        tmp_chain2 <- chain[,which(colnames(chain) == "dp"), drop=FALSE]
        for(i in 1:ncol(tmp_chain1)) tmp_chain1[,i] <- tmp_chain1[,i] * (1-tmp_chain2[,i])
        tmp_chain <- tmp_chain1
        parName <- "mu"
    }
    
    melted_chain <- reshape2::melt(tmp_chain)
    colnames(melted_chain) <- c("parName","value")
    melted_chain$parName <- as.character(melted_chain$parName)
    
    types_table <- unique(parTab[parTab$names %in% parName,c("names","type")])

    if(parName == "mod"){
        types_table <- rep(types_table, 4)
        types_table$names <- c("mod","mod.1","mod.2","mod.3")
    }

   
    types_table$names <- colnames(tmp_chain)

    
    if(options$priming) types_table[types_table$names %in% c("beta","c"),"type"] <- "Priming"
    
    colnames(types_table) <- c("parName","type")
    types_table$type <- type_names[types_table$type]
    
    melted_chain <- merge(melted_chain,types_table,by="parName")
    
    quantiles <- plyr::ddply(melted_chain, c("parName","type"), function(x) signif(quantile(x$value,c(0.025,0.5,0.975)),3))
    means <- plyr::ddply(melted_chain, c("parName","type"), function(x) signif(mean(x$value),3))
    modes <- plyr::ddply(melted_chain, c("parName","type"), function(x) signif(estimate_mode(x$value),3))
    
    res <- merge(quantiles, means, id.vars=c("parName","type"))
    res <- merge(res, modes, by=c("parName","type")) 
    
    colnames(res) <- c("parName","type","lower","median","upper","mean","mode")
    
    ## Code to reorder the exposure type factor levels
    present_names <- ordered_names[ordered_names %in% as.character(unique(melted_chain$type))]
    if(add_priming_blank) present_names <- c(present_names,"Priming")
    melted_chain$type = factor(melted_chain$type,present_names)
    melted_chain[melted_chain$parName %in% skip_pars,"value"] <- NA

    
    use_var <- "type"
    if(parName == "mod"){
      use_var <- "parName"
      convert <- c("mod"="ρ1","mod.1"="ρ2","mod.2"="ρ3","mod.3"="ρ4")
      melted_chain$parName <- convert[melted_chain$parName]
    }
    if(plot_intervals){
        p <- ggplot(melted_chain) + 
            geom_violin(aes_string(x=use_var,y="value"),fill="grey",draw_quantiles=c(0.025,0.5,0.975),trim=TRUE,scale="width") + 
            coord_cartesian(ylim=c(ymin,ymax)) + xlab("") + ylab(ylab) +
            geom_hline(yintercept = c(ymin,yupper),linetype="dashed", alpha=0.3) +
             scale_x_discrete(drop=FALSE) +
            theme(axis.text.x = element_text(size=8,family="Arial",colour="black"),
                  axis.text.y = element_text(size=8,family="Arial",colour="black"),
                  axis.title.y=element_text(size=10,family="Arial",colour="black"))  + theme_classic()
    } else {
        p <- ggplot(melted_chain) + 
            geom_violin(aes_string(x=use_var,y="value"),fill="grey",draw_quantiles=c(0.5),trim=TRUE,scale="width") + 
            coord_cartesian(ylim=c(ymin,ymax)) + xlab("") + ylab(ylab) +
            geom_hline(yintercept = c(ymin,yupper),linetype="dashed", alpha=0.3) +
            scale_x_discrete(drop=FALSE) +
            theme(axis.text.x = element_text(size=8,family="Arial",colour="black"),
                  axis.text.y = element_text(size=8,family="Arial",colour="black"),
                  axis.title.y=element_text(size=10,family="Arial",colour="black"))  + theme_classic()
    }
    return(list(res, p))
}


generate_cr_profiles <- function(chain, options, parTab, x_max){
    mu_chain <- chain[,c(1,which(colnames(chain) == "mu"),ncol(chain)),drop=FALSE]
    types_table_mu<- unique(parTab[parTab$names =="mu",c("names","type")])
    types_table_mu$names <- colnames(mu_chain)[2:(ncol(mu_chain)-1)]
    mu_chain <- reshape2::melt(mu_chain, c("sampno","chain"))
    colnames(mu_chain) <- c("sampno","chain","names","mu")
    mu_chain <- join(mu_chain,types_table_mu, by="names")
    mu_chain <- mu_chain[,c("sampno","chain","mu","type")]
    
    sigma_chain <- chain[,c(1,which(colnames(chain) == "sigma"),ncol(chain)),drop=FALSE]
    types_table_sigma <- unique(parTab[parTab$names =="sigma",c("names","type")])
    types_table_sigma$names <- colnames(sigma_chain)[2:(ncol(sigma_chain)-1)]
    sigma_chain <- reshape2::melt(sigma_chain, c("sampno","chain"))
    colnames(sigma_chain) <- c("sampno","chain","names","sigma")
    sigma_chain <- join(sigma_chain,types_table_sigma, by="names")
    sigma_chain <- sigma_chain[,c("sampno","chain","sigma","type")]

    if(!options$typed_cr){
        sigma_chain <- rep(sigma_chain,nrow(types_table_mu))
        sigma_chain$type <- mu_chain$type
    }
    cr_chain <- merge(mu_chain, sigma_chain, by=c("sampno","type","chain"))

    if(options$priming){
        c_chain <- chain[,c(1,which(colnames(chain) == "c"),ncol(chain)),drop=FALSE]
        c_chain <- reshape2::melt(c_chain, c("sampno","chain"))
        colnames(c_chain) <- c("sampno","chain","names","mu")
        c_chain$type <- "priming"
        c_chain <- c_chain[,c("sampno","chain","mu","type")]
        
        beta_chain <- chain[,c(1,which(colnames(chain) == "beta"),ncol(chain)),drop=FALSE]
        beta_chain <- reshape2::melt(beta_chain, c("sampno","chain"))
        colnames(beta_chain) <- c("sampno","chain","names","sigma")
        beta_chain$type <- "priming"
        beta_chain <- beta_chain[,c("sampno","chain","sigma","type")]
        
        priming_chain <- merge(c_chain, beta_chain, by=c("sampno","type","chain"))
        cr_chain <- rbind(cr_chain, priming_chain)
    }

    for(i in 1:((x_max*10)+1)){
        x <- (i-1)/10
        y <- cr_chain[,"mu"] - x*cr_chain[,"sigma"]
        y[y < 0] <- 0
        cr_chain <- cbind(cr_chain,y)
    }
    colnames(cr_chain) <- c("sampno","type","chain","mu","sigma",as.character(seq(0,x_max,by=1/x_max)))
    
    quantiles <- plyr::ddply(cr_chain[,c("type",as.character(seq(0,x_max,by=1/x_max)))],"type",
                             function(x){
                                 y <- apply(x[,2:ncol(x)],2,quantile,probs=c(0.025,0.5,0.975))
                                 y1 <- apply(x[,2:ncol(x)],2,mean)
                                 y <- rbind(y, y1)
                                 colnames(y) = seq(0,x_max,by=1/x_max)
                                 res <- data.frame("quantile"=c("lower","mid","upper","mean"),y)
                                 res
                             })

    lower <- quantiles[quantiles$quantile=="lower",c(1,3:ncol(quantiles))]
    mid <- quantiles[quantiles$quantile=="mid",c(1,3:ncol(quantiles))]
    upper <- quantiles[quantiles$quantile=="upper",c(1,3:ncol(quantiles))]
    mean <- quantiles[quantiles$quantile == "mean",c(1,3:ncol(quantiles))]

    melted_lower <- reshape2::melt(lower,id.vars="type")
    melted_mid <- reshape2::melt(mid,id.vars="type")
    melted_upper <- reshape2::melt(upper,id.vars="type")
    melted_mean <- reshape2::melt(mean,id.vars="type")
    colnames(melted_lower) <- c("type","x","lower")
    colnames(melted_mid) <- c("type","x","mid")
    colnames(melted_upper) <- c("type","x","upper")
    colnames(melted_mean) <- c("type","x","mean")

    all_melted <- join(melted_lower,melted_mid,by=c("type","x"))
    all_melted <- join(all_melted,melted_upper,by=c("type","x"))
    all_melted <- join(all_melted,melted_mean,by=c("type","x"))
    all_melted$x <- as.numeric(all_melted$x)/10
    all_melted$x <- as.numeric(all_melted$x) - 0.1

    all_melted$type <- type_names[all_melted$type]
    all_melted$type = factor(all_melted$type,ordered_names)
    p <- ggplot(all_melted) + 
        geom_line(aes(x=x,y=mid)) + 
        geom_ribbon(aes(x=x,ymin=lower,ymax=upper),alpha=0.2) + 
        facet_wrap(~type)+
        theme_bw() +
        theme(legend.position="none") +
        ylab("Realised boost") +
        xlab("Antigenic distance, x")

    xs <- c("V0"=0,"V1"=6.226721038, "V2"=6.1271143591, "V3"=0.6959541652, "V4"=1.5495453075)
    xs_notrounded <- xs <- xs[order(xs)]
    xs <- round(xs,digits = 1)
    all_melted$x <- round(all_melted$x,1)
    ys <- ddply(all_melted,"type",function(x) x[x$x %in% xs,"mid"])
    results_table <- as.data.frame(t(round(ys[,2:ncol(ys)],3)))
    results_table <- cbind(round(xs_notrounded,3),results_table)
    colnames(results_table) <- c("Antigenic distance",as.character(ys$type))
    
    ys <- reshape2::melt(ys)
    ys$variable <- as.numeric(ys$variable)
    ys$variable <- xs[ys$variable]

    ys$type <- type_names[ys$type]
    ys$type <- factor(ys$type, ordered_names)
    
    p_final <- p + 
        geom_segment(data=ys,aes(x=variable,y=value,xend=variable,yend=0,group=type),linetype="dashed",alpha=0.5)+ 
        geom_segment(data=ys,aes(x=variable,y=value,xend=0,yend=value,group=type),linetype="dashed",alpha=0.5) +
        scale_x_continuous(breaks=seq(0,10,by=1),labels=seq(0,10,by=1),expand=c(0,0)) +
        theme(text=element_text(family="Arial"),
              panel.grid.minor=element_blank(),
              panel.spacing=unit(1,"lines"),
              axis.text=element_text(size=8),
              strip.text=element_text(size=8))+
        scale_y_continuous(breaks=seq(0,15,by=2),labels=seq(0,15,by=2),expand=c(0,0))
    return(list(results_table, p_final))
}


y0_mod_plot <- function(chain, x_max, mu, obs_max=12, xby=10, nsamp=100,ymax=12){
##############
    ## y0 mod plot
##############
    x <- seq(0,x_max,by=0.1)
    
    prime_mu <- function(mu, y0, y0_mod, boost_limit){
        y0[y0 >= boost_limit] <- boost_limit
        tmp <- y0_mod*y0 + mu
        tmp[tmp < 0] <- 0
        return(tmp)
    }
    xbreaks <- c(seq(0,x_max,by=xby),obs_max)
    xbreaks <- xbreaks[xbreaks != floor(obs_max/xby)*xby]
    xbreaks <- c(xbreaks,obs_max)
    tmp <- seq(0,x_max,by=xby)
    tmp <- tmp[tmp != floor(obs_max/xby)*xby]
    xlabels <- c(tmp,paste0("Max observable\ntitre=",obs_max))

    samps <- sample(1:nrow(chain), nsamp)
    
    mus <- matrix(nrow=nsamp,ncol=length(x))
    for(i in seq_along(samps)){
        index <- samps[i]
        mus[i,] <- prime_mu(mu, x, chain[index,"y0_mod"],chain[index,"boost_limit"])
    }
    dat <- apply(mus, 2, function(x) c(mean(x),quantile(x,c(0.025,0.975))))
    dat <- as.data.frame(t(dat))
    colnames(dat) <- c("y","lower","upper")
    dat$x <- x
    dat2 <- dat[dat$x==obs_max,]

    ylabels <- unique(c(seq(0,ymax,by=2),signif(dat2[1,"lower"],3),signif(dat2[1,"y"],3),
                signif(dat2[1,"upper"],3),mu))
 
    ylabels <- ylabels[!duplicated(round(ylabels))]
    print(ylabels)
    ylabels <- c(ylabels, mu)
    ybreaks <- ylabels
    ylabels[which(ylabels==mu)] <- paste0("μ=",mu)
    
    y0_p2 <- ggplot(dat) + 
        geom_line(aes(x=x,y=y)) + 
        geom_ribbon(aes(x=x,ymax=upper,ymin=lower),alpha=0.3) +
        theme_bw() +
        theme(axis.text.x=element_text(size=10,colour="black",family="Arial"),
              axis.text.y=element_text(size=10,colour="black",family="Arial"))+
        scale_x_continuous(expand=c(0,0), breaks=xbreaks,labels=xlabels)+
        geom_segment(data=dat2,aes(x=x,y=upper,xend=x,yend=0),linetype="dashed",alpha=0.5)+ 
        geom_segment(data=dat2,aes(x=0,y=lower,xend=x,yend=lower),linetype="dashed",alpha=0.5)+
        geom_segment(data=dat2,aes(x=0,y=upper,xend=x,yend=upper),linetype="dashed",alpha=0.5)+
                                        #geom_segment(data=dat2,aes(x=x,y=y,xend=x,yend=0),linetype="dashed",alpha=0.5)+ 
        geom_segment(data=dat2,aes(x=0,y=y,xend=x,yend=y),linetype="dashed",alpha=0.5)+ 
        scale_y_continuous(expand=c(0,0),limits=c(0,ymax),breaks=ybreaks,labels=ylabels)+
        ylab("Realised boost") +
        xlab("Initial titre at time of exposure, y0") 
        #ggtitle(paste0("Estimated titre dependent boosting relationship, μ = ",mu))
    return(y0_p2)
}

#########################
## Combined density plot
#########################
combined_density_plot <- function(estimates, parName, title_parName, ymin,ymax, priormin,priormax,
                                  saveDir="~/",savePNG=FALSE,saveEPS=FALSE,
                                  blueLine=-1000,redLine=-1000, fillBy=NULL){
  p <- ggplot(estimates[estimates$Parameter.name == parName,]) + 
    geom_pointrange(aes_string(x="runID",y="Median",ymax="X97.5..CI",ymin="X2.5..CI",col=fillBy),size=0.25,fatten=0.5) + 
    facet_wrap(~Exposure.Type) +
    theme_bw() +
      theme(text=element_text(size=8,family="Arial"),
            strip.text=element_text(size=8,family="Arial"),
          axis.text.x=element_text(angle=45,hjust=1,size=8),
          axis.text.y=element_text(size=8),
          legend.position="bottom") +
    ylab(paste0("Estimates for ",title_parName)) +
    coord_cartesian(ylim=c(ymin,ymax)) + 
    geom_hline(yintercept=c(priormin,priormax),linetype="dashed",alpha=0.5) + 
    geom_hline(yintercept=blueLine,linetype="dashed",col="blue",alpha=0.5) + 
    geom_hline(yintercept=redLine,linetype="dashed",col="red",alpha=0.5)

  if(savePNG){
    png(paste0(saveDir,title_parName,"_densities.png"),width=6,height=5,res=300,unit="in",family="Arial")
    print(p)
    dev.off()
  }
  if(saveEPS){
    cairo_ps(paste0(saveDir,title_parName,"_densities.eps"),width=6,height=5,family="Arial")
    print(p)
    dev.off()
  }
  p
}
