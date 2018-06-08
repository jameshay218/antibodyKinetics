###########################################################
## PLOTS
###########################################################
output$protocol_plot <- renderPlot({
    if(!is.null(parameters[["exposureTab"]]) && nrow(parameters[["exposureTab"]]) > 0){
        tmpTab <- parameters[["exposureTab"]]
        ggplot(tmpTab) +
            geom_vline(aes(xintercept=values,col=exposure,group=exposure)) +
            geom_rug(sides="b") +
            geom_text(aes(x=values+2,label=type,col=exposure,y=2.5),angle=90)+
            scale_x_continuous(limits=c(0,inputs$tmax)) +
            facet_wrap(~group,ncol=2) +
            ylab("")+
            xlab("Time (days)")+
            theme_bw() +
            theme(axis.text.y=element_blank())
    }
})
output$main_plot <- renderPlot({
    if(!is.null(parameters$parTab) & !is.null(parameters$exposureTab)){
        ## Check if we have the data to make the plot
        top_parTab <- data.frame(names=c("lower_bound","S","EA","MAX_TITRE"), id="all",
                                 values=c(inputs$lower_bound,0.79,0.2,inputs$max_titre),
                                 type="all",
                                 exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                 lower_bound=c(-1000,0,0,0),upper_bound=c(0,1,1,100),stringsAsFactors=FALSE)
        print(get_available_exposure_types_cr())

        ## If flag is checked for titre dependent boosting, add these entries. Otherwise, use
        ## values that remove this mechanism
        if(inputs$titre_dependent_boosting){
            boost_limit <- inputs$boost_limit
            y0_mod <- inputs$y0_mod
            boost_limit_lower <- 0
            boost_limit_upper <- 12

            y0_mod_lower <- -1
            y0_mod_upper <- 1
        } else {
            boost_limit <- -1
            y0_mod <- -1000

            boost_limit_lower <- -1
            boost_limit_upper <- 12
            
            y0_mod_lower <- -1000
            y0_mod_upper <- 1
            
        }
        
        if(inputs$cr_flags != 0){
            tmpCrTab <- parameters$crTab[parameters$crTab$names %in% get_available_exposure_types_cr(),]
            cr_values <- tmpCrTab$values
            cr_names <- tmpCrTab$names
            print(tmpCrTab)
            bot_parTab <- data.frame(names=c("beta","c",rep("sigma",length(cr_names)),"y0_mod","boost_limit"),id="all",
                                     values=c(inputs$beta,inputs$c,cr_values,y0_mod,boost_limit),
                                     type=c("all","all",cr_names,"all","all"),
                                     exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                     lower_bound=c(0,0,rep(0,length(cr_names)),y0_mod_lower,boost_limit_lower),upper_bound=c(100,20,rep(100,length(cr_names)),y0_mod_upper,boost_limit_upper),stringsAsFactors=FALSE)
        } else {
            bot_parTab <- data.frame(names=c("beta","c","sigma","y0_mod","boost_limit"),id="all",
                                     values=c(inputs$beta,inputs$c,0,y0_mod,boost_limit),
                                     type=c("all","all","all","all","all"),
                                     exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                     lower_bound=c(0,0,0,y0_mod_lower,boost_limit_lower),upper_bound=c(100,20,100,y0_mod_upper,boost_limit_upper),stringsAsFactors=FALSE)
        }

        mod_parTab <- data.frame(names="mod",id=NA,values=c(inputs$mod1,inputs$mod2,inputs$mod3,inputs$mod4),
                                 type="all",exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                 lower_bound=0,upper_bound=1,stringsAsFactors=FALSE)

        distance_parTab <- data.frame(names="x",id=NA,values=parameters$antigenicDistTab$Distance,
                                      type="all",exposure=parameters$antigenicDistTab$Strain.1,
                                      strain=parameters$antigenicDistTab$Strain.2,
                                      order=NA,fixed=1,steps=0.1,lower_bound=0,upper_bound=13,
                                      stringsAsFactors=FALSE)
        tmpTab <- parameters$parTab
        tmpTab[tmpTab$names == "m","values"] <- exp(tmpTab[tmpTab$names == "m","values"])
        
        parTab <- rbind(top_parTab,tmpTab,bot_parTab,distance_parTab,mod_parTab)
        #print(parTab)
        typing <- inputs$typing_flags != 0
        cross_reactivity <- inputs$cr_flags != 0
        print(inputs$form)

        f <- create_model_group_func_cpp(parTab,parameters$exposureTab,version="model",form=as.character(inputs$form),typing=typing,cross_reactivity=cross_reactivity)

        times <- seq(0,100,by=0.1)

        y <- f(parTab$values, times)
        y <- as.data.frame(y)

        n_strains <- length(unique(parameters$exposureTab$strain))
        n_groups <- length(unique(parameters$exposureTab$group))

        y$group <- rep(1:n_groups,each=n_strains)
        y$strain <- rep(1:n_strains,n_groups)

        colnames(y) <- c(times,"group","strain")
        dat <- reshape2::melt(y,id.vars=c("group","strain"))
        colnames(dat) <- c("group","strain","times","value")
        dat$times <- as.numeric(as.character(dat$times))
        dat$strain <- as.factor(dat$strain)
        dat$group <- as.factor(dat$group)
        ggplot(dat) + geom_line(aes(x=times,col=strain,y=value)) + facet_wrap(~group) + theme_bw()
    }
})
