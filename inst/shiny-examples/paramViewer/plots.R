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
                                 exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                 lower_bound=c(-1000,0,0,0),upper_bound=c(0,1,1,100),stringsAsFactors=FALSE)
        tmpCrTab <- parameters$crTab[parameters$crTab$names %in% get_available_exposure_types_cr(),]
        cr_values <- tmpCrTab$values
        cr_names <- tmpCrTab$names
        bot_parTab <- data.frame(names=c("beta","c",cr_names,"y0_mod"),id="all",
                             values=c(inputs$beta,inputs$c,cr_values,inputs$y0_mod),
                             exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                             lower_bound=c(-20,0,-20,-20),upper_bound=c(2,20,2,2),stringsAsFactors=FALSE)

        mod_parTab <- data.frame(names="mod",id=NA,values=c(inputs$mod1,inputs$mod2,inputs$mod3,inputs$mod4),
                                 exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                 lower_bound=0,upper_bound=1,stringsAsFactors=FALSE)
        
    }
})
