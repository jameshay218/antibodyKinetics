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
    ## Check if we have the data to make the plot
    
})
