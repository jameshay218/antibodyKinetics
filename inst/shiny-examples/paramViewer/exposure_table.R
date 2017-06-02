##################################
## Exposure Handsontable
##################################
## This bit of code makes the exposureTab an interactive table that we can edit.
observe({
    if(!is.null(inputs$exposure_table)){
        parameters[["previous"]] <- isolate(parameters[["exposureTab"]])
        exposureTab = hot_to_r(inputs$exposure_table)
    } else {
        if(is.null(parameters[["exposureTab"]]))
            exposureTab <- NULL
        else
            exposureTab <- parameters[["exposureTab"]]
    }
    parameters[["exposureTab"]] <- exposureTab
})

output$exposure_table <- renderRHandsontable({
    exposureTab <- parameters$exposureTab
    if(!is.null(exposureTab)){
        rhandsontable(exposureTab, useTypes = FALSE, stretchH="all")
    }})

output$export_exposures <- downloadHandler(
    filename = "exposureTab.csv",
    content=function(file){
        write.csv(parameters$exposureTab,file,row.names=FALSE)
    })

observeEvent(inputs$upload_exposures,{
    newTab <- read.csv(isolate(inputs$exposure_tab_input$datapath),stringsAsFactors=FALSE)
    if(!is.null(newTab)){
        ## If the new parameter table has types in it, need to set the flags etc correctly
        if(any(strong_types%in% newTab$type)){
            updateSelectInput(session,"typing_flags",selected=2)
        } else if(any(weak_types %in% newTab$type)){
            updateSelectInput(session,"typing_flags",selected=1)
        } else {
            updateSelectInput(session,"typing_flags",selected=0)
        }
        updateNumericInput(session,"tmax",value=max(newTab$end))
        updateNumericInput(session,"n_strains",value=length(unique((newTab$strain))))
        parameters$exposureTab <- newTab

        newParTab <- NULL
        
        ## THEN NEED TO UPDATE THE OVERALL PARAMETER TABLE
        if(inputs$typing_flags == 0 & inputs$cr_flags == 0){
            ## If no typing and no CR add entry for each ID and strain
            for(id in unique(newTab$id)){
                tmp <- newTab[newTab$id == id,]
                exposure <- tmp$exposure[1]
                for(strain in unique(tmp$strain)){
                    mu <- runif(1,min_mu,max_mu)
                    dp <- runif(1,min_dp,max_dp)
                    tp <- runif(1,min_tp,max_tp)
                    ts <- runif(1,min_ts,max_ts)
                    m <- runif(1,min_m,max_m)
                    
                    newParTab <- rbind(newParTab, data.frame(names=c("mu","tp","dp","ts","m"),
                                                             id=id,
                                                             values=c(mu,tp,dp,ts,m),
                                                             type="all",
                                                             exposure=exposure,
                                                             strain=strain,
                                                             order=NA,
                                                             fixed=1,steps=0.1,lower_bound=0,
                                                             upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                                             stringsAsFactors=FALSE
                                                             ))
                    
                }
            }
        } else if(inputs$typing_flags == 0 & inputs$cr_flags != 0){
            ## If no typing and CR, add entry for each ID
            for(id in unique(newTab$id)){
                     mu <- runif(1,min_mu,max_mu)
                    dp <- runif(1,min_dp,max_dp)
                    tp <- runif(1,min_tp,max_tp)
                    ts <- runif(1,min_ts,max_ts)
                    m <- runif(1,min_m,max_m)
                    
                tmp <- newTab[newTab$id == id,]
                exposure <- tmp$exposure[1]
                newParTab <- rbind(newParTab, data.frame(names=c("mu","tp","dp","ts","m"),
                                                         id=id,
                                                         values=c(mu,tp,dp,ts,m),
                                                         type="all",
                                                         exposure=exposure,
                                                         strain="all",
                                                         order=NA,
                                                         fixed=1,steps=0.1,lower_bound=0,
                                                         upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                                         stringsAsFactors=FALSE
                                                         ))
            }
        } else if(inputs$typing_flags != 0 & inputs$cr_flags == 0){
            ## If typing but no CR, add entry for each type, exposure and strain combo
            tmp <- unique(newTab[,c("type","exposure","strain")])
            for(i in 1:nrow(tmp)){
                     mu <- runif(1,min_mu,max_mu)
                    dp <- runif(1,min_dp,max_dp)
                    tp <- runif(1,min_tp,max_tp)
                    ts <- runif(1,min_ts,max_ts)
                    m <- runif(1,min_m,max_m)
                    
                newParTab <- rbind(newParTab, data.frame(names=c("mu","tp","dp","ts","m"),id="all",
                                                         values=c(mu,tp,dp,ts,m),type=tmp[i,"type"],
                                                         exposure=tmp[i,"exposure"],strain=tmp[i,"strain"],
                                                         order=NA,fixed=1,steps=0.1,lower_bound=0,
                                                         upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                                         stringsAsFactors=FALSE
                                                         ))
                     
            }
        } else {
            ## Otherwise, just add entry for each type
            for(type in unique(newTab$type)){
                mu <- runif(1,min_mu,max_mu)
                dp <- runif(1,min_dp,max_dp)
                tp <- runif(1,min_tp,max_tp)
                ts <- runif(1,min_ts,max_ts)
                m <- runif(1,min_m,max_m)
                
                newParTab <- rbind(newParTab, data.frame(names=c("mu","tp","dp","ts","m"),id="all",
                                                         values=c(mu,tp,dp,ts,m),type=type,
                                                         exposure="all",strain="all",
                                                         order=NA,fixed=1,steps=0.1,lower_bound=0,
                                                         upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                                         stringsAsFactors=FALSE
                                                         ))
            }
        }
    }
    parameters$parTab <- newParTab
})
