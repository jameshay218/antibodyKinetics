###########################################################
## EXPOSURE TABLE MANAGEMENT
## This file is needed to allow users to add, update, remove and delete the exposure table
## This all relates to the "Exposures" tab
###########################################################
## When "add" is clicked, create the appropriate entries and put this into the exposure table
observeEvent(inputs$add_exposure,{
    ## Was this a primed exposure?
    if(inputs$is_primed){ primed <- 1
    }else{ primed <- 0}
    
    ## The new entry needs the right ID and parameter values. If the given ID exists,
    ## then need to just update rather than add.
    ## The ID is created from the inputs
    mu <- runif(1,min_mu,max_mu)
    tp <- runif(1,min_tp,max_tp)
    dp <- runif(1,min_dp,max_dp)
    ts <- runif(1,min_ts,max_ts)
    m <- runif(1,min_m,max_m)

    if(any(c(length(inputs$exposure_affects) == 0),
           is.na(inputs$exposure_ti),
           is.na(inputs$exposure_type),
           is.na(inputs$exposure_strain),
           is.na(inputs$exposure_group))) {
        print("Error - missing inputs")
    }
    else {
        ## Get the time of the next infection, if any, and get the
        ## number of exposures before this in the group
        new_pars <- add_order_nextt(inputs,parameters)
        next_t <- new_pars[["next_t"]]
        order <- new_pars[["order"]]
        newTab <- new_pars[["newTab"]]
        new_id <- new_pars[["new_id"]]
        
        exposureRow <- isolate(data.frame(id=new_id,values=inputs$exposure_ti,type=inputs$exposure_type,
                                          exposure=inputs$exposure_strain,
                                          strain=inputs$exposure_affects,order=order,primed=primed,
                                          group=inputs$exposure_group,end=inputs$tmax,next_t=next_t,
                                          stringsAsFactors=FALSE))
        parameters$exposureTab <- newTab

        ## If no typing or cross reactivity, need to add an entry for each unique strain/exposure/time
        ## combination
        if(inputs$typing_flags == 0 & inputs$cr_flags == 0){
            parTabRow <- do.call("rbind",lapply(inputs$exposure_affects, function(affects){
                data.frame(names=c("mu","tp","dp","ts","m"),
                           id=new_id,
                           values=c(mu,tp,dp,ts,m),
                           type=inputs$exposure_type,
                           exposure=inputs$exposure_strain,
                           strain=affects,
                           order=NA,
                           fixed=1,steps=0.1,lower_bound=0,upper_bound=c(max_mu,max_tp,1,max_ts,1),
                           stringsAsFactors=FALSE
                           )
            }))
            ## If typing but no cross reactivity, add entry for each unique
            ## strain/exposure/type combo
        } else if(inputs$typing_flags != 0 & inputs$cr_flags == 0){
            parTabRow <- do.call("rbind",lapply(inputs$exposure_affects, function(affects){
                ## For each strain this affects
                newRow <- NULL
                if(is.null(parameters$parTab) || nrow(parameters$parTab[parameters$parTab$exposure == inputs$exposure_strain &
                                                                        parameters$parTab$strain == affects &
                                                                        parameters$parTab$type == inputs$exposure_type,]) == 0){
                    newRow <- data.frame(names=c("mu","tp","dp","ts","m"),
                                         id=new_id,
                                         values=c(mu,tp,dp,ts,m),
                                         type=inputs$exposure_type,
                                         exposure=inputs$exposure_strain,
                                         strain=affects,
                                         order=NA,
                                         fixed=1,steps=0.1,lower_bound=0,upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                         stringsAsFactors=FALSE
                                         )
                }
                newRow
            }))
            ## If no typing but cross reactivity, need entry for each unique
            ## exposure/time combination
        } else if(inputs$typing_flags == 0 & inputs$cr_flags != 0){
            parTabRow <- NULL
            if(is.null(parameters$parTab) || nrow(parameters$parTab[parameters$parTab$exposure == inputs$exposure_strain &
                                                                    parameters$parTab$values == inputs$exposure_ti,]) == 0){
                parTabRow <- data.frame(names=c("mu","tp","dp","ts","m"),
                                        id=new_id,
                                        values=c(mu,tp,dp,ts,m),
                                        type=inputs$exposure_type,
                                        exposure=inputs$exposure_strain,
                                        strain="all",
                                        order=NA,
                                        fixed=1,steps=0.1,lower_bound=0,upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                        stringsAsFactors=FALSE
                                        )
            }
            ## Otherwise, we have tying and cross reactivity, in which case we need entry for each exposure type
        } else {
            parTabRow <- NULL
            if(is.null(parameters$parTab) || nrow(parameters$parTab[parameters$parTab$type == inputs$exposure_type,]) == 0){
                parTabRow <- data.frame(names=c("mu","tp","dp","ts","m"),
                                        id=new_id,
                                        values=c(mu,tp,dp,ts,m),
                                        type=inputs$exposure_type,
                                        exposure="all",
                                        strain="all",
                                        order=NA,
                                        fixed=1,steps=0.1,lower_bound=0,upper_bound=c(max_mu,max_tp,1,max_ts,1),
                                        stringsAsFactors=FALSE
                                        )
            }
            
        }
        parameters$exposureTab <- isolate(rbind(parameters$exposureTab, exposureRow))
        parameters$exposureTab <- parameters$exposureTab[order(parameters$exposureTab$group, parameters$exposureTab$values),]
        parameters$parTab <- isolate(rbind(parameters$parTab, parTabRow))
    }
})

observeEvent(inputs$remove_exposure,{
    removed_exposure <- parameters$exposureTab[parameters$exposureTab$id == inputs$exposure_id,]
    parameters$exposureTab <- remove_order_nextt(inputs,parameters)
    newParTab <- parameters$parTab
    if(is.null(parameters$exposureTab) | nrow(parameters$exposureTab) == 0){
        newParTab <- NULL
    } else if(inputs$typing_flags == 0 & inputs$cr_flags == 0){
        ## If no typing and no cross-reactivity, just remove this ID
        newParTab <- isolate(parameters$parTab[parameters$parTab$id != inputs$exposure_id,])
    } else if(inputs$typing_flags != 0 & inputs$cr_flags == 0){
        ## If typing but no cross reactivity, remove entry if this is the last
        ## strain/exposure/type combo
        newParTab <- dplyr::semi_join(newParTab,parameters$exposureTab,by=c("exposure","type","strain"))
    } else if(inputs$typing_flags == 0 & inputs$cr_flags != 0){
        ## If no typing but cross reactivity, remove entry entry is this is the last
        ## exposure of this strain
        newParTab <- isolate(parameters$parTab[parameters$parTab$id != inputs$exposure_id,])
    } else {
        ## Otherwise, we have typing and cross reactivity, in which case we need entry remove
        ## only if this is the last of this exposure type
        if(nrow(parameters$exposureTab[parameters$exposureTab$type == removed_exposure[1,"type"],])==0){
            newParTab <- parameters$parTab[!(parameters$parTab$type == removed_exposure[1,"type"]),]
        }
    }
    if(is.null(newParTab))
        parameters$parTab <- NULL
    else
        parameters$parTab <- newParTab
})

observeEvent(inputs$clear_exposures,{
    parameters$exposureTab <- NULL
    parameters$parTab <- NULL
})
