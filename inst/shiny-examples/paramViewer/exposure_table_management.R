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
    mu <- runif(1,0,max_mu)
    tp <- runif(1,0,max_tp)
    dp <- runif(1,0,1)
    ts <- runif(1,0,max_ts)
    m <- runif(1,min_m,0)

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
        print(new_pars)
        next_t <- new_pars[["next_t"]]
        order <- new_pars[["order"]]
        newTab <- new_pars[["newTab"]]
        new_id <- new_pars[["new_id"]]

        parameters$exposureTab <- newTab
        exposureRow <- isolate(data.frame(id=new_id,values=inputs$exposure_ti,type=inputs$exposure_type,
                                          exposure=inputs$exposure_strain,
                                          strain=inputs$exposure_affects,order=order,primed=primed,
                                          group=inputs$exposure_group,end=inputs$tmax,next_t=next_t,
                                          stringsAsFactors=FALSE))
        ## If no typing or cross reactivity, need to add an entry for each unique strain/exposure/time
        ## combination       
        if(inputs$typing_flag == 0 & inputs$cr_flags == 0){
            parTabRow <- do.call("rbind",lapply(inputs$exposure_affects, function(affects){
                data.frame(names=c("mu","tp","dp","ts","m"),
                           id=new_id,
                           values=c(mu,tp,dp,ts,m),
                           type=inputs$exposure_type,
                           exposure=inputs$exposure_strain,
                           strain=affects,
                           order=1,
                           fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1),
                           stringsAsFactors=FALSE
                           )
            }))
            parameters$parTab <- isolate(rbind(parameters$parTab, parTabRow))
            ## If typing but no cross reactivity, add entry for each unique
            ## strain/exposure/type combo
        } else if(inputs$typing_flag != 0 & inputs$cr_flag == 0){
            parTabRow <- do.call("rbind",lapply(inputs$exposure_affects,
                                                ## For each strain this affects
                                                function(affects){
                                                    newRow <- NULL
                                                    if(nrow(parameters$parTab[parameters$parTab$exposure == inputs$exposure_strain &
                                                                              parameters$parTab$strain == affects &
                                                                              parameters$parTab$type == inputs$exposure_type,]) == 0){
                                                        newRow <- data.frame(names=c("mu","tp","dp","ts","m"),
                                                                             id=new_id,
                                                                             values=c(mu,tp,dp,ts,m),
                                                                             type=inputs$exposure_type,
                                                                             exposure=inputs$exposure_strain,
                                                                             strain=affects,
                                                                             order=1,
                                                                             fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1),
                                                                             stringsAsFactors=FALSE
                                                                             )
                                                    }
                                                    newRow
                                                }))
            ## If no typing but no cross reactivity, need entry for each unique
            ## exposure/time combination
        } else if(inputs$typing_flag == 0 & inputs$cr_flag != 0){
            parTabRow <- NULL
            if(nrow(parameters$parTab[parameters$parTab$exposure == inputs$exposure_strain &
                                      parameters$parTab$values == inputs$exposure_ti,]) == 0){
                parTabRow <- data.frame(names=c("mu","tp","dp","ts","m"),
                                        id=new_id,
                                        values=c(mu,tp,dp,ts,m),
                                        type=inputs$exposure_type,
                                        exposure=inputs$exposure_strain,
                                        strain="all",
                                        order=1,
                                        fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1),
                                        stringsAsFactors=FALSE
                                        )
            }
            ## Otherwise, we have tying and cross reactivity, in which case we need entry for each exposure type
        } else {
            parTabRow <- NULL
            if(nrow(parameters$parTab[parameters$parTab$type == inputs$exposure_type,]) == 0){
                parTabRow <- data.frame(names=c("mu","tp","dp","ts","m"),
                                        id=new_id,
                                        values=c(mu,tp,dp,ts,m),
                                        type=inputs$exposure_type,
                                        exposure=inputs$exposure_strain,
                                        strain="all",
                                        order=1,
                                        fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1),
                                        stringsAsFactors=FALSE
                                        )
            }
            
        }
        
        
        print(parTabRow)
        parameters$exposureTab <- isolate(rbind(parameters$exposureTab, exposureRow))
        parameters$exposureTab <- parameters$exposureTab[order(parameters$exposureTab$group, parameters$exposureTab$values),]
        
    }
})

observeEvent(inputs$remove_exposure,{
    print("Removing exposure")
    print(inputs$exposure_id)
    parameters$exposureTab <- remove_order_nextt(inputs,parameters)
    newParTab <- isolate(parameters$parTab[parameters$parTab$id != inputs$exposure_id,])
    parameters$parTab <- newParTab
})

observeEvent(inputs$update_exposure,{
    primed <- 0
    if(inputs$primed) primed <- 1
    exposureRow <- isolate(data.frame(id=inputs$ID,
                                      values=inputs$ti,
                                      type=inputs$type,
                                      exposure=inputs$exposure_strain,
                                      strain=inputs$exposure_affects,
                                      order=1,
                                      primed=primed,
                                      group=inputs$group,
                                      end=inputs$tmax,
                                      next_t=inputs$tmax))

    parTabRow <- isolate(data.frame(names=c("mu","tp","dp","ts","m"),
                                    id=inputs$ID,
                                    values=c(inputs$mu,inputs$tp,inputs$dp,inputs$ts,inputs$m),
                                    type=inputs$type,
                                    exposure=inputs$exposure_strain,
                                    strain=rep(inputs$exposure_affects,each=5),
                                    order=1,
                                    fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1)))
    parameters$exposureTab <- isolate(rbind(parameters$exposureTab, exposureRow))
    parameters$parTab <- isolate(rbind(parameters$parTab, parTabRow))
})

observeEvent(inputs$clear_exposures,{
    parameters$exposureTab <- NULL
    parameters$parTab <- parameters$parTab[!parameters$parTab %in% get_ids(parameters),]
})
