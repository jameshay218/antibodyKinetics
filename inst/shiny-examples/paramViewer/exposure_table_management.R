###########################################################
## EXPOSURE TABLE MANAGEMENT
## This file is needed to allow users to add, update, remove and delete the exposure table
## This all relates to the "Exposures" tab
###########################################################
## When "add" is clicked, create the appropriate entries and put this into the exposure table
observeEvent(inputs$add_exposure,{
    ## Was this a primed exposure?
    if(inputs$is_primed)
        primed <- 1
    else
        primed <- 0
        
    ## The new entry needs the right ID and parameter values. If the given ID exists,
    ## then need to just update rather than add.
    ## The ID is created from the inputs
    order <- 1
    new_id <- isolate(paste0("G",inputs$exposure_group,"E",order,"S",inputs$exposure_strain))
    print(new_id)

    mu <- runif(1,0,max_mu)
    tp <- runif(1,0,max_tp)
    dp <- runif(1,0,1)
    ts <- runif(1,0,max_ts)
    m <- runif(1,min_m,0)
    
    exposureRow <- isolate(data.frame(id=new_id,values=inputs$ti,type=inputs$type,exposure=inputs$exposure_strain,
                                      strain=inputs$exposure_affects,order=order,primed=primed,
                                      group=inputs$exposure_group,end=inputs$tmax,next_t=inputs$tmax))

    #parTabRow <- isolate(data.frame(names=c("mu","tp","dp","ts","m"),
    #                                id=inputs$ID,
    #                                values=c(inputs$mu,inputs$tp,inputs$dp,inputs$ts,inputs$m),
    #                                type=inputs$type,
    #                                exposure=inputs$exposure_strain,
    #                                strain=rep(inputs$exposure_affects,each=5),
    #                                order=1,
    #                                fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1)))
    parameters$exposureTab <- isolate(rbind(parameters$exposureTab, exposureRow))
    #parameters$parTab <- isolate(rbind(parameters$parTab, parTabRow))
})

observeEvent(inputs$remove_exposure,{
    newDat <- isolate(parameters$exposureTab[parameters$exposureTab$id != inputs$ID,])
    newParTab <- isolate(parameters$parTab[parameters$parTab$id != inputs$exposure_select,])
    parameters$exposureTab <- newDat
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

