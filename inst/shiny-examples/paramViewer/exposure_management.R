###########################################################
## EXPOSURE INPUT MANAGEMENT
## This file is needed to allow users to modify values relating to exposures
## This all relates to the "Exposures" tab
###########################################################
#############
## INPUTS
#############
## Choose exposure ID input. This should be a list of all of the IDs in the
## current exposure table. If nothing else in there, then this should be set to
## "New". Possibly use this in other tab too
output$choose_exposure_id <- renderUI({
    choices <- c("New",get_ids(parameters))
    selectInput("exposure_id","ID",
                choices=choices)
})

## Choose exposure type
output$choose_exposure_type <- renderUI({
    types <- get_types(inputs)
    selectInput("exposure_type","Type",choices=types)
})

## Choose the exposure infection time
output$choose_exposure_ti <- renderUI({numericInput("exposure_ti","Time",0,min=0,max=inputs$tmax)})

## Choose the exposure group from 1 to 10
output$choose_exposure_group <- renderUI({numericInput("exposure_group","Group",1,min=1,max=10)})

## When choosing which exposure strain to use, only display the number of strains we have
output$choose_exposure_strain <- renderUI({
    selectInput("exposure_strain","Exposure",
                choices=as.list(exposure_strains[1:inputs$n_strains]),
                selected=1)
})

## Similarly, only display the present strains when determining which titres an exposure
## affects
output$choose_exposure_affects <- renderUI({
    checkboxGroupInput("exposure_affects","Affects",
                       choices=as.list(exposure_strains[1:inputs$n_strains]),
                       selected=1,
                       inline=TRUE)
})

## Check box - if ticked, this exposure is primed
output$choose_is_primed <- renderUI({checkboxInput("is_primed","Primed",FALSE)})

## Sliders for exposure specific inputs
output$select_mu <- renderUI({sliderInput("mu", "Boost", value=8, min=0, max=max_mu,step=0.1)})
output$select_tp <- renderUI({sliderInput("tp", "Time to peak", value=12, min=0, max=max_tp,step=0.1)})
output$select_dp <- renderUI({sliderInput("dp", "Proportional drop", value=0.5, min=0, max=1,step=0.001)})
output$select_ts <- renderUI({sliderInput("ts", "Time to 2nd waning phase", value=10, min=0, max=max_ts,step=0.1)})
output$select_m <- renderUI({sliderInput("m", "Long term waning (log)", value=log(0.003), min=min_m, max=0,step=0.01)})


## Management of parameter table when sliders change for "mu"
observeEvent(inputs$mu,{
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "mu","values"] <- inputs$mu)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "mu","values"] <- inputs$mu)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "mu","values"] <- inputs$mu)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "mu","values"] <- inputs$mu)
            }
        }
    }
})

## Management of parameter table when sliders change for "tp"
observeEvent(inputs$tp,{
    if(!is.null(parameters$parTab)){
    if(inputs$cr_flags == 0){
        if(inputs$typing_flags == 0){
            isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "tp","values"] <- inputs$tp)
        } else {
            isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "tp","values"] <- inputs$tp)
        }
    } else {
        if(inputs$typing_flags == 0){
            isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "tp","values"] <- inputs$tp)
        } else {
            isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$names == "tp","values"] <- inputs$tp)
        }
    }
    }
})



## Management of parameter table when sliders change for "dp"
observeEvent(inputs$dp,{
    if(!is.null(parameters$parTab)){
    if(inputs$cr_flags == 0){
        if(inputs$typing_flags == 0){
            isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "dp","values"] <- inputs$dp)
        } else {
            isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "dp","values"] <- inputs$dp)
        }
    } else {
        if(inputs$typing_flags == 0){
            isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "dp","values"] <- inputs$dp)
        } else {
            isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$names == "dp","values"] <- inputs$dp)
        }
    }
    }
})

## Management of parameter table when sliders change for "m"
observeEvent(inputs$m,{
    if(!is.null(parameters$parTab)){
    if(inputs$cr_flags == 0){
        if(inputs$typing_flags == 0){
            isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "m","values"] <- inputs$m)
        } else {
            isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "m","values"] <- inputs$m)
        }
    } else {
        if(inputs$typing_flags == 0){
            isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$names == "m","values"] <- inputs$m)
        } else {
            isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                      parameters$parTab$exposure == inputs$exposure_exposure_select &
                                      parameters$parTab$strain == inputs$exposure_strain_select &
                                      parameters$parTab$names == "m","values"] <- inputs$m)
        }
    }
    }
})


## Management of parameter table when sliders change for "ts"
observeEvent(inputs$ts,{
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "ts","values"] <- inputs$ts)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "ts","values"] <- inputs$ts)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "ts","values"] <- inputs$ts)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "ts","values"] <- inputs$ts)
            }
        }
    }
})


##########
## EVENTS
##########
## Related to the main panel exposure selection
observeEvent(inputs$exposure_select,{
    ## Get the parameter table subset related to this ID/type
    id_choices <- get_ids(parameters)
    #type_choices <- get_types(inputs)
    type_choices <- unique(parameters$parTab$type)

    if(inputs$exposure_select != "none"){
        tmpTab <- parameters$parTab
        if(inputs$typing_flags == 0){
            tmpTab <- tmpTab[tmpTab$id == inputs$exposure_select,]
            strains <- unique(tmpTab$strain)
            exposures <- unique(tmpTab$exposure)
            updateSelectInput(session,"exposure_strain_select",choices=strains)
            updateSelectInput(session,"exposure_exposure_select",choices=exposures)
            if(inputs$cr_flags == 0){
                tmpTab <- tmpTab[tmpTab$strain == strains[1] & tmpTab$exposure == exposures[1],]
            }
        } else {
            tmpTab <- tmpTab[tmpTab$type == inputs$exposure_select,]
            exposures <- unique(tmpTab$exposure)
            strains <- unique(tmpTab$strain)
            updateSelectInput(session,"exposure_exposure_select",choices=exposures)
            updateSelectInput(session,"exposure_strain_select",choices=strains)
            if(inputs$cr_flags == 0){
                tmpTab <- tmpTab[tmpTab$strain == strains[1]  & tmpTab$exposure == exposures[1],]
            }
        }
        mu <- tmpTab[tmpTab$names == "mu","values"]
        tp <- tmpTab[tmpTab$names == "tp","values"]
        dp <- tmpTab[tmpTab$names == "dp","values"]
        ts <- tmpTab[tmpTab$names == "ts","values"]
        m <- tmpTab[tmpTab$names == "m","values"]
        isolate(updateSliderInput(session,inputId="mu",value=mu))
        isolate(updateSliderInput(session,inputId="tp",value=tp))
        isolate(updateSliderInput(session,inputId="dp",value=dp))
        isolate(updateSliderInput(session,inputId="ts",value=ts))
        isolate(updateSliderInput(session,inputId="m",value=m))
    }
})

observeEvent(inputs$exposure_strain_select,{
    tmpTab <- parameters$parTab
    if(inputs$typing_flags == 0){
        tmpTab <- tmpTab[tmpTab$id == inputs$exposure_select & tmpTab$exposure == inputs$exposure_exposure_select,]
    } else {
        tmpTab <- tmpTab[tmpTab$type == inputs$exposure_select & tmpTab$exposure == inputs$exposure_exposure_select,]
    }
    tmpTab <- tmpTab[tmpTab$strain == inputs$exposure_strain_select,]
    mu <- tmpTab[tmpTab$names == "mu","values"]
    tp <- tmpTab[tmpTab$names == "tp","values"]
    dp <- tmpTab[tmpTab$names == "dp","values"]
    ts <- tmpTab[tmpTab$names == "ts","values"]
    m <- tmpTab[tmpTab$names == "m","values"]
    isolate(updateSliderInput(session,inputId="mu",value=mu))
    isolate(updateSliderInput(session,inputId="tp",value=tp))
    isolate(updateSliderInput(session,inputId="dp",value=dp))
    isolate(updateSliderInput(session,inputId="ts",value=ts))
    isolate(updateSliderInput(session,inputId="m",value=m))
})

observeEvent(inputs$exposure_exposure_select,{
    tmpTab <- parameters$parTab
    if(inputs$typing_flags == 0){
        tmpTab <- tmpTab[tmpTab$id == inputs$exposure_select & tmpTab$exposure == inputs$exposure_exposure_select,]
    } else {
        tmpTab <- tmpTab[tmpTab$type == inputs$exposure_select & tmpTab$exposure == inputs$exposure_exposure_select,]
    }
    tmpTab <- tmpTab[tmpTab$strain == inputs$exposure_strain_select,]
    mu <- tmpTab[tmpTab$names == "mu","values"]
    tp <- tmpTab[tmpTab$names == "tp","values"]
    dp <- tmpTab[tmpTab$names == "dp","values"]
    ts <- tmpTab[tmpTab$names == "ts","values"]
    m <- tmpTab[tmpTab$names == "m","values"]
    isolate(updateSliderInput(session,inputId="mu",value=mu))
    isolate(updateSliderInput(session,inputId="tp",value=tp))
    isolate(updateSliderInput(session,inputId="dp",value=dp))
    isolate(updateSliderInput(session,inputId="ts",value=ts))
    isolate(updateSliderInput(session,inputId="m",value=m))
})
    

    
## When we update the selected exposure ID here, we need to update the parameter entries
observeEvent(inputs$exposure_id,{
    ## Get the exposure table subset related to this ID
    id_choices <- get_ids(parameters)
    tmpTab <- parameters$exposureTab
    tmpTab <- tmpTab[tmpTab$id == inputs$exposure_id,]
    parameters$selectedID <- which(id_choices == inputs$exposure_id)
    
    ## Get available type choices - note that each ID should be only one type!
    type_choices <- get_types(inputs)
    cur_type <- unique(tmpTab$type)[1]
    parameters$selectedType <- which(type_choices == cur_type)

    
    ti <- unique(tmpTab$values)[1] ## Get exposure time - note that each ID should be one exposure time
    exposure_strain <- unique(tmpTab$exposure)[1]  ## Get exposure strain - note that this should also be only one strain
    affected_strains <- unique(tmpTab$strain)  ## Get strains affected
    group <- unique(tmpTab$group) ## Get group
    primed <- unique(tmpTab$primed)            
    
    updateSelectInput(session,inputId="exposure_type",selected=cur_type)
    updateNumericInput(session,inputId="exposure_ti",value=ti)
    updateNumericInput(session,inputId="exposure_group",value=group)
    updateSelectInput(session,inputId="exposure_strain",selected=exposure_strain)
    updateCheckboxGroupInput(session,inputId="exposure_affects",selected=affected_strains)
    updateCheckboxInput(session,"is_primed",value=primed)
})

## Randomise mu
observeEvent(inputs$randomise_mu,{
    new_mu <- runif(1,0,max_mu)
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "mu","values"] <- new_mu)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "mu","values"] <- new_mu)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "mu","values"] <- new_mu)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "mu","values"] <- new_mu)
            }
        }
    }
    isolate(updateSliderInput(session,inputId="mu",value=new_mu))
})
    
## Randomise tp
observeEvent(inputs$randomise_tp,{
    new_tp <- runif(1,0,max_tp)
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "tp","values"] <- new_tp)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "tp","values"] <- new_tp)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "tp","values"] <- new_tp)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "tp","values"] <- new_tp)
            }
        }
    }
    isolate(updateSliderInput(session,inputId="tp",value=new_tp))
})
    
## Randomise mu
observeEvent(inputs$randomise_dp,{
    new_dp <- runif(1,0,1)
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "dp","values"] <- new_dp)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "dp","values"] <- new_dp)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "dp","values"] <- new_dp)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "dp","values"] <- new_dp)
            }
        }
    }
    isolate(updateSliderInput(session,inputId="dp",value=new_dp))
})
    
## Randomise mu
observeEvent(inputs$randomise_ts,{
    new_ts <- runif(1,0,max_ts)
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "ts","values"] <- new_ts)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "ts","values"] <- new_ts)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "ts","values"] <- new_ts)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "ts","values"] <- new_ts)
            }
        }
    }
    isolate(updateSliderInput(session,inputId="ts",value=new_ts))
})
    
## Randomise m
observeEvent(inputs$randomise_m,{
    new_m <- runif(1,min_m,0)
    if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "m","values"] <- new_m)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "m","values"] <- new_m)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$names == "m","values"] <- new_m)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$names == "m","values"] <- new_m)
            }
        }
    }
    isolate(updateSliderInput(session,inputId="m",value=new_m))
})
    
observeEvent(inputs$randomise_all,{
    new_mu <- runif(1,0,max_mu)
    new_dp <- runif(1,0,1)
    new_tp <- runif(1,0,max_tp)
    new_ts <- runif(1,0,max_ts)
    new_m <- runif(1,min_m,0)

    new_values <- c(new_mu, new_tp,new_dp,new_ts,new_m)
if(!is.null(parameters$parTab)){
        if(inputs$cr_flags == 0){
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select,"values"] <- new_values)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select,"values"] <- new_values)
            }
        } else {
            if(inputs$typing_flags == 0){
                isolate(parameters$parTab[parameters$parTab$id == inputs$exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select,"values"] <- new_values)
            } else {
                isolate(parameters$parTab[parameters$parTab$type == inputs$exposure_select &
                                          parameters$parTab$exposure == inputs$exposure_exposure_select &
                                          parameters$parTab$strain == inputs$exposure_strain_select,"values"] <- new_values)
            }
        }
    }
     isolate(updateSliderInput(session,inputId="mu",value=new_mu))
     isolate(updateSliderInput(session,inputId="tp",value=new_tp))
     isolate(updateSliderInput(session,inputId="dp",value=new_dp))
     isolate(updateSliderInput(session,inputId="ts",value=new_ts))
     isolate(updateSliderInput(session,inputId="m",value=new_m))
    
})
