## Function to display the available exposure types depending on the selected type option
get_types <- function(inputs){
    types <- NULL
    if(inputs$typing_flags == 0){
        types <- c("all"="all")
    } else if(inputs$typing_flags == 1){
        types <- weak_types
    } else {
        types <- strong_types
    }
    return(types)
}

## Display available exposure IDs
get_ids <- function(parameters){
    unique(parameters$exposureTab$id)
}

## Correctly update exposure table when new exposures added
add_order_nextt<- function(inputs, parameters){
    ## get exposures at same time, before or after
    next_exposure_indices <- which(parameters$exposureTab$values > inputs$exposure_ti &
                                    parameters$exposureTab$group == inputs$exposure_group)
    prev_exposure_indices <- which(parameters$exposureTab$values < inputs$exposure_ti &
                                    parameters$exposureTab$group == inputs$exposure_group)
    same_exposure_indices <- which(parameters$exposureTab$values == inputs$exposure_ti &
                                   parameters$exposureTab$group == inputs$exposure_group)
    
    next_infection <- parameters$exposureTab[next_exposure_indices,]
    previous_infection <- parameters$exposureTab[prev_exposure_indices,]
    same_infection <- parameters$exposureTab[same_exposure_indices,]

    tmpTab <- parameters$exposureTab
    
    ## If instead we get a data frame of 0 rows, just set to NULL
    if(!is.null(next_infection) && nrow(next_infection)== 0) next_infection <- NULL
    if(!is.null(previous_infection) && nrow(previous_infection)== 0) previous_infection <- NULL
    if(!is.null(same_infection) && nrow(same_infection)== 0) same_infection <- NULL
    
    ## If this is the only exposure, this is number 1
    if(is.null(same_infection) & is.null(next_infection) & is.null(previous_infection)){
        ## Just add in naively
        new_order <- 1
        next_t <- inputs$tmax
        new_id <- paste0("G",inputs$exposure_group,"E",new_order,"S",inputs$exposure_strain)
        return(list("next_t"=next_t,"order"=new_order,"new_id"=new_id,"newTab"=tmpTab))
    }

    ## Otherwise, there are other exposures. If this is at the same time as another
    ## exposure, then just use the same properties.
    if(!is.null(same_infection)){
        new_order <- same_infection[1,"order"]
        next_t <- same_infection[1,"next_t"]
        new_id <- paste0("G",inputs$exposure_group,"E",new_order,"S",inputs$exposure_strain)
        return(list("next_t"=next_t,"order"=new_order,"new_id"=new_id,"newTab"=tmpTab))
    }

    ## Otherwise, there might be exposures before or after
    ## If this is the first exposure, then set this order to 1 and increase
    ## the other orders by 1. Next t is next exposure time
    if(is.null(previous_infection)){
        new_order <- 1
        next_t <- inputs$tmax
        ## If there are exposures afterwards
        if(!is.null(next_infection)){
            next_t <- next_infection[1,"values"]
            tmpTab[next_exposure_indices,"order"] <- tmpTab[next_exposure_indices,"order"] + 1
            for(i in next_exposure_indices) tmpTab[i,"id"] <- paste0("G",tmpTab[i,"group"],"E",tmpTab[i,"order"],"S",tmpTab[i,"exposure"])
        }
        new_id <- paste0("G",inputs$exposure_group,"E",new_order,"S",inputs$exposure_strain)
    } else {
        ## Otherwise, if there are previous exposure, new order is max of previous exposure
        ## orders +1. And previous exposure "next_t" is set to current exposure time       
        new_order <- max(previous_infection$order) + 1
        next_t <- inputs$tmax

        tmpTab[tmpTab$values < inputs$exposure_ti &
               tmpTab$order == (new_order - 1) &
               tmpTab$group == inputs$exposure_group,
               "next_t"] <- inputs$exposure_ti
        ## If there are exposures afterwards
        if(!is.null(next_infection)){
            next_t <- next_infection[1,"values"]
            tmpTab[next_exposure_indices,"order"] <- tmpTab[next_exposure_indices,"order"] + 1
            for(i in next_exposure_indices) tmpTab[i,"id"] <- paste0("G",tmpTab[i,"group"],"E",tmpTab[i,"order"],"S",tmpTab[i,"exposure"])
        }
        new_id <- paste0("G",inputs$exposure_group,"E",new_order,"S",inputs$exposure_strain)
    }
    return(list("next_t"=next_t,"order"=new_order,"new_id"=new_id,"newTab"=tmpTab))
}


## Correctly update exposure table when new exposures removed
remove_order_nextt<- function(inputs, parameters){
    ## get exposures at same time, before or after
    removed_exposure <- parameters$exposureTab[parameters$exposureTab$id == inputs$exposure_id,]
    ti <- removed_exposure[1,"values"]
    order <- removed_exposure[1,"order"]
    ## Get all exposures apart from the selected one

    tmpTab <- parameters$exposureTab[!(parameters$exposureTab$id == inputs$exposure_id),]
    next_exposure_indices <- which(tmpTab$values > ti &
                                     tmpTab$group == inputs$exposure_group)
    prev_exposure_indices <- which(tmpTab$values < ti &
                                   tmpTab$group == inputs$exposure_group)
    
    next_infection <- tmpTab[next_exposure_indices,]
    previous_infection <- tmpTab[prev_exposure_indices,]
    
    ## If instead we get a data frame of 0 rows, just set to NULL
    if(!is.null(next_infection) && nrow(next_infection)== 0) next_infection <- NULL
    if(!is.null(previous_infection) && nrow(previous_infection)== 0) previous_infection <- NULL
    
    ## If this is the only exposure, exposure table is now empty
    if(is.null(next_infection) & is.null(previous_infection)){
        ## Just remove naively
        return(tmpTab)
    }
    ## If there are subsequent exposures, need to reduce their order and ID
    next_t <- inputs$tmax
    if(!is.null(next_infection)){
        tmpTab[next_exposure_indices,"order"] <- tmpTab[next_exposure_indices,"order"] - 1
        for(i in next_exposure_indices) tmpTab[i,"id"] <- paste0("G",tmpTab[i,"group"],"E",tmpTab[i,"order"],"S",tmpTab[i,"exposure"])
        next_t <- tmpTab[next_exposure_indices[1],"values"]
    }

    ## If there were previous exposures, need to change their next_t
    if(!is.null(previous_infection)){
        tmpTab[tmpTab$values < inputs$exposure_ti &
               tmpTab$order == (order - 1) &
               tmpTab$group == inputs$exposure_group,
               "next_t"] <- next_t
    }
    return(tmpTab)
}


