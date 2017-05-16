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
