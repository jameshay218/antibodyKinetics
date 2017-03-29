#' @export
posterior_original <- function(params, names,data,PRIOR_FUNC=NULL){
    y <- predict_titres(params[!(names %in% c("S","EA","MAX_TITRE"))], data[,1])
    ln <- obs_likelihood(y[,2],data[,2],params[names %in% c("S","EA","MAX_TITRE")])
    if(!is.null(PRIOR_FUNC)) ln <- ln + PRIOR_FUNC(params, names)
    return(ln)
    
}

#' @export
create_posterior_original <- function(parTab,data,PRIOR_FUNC=NULL){
    names <- parTab$names
    f <- function(pars){
        return(posterior_original(pars,names,data,PRIOR_FUNC))
    }
    return(f)
    
}
