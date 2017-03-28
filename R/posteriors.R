#' @export
create_test_posterior <- function(parTab, data, PRIOR_FUNC=NULL){
    names <- parTab$names
    f <- function(pars){
        return(test_posterior(pars, names, data))
    }
    return(f)    
}

#' @export
test_posterior <- function(pars, names, data){
    names(pars) <- names
    mu <- pars["mu"]
    sd <- pars["sd"]

    lnlik <- sum(dnorm(data, mean=mu,sd=sd, log=TRUE))
    lnlik
}


#' @export
posterior_2 <- function(params, names,data){
    y <- predict_titres(params[!(names %in% c("S","EA"))], data[,1])
    ln <- posterior(y[,2],data[,2],params[names %in% c("S","EA")])
    return(ln)
    
}

#' @export
create_test_posterior2 <- function(parTab,data,PRIOR_FUNC=NULL){
    names <- parTab$names
    f <- function(pars){
        return(posterior_2(pars,names,data))
    }
    return(f)
       
}
