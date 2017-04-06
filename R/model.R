#' Antigenic distance calc
#'
#' Calculates the necessary antigenic distance needed to give a given
#' percentage cross reactivity with gradient m and percentage y
#' @param y the desired proportion cross reactivity
#' @param m the rate of the exponential function for cross reactivity
#' @return a single value, x#'
#' @export
#' @useDynLib antibodyKinetics
calculate_x <- function(y, m){
  return(-log(y)/m)
}

#' Model trajectory calc
#'
#' Calculates the ferret model trajectory for a single infection event.
#' Uses an R implementation so easy to code and has named parameter vectors.
#' @param pars the named vector of model parameters
#' @param times the vector of time in days to solve over
#' @param logSigma if TRUE, uses the exponential of pars["sigma"]
#' @return a vector of antibody titres
#' @export
#' @family model functions
#' @useDynLib antibodyKinetics
#' @examples
#' pars <- c("mu"=8,"dp"=0.5,"tp"=12,"ts"=10,"m"=0.003,"y0"=0,"t_i"=10,
#'           "sigma"=0.01, "beta"=0.02,"c"=4,"x"=0,"primed"=0,"mod"=1)
#' times <- seq(0,100,by=10)
#' y <- model_trajectory(pars,times,FALSE)
model_trajectory <- function(pars, times, logSigma=TRUE){
    ## Calculate modified mu from cross reactivity and modifiers
    mu <- pars["mu"]
    dp <- pars["dp"]
    tp <- pars["tp"]
    ts <- pars["ts"]
    m <- pars["m"]
    y0 <- pars["y0"]
    t_i <- pars["t_i"]

    if(logSigma){
        sigma <- exp(pars["sigma"])
        beta <- exp(pars["beta"])
    } else {
        sigma <-  pars["sigma"]
        beta <- pars["beta"]
    }
    
    cr <- exp(-sigma*pars["x"])
    prime_cr <- pars["c"]*exp(-beta*pars["x"])*pars["primed"]

    mu <- mu*cr*pars["mod"] + prime_cr
    
    y <- numeric(length(times))
    i <- 1
    
    ## Loops through all times and calculate titre based on time relative to time of infection
    for(t in times){
        y[i] <- (
            (t <= t_i)*0
        ) +
            (
                (t > t_i)*(t <= (tp + t_i))*((mu/tp)*t-(mu/tp)*t_i)
            ) +
            (
                (t > (tp + t_i))*(t <= (ts + tp + t_i))*((-(dp*mu)/ts)*t+((mu*dp)/ts)*(t_i+tp) + mu)
            ) +
            (
                (t > (ts + tp + t_i))*(-m*t+m*(t_i+tp+ts)+(1-dp)*mu)
            ) + y0
        i <- i + 1
    }
    return(y)
}

#' Multiple group model trajectory
#'
#' Calculates model trajectories for multiple groups
#' @param parTab the parameter table containing at least a values
#' column and a names column which can be used to solve \code{\link{model_trajectory}}
#' @param cr_table the cross reactivity part of the parameter table, with values and names for the exposure strain, the measured strain, and the antigenic distance (x)
#' @param order_tab the table for parameters modifying boosting based on infection order (antigenic seniority)
#' @param exposures the table for exposure types and times
#' @param strains a vector with the names of all of the strains involved in the model
#' @param times a vector of times to solve the model over
#' @return a table of antibody titres for the given times, with a column for times, group and colnames of the measured strain
#' @export
#' @useDynLib antibodyKinetics
#' @family model functions
#' @seealso \code{\link{model_func}} for next level, \code{link{model_trajectory}} for model solving code
model_func_groups <- function(parTab, cr_table, order_tab, exposures, strains, times){
    ## Get unique groups (groups of exposures)
    groups <- unique(exposures$group)

    ## Precompute matrix to save antibody titres
    y <- matrix(nrow=length(times)*length(groups),ncol=length(strains))
    colnames(y) <- strains
    y <- data.frame(times=rep(times, length(groups)),y, group=rep(groups,each=length(times)))
    ## For each group, isolate group specific exposures and solve model
    for(group in groups){
        tmpExposures <- exposures[exposures$group == group,] 
        y[y$group==group,strains] <- model_func(parTab, cr_table, order_tab, tmpExposures, strains, times)
    }
    return(y)
}

#' Model solver for one group
#'
#' Solves the antibody kinetics model for a single group, which may have multiple measured and exposure strains
#' @param parTab the parameter table containing at least a values
#' column and a names column which can be used to solve \code{\link{model_trajectory}}
#' @param cr_table the cross reactivity part of the parameter table, with values and names for the exposure strain, the measured strain, and the antigenic distance (x)
#' @param order_tab the table for parameters modifying boosting based on infection order (antigenic seniority)
#' @param exposures the table for exposure types and times
#' @param strains a vector with the names of all of the strains involved in the model
#' @param times a vector of times to solve the model over
#' @return a table of antibody titres for the given times, with a column for times and colnames of the measured strain
#' @export
#' @useDynLib antibodyKinetics
#' @family model functions
model_func <- function(parTab, cr_table, order_tab, exposures, strains, times){
    trajectories <- matrix(0, ncol=length(strains),nrow=length(times))
    index <- 1
    ## For each strain to be measured
    for(strain in strains){

        ## For each exposure
        for(i in 1:nrow(exposures)){

            ## Get exposure time, type and exposure strain
            t_i <- exposures[i,"values"]
            type <- exposures[i,"type"]
            exposure <- exposures[i,"exposure"]

            ## Is this first, second etc exposure?
            order <- exposures[i,"order"]

            ## Was there priming?
            isPrimed <- exposures[i,"primed"]

            ## Get model parameters for this type of exposure
            tmpParTab <- parTab[parTab$type %in% c(type,"all","priming"),]
            pars <- tmpParTab$values
            names(pars) <- tmpParTab$names

            ## Get the correct modifier for this exposure order
            mod <- order_tab[order_tab$order == order,"values"]
            ind <- sort(c(exposure,strain))

            ## Find the antigenic distance between the measured and exposure strain
            x <- cr_table[cr_table$exposure==ind[1] & cr_table$strain == ind[2],"values"]

            ## Solve model
            pars <- c("t_i"=t_i,pars,"mod"=mod,"x"=x,y0=0,"primed"=isPrimed)
            y <- model_trajectory(pars, times)
            trajectories[,index] <-  trajectories[,index] +y
            
        }
        index <- index + 1
    }
  colnames(trajectories) <- strains
  return(trajectories)
}
