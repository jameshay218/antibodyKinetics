#' Deviance
#'
#' Calculates deviance of a vector with a given likelihood function
#' @param x a vector of parameters used to solve the likelihood function
#' @param likelihood a likelihood function taking just a vector of parameters
#' @return a single deviance value
#' @export
calc_deviance  <- function(x, likelihood){
    return(-2*(likelihood(x)))
}

#' Posterior mean
#'
#' Calculates posterior mean of a given MCMC chain
#' @param chain the MCMC chain with a lnlike column
#' @return the posterior mean
#' @export
calc_post_mean <- function(chain){
    return(-2*mean(chain$lnlike))
}

#' Posterior variance
#'
#' Calculates posterior variance from a given MCMC chain
#' @param chain the MCMC chain with a lnlike column
#' @return the posterior variance
#' @export
calc_post_var <- function(chain){
    meanBar <- calc_post_mean(chain)
    tmp <- 0
    for(i in 1:nrow(chain)){
        tmp <- tmp + (-2*chain[i,"lnlike"] - meanBar)^2
    }
  varBar <- (1/(nrow(chain)-1)) * sum(tmp)
  return(varBar)
}

#' DIC
#'
#' Calculates DIC from a given MCMC chain
#' @param chain the MCMC chain with a lnlike column and all columns (number of columns will be used
#' @return a single DIC value
#' @export
calculate_DIC <- function(chain){
    DIC <- calc_post_mean(chain) + calc_post_var(chain)/2
    return(DIC)
}

#' BIC calculation
#' 
#' Given an MCMC chain, calculates the BIC
#' @param chain the MCMC chain to be tested
#' @param parTab the parameter table used for this chain
#' @param dat the data
#' @return a single BIC value
#' @export
calculate_BIC <- function(chain, parTab, dat){
    n <- nrow(dat)*ncol(dat)
    maxLik <- -2*max(chain$lnlike)
    B <- length(parTab[parTab$fixed==0,"values"])*log(n)
    return(maxLik + B)
}


#' @export
calc_dic_2 <- function(lik.fun,chain){
    D.bar <- -2*mean(chain$lnlike)
    theta.bar <- as.numeric(summary(as.mcmc(chain[,2:(ncol(chain)-1)]))$statistics[,"Mean"])
    D.hat <- -2*lik.fun(theta.bar)
    pD <- D.bar - D.hat
    pV <- var(-2*chain$lnlike)/2
    list(DIC=2*pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}

#' @export
calculate_AIC <- function(chain, parTab){
    k <- nrow(parTab[parTab$fixed == 0,])
    AIC <- 2*k - 2*(max(chain$lnlike))
    return(AIC)
}


#' WAIC using the pwaic2 in Gelman 2013
#'
#' @export
calculate_WAIC <- function(chain, parTab, dat, f, N){
    ## Do this for a subsample of the MCMC chain to save time
    samps <- sample(1:nrow(chain),size = N,replace=FALSE)
    chain <- chain[samps,]

    ## Separate out data - times and measurements
    times <- dat[1,]
    dat <- dat[2:nrow(dat),]
    
    expectation_likelihood <- 0
    tmp <- matrix(nrow=nrow(chain),ncol=nrow(dat)*ncol(dat))

    ## For each sample
    for(i in 1:nrow(chain)){
        pars <- zikaProj::get_index_pars(chain,i)

        ## Get parameters
        y <- f(pars,times)
        index <- 1

        ## For each strain/group pairing
        for(j in 1:nrow(y)){
            ## For each time point
            for(x in 1:ncol(y)){
                ## We have 3 of each individual
                for(q in 1:3){
                    ## Get likelihood for this point
                    wow <- norm_error(y[j,x],dat[(3*(j-1)+1)+(q-1),x],pars["S"],pars["MAX_TITRE"])
                    ## expectation posterior
                    expectation_likelihood <- expectation_likelihood + log(wow)

                    ## Store this likelihood
                    tmp[i, index] <- wow
                    index <- index + 1                
                }
            }
        }
    }

    ## Get mean of log likelihoods
    expectation_likelihood <- expectation_likelihood/nrow(chain)

    ## Get variance of all log likelihoods for each data point, y
    vars <- apply(log(tmp),2,var)

    ## lppd is the log posterior predictive density
    ## Taken by taking the log of the mean of each set of likelihoods
    ## for each observation y, and then sum this.
    lppd <- sum(log(colMeans(tmp)))

    ## PWAIC is measure of the number of effective free parameters in the model
    ## It is the sum of the variances of log likelihoods
    pwaic <- sum(vars)

    ## Alternative way of measuring PWAIC
    ## a <- log(apply(tmp,2, mean))
    ## b <- apply(log(tmp),2,mean)
    #pwaic1 <- 2*(a - b)
    
    return(list(WAIC=-2*(lppd-pwaic),pwaic=pwaic))
}
