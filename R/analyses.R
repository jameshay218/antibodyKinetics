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
#' Calculates WAIC using the second method described in Gelman 2013
#' @param chain the MCMC chain with each row corresponding to an iteration
#' @param parTab the parameter table used for this chain
#' @param dat the data frame of data used to solve the likelihood
#' @param f pointer to the model function used in the MCMC procedure
#' @param N number of samples to take from the MCMC chain
#' @return list with WAIC for this model calculated from the MCMC chain, and pwaic which is a constituent of WAIC
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

    data("ferret_titres")
    
   ## For each sample
    for(i in 1:nrow(chain)){
        pars <- zikaInfer::get_index_pars(chain,i)

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

#' Calculate per data point likelihood matrix
#'
#' From the output of an MCMC chain, calculates the likelihood of observing each data point. This can then be used in loo::loo
#' @inheritParams calculate_WAIC
#' @param nindiv the number of individuals in each group/strain combination
#' @return a matrix of same dimensions as dat
#' @export
calculate_lik_matrix <- function(chain, parTab, dat, f, N, nindiv=3){
    ## Do this for a subsample of the MCMC chain to save time
    samps <- sample(1:nrow(chain),size = N,replace=FALSE)
    chain <- chain[samps,]

    ## Separate out data - times and measurements
    times <- dat[1,]
    dat <- dat[2:nrow(dat),]
    
    tmp <- matrix(nrow=nrow(chain),ncol=nrow(dat)*ncol(dat))
    data(ferret_titres)

    groups <- strains <- character(nrow(dat)*ncol(dat))
    all_times <- numeric(nrow(dat)*ncol(dat))
    indiv <- character(nrow(dat)*ncol(dat))
    times <- c(0,21,37,49,70)
    data_points <- numeric(nrow(dat)*ncol(dat))
 
    index <- 1
    ## For each strain/group pairing
    i <- 1
    pars <- zikaInfer::get_index_pars(chain,i)

    ## Get parameters
    y <- f(pars,times)
    
    for(j in 1:nrow(y)){
        ## For each time point
        for(x in 1:ncol(y)){
            ## We have 3 of each individual
            for(q in 1:nindiv){
                groups[index] <- as.character(ferret_titres[(3*(j-1)+1)+(q-1),"group"])
                strains[index] <- as.character(ferret_titres[(3*(j-1)+1)+(q-1),"strain"])
                indiv[index] <- as.character(ferret_titres[(3*(j-1)+1)+(q-1),"indiv"])
                all_times[index] <- times[x]
                data_points[index] <- dat[(3*(j-1)+1)+(q-1),x]
                index <- index + 1
            }
        }
    }
    ## For each sample
    for(i in 1:nrow(chain)){
        pars <- zikaInfer::get_index_pars(chain,i)

        ## Get parameters
        y <- f(pars,times)
        index <- 1

        ## For each strain/group pairing
        for(j in 1:nrow(y)){
            ## For each time point
            for(x in 1:ncol(y)){
                ## We have 3 of each individual
                for(q in 1:nindiv){                   
                    ## Get likelihood for this point
                    wow <- norm_error(y[j,x],dat[(3*(j-1)+1)+(q-1),x],pars["S"],pars["MAX_TITRE"])
                    ## Store this likelihood
                    tmp[i, index] <- wow
                    index <- index + 1                
                }
            }
        }
    }
    all_labels <- data.frame(groups,strains,indiv,all_times, data_points)
    return(list(tmp, all_labels))
}



#' Calculate LOOIC using loo package
#'
#' Calculates PSIS-LOO CV as in Vehtari et al. 2017a/b, using the loo package.
#' @inheritParams calculate_lik_matrix
#' @export
calculate_loo <- function(chain, parTab, dat, f, N=1000, nindiv=3){
  samps <- sample(1:nrow(chain),size = N,replace=FALSE)
  chain <- chain[samps,]
  tmp <- calculate_lik_matrix(chain, parTab, dat, f, N,nindiv=nindiv)
  all_liks <- log(tmp[[1]])
  labels <- tmp[[2]]
  rel_n_eff <- loo::relative_eff(exp(all_liks),chain_id=rep(1,each=N))
  loo1 <- loo::loo(all_liks,r_eff=rel_n_eff,cores=1)
  pareto_k <- data.frame(labels, pareto_k=loo1$diagnostics$pareto_k)
  return(list(loo1,pareto_k))
}

#' Get titre data labels
#'
#' @export
get_titre_labels <- function(nindiv=3){
    data("ferret_titres")
    groups <- NULL
    strains <- NULL
    all_times <- NULL
    times <- c(0,21,37,49,70)
    index <- 1
    for(j in 1:nrow(ferret_titres)){
        for(x in 1:length(times)){
            groups <- c(groups, as.character(ferret_titres[index,"group"]))
            strains <- c(strains, as.character(ferret_titres[index,"strain"]))
            all_times<- c(all_times, times[x])
        }
    }
    all_labels <- data.frame(groups,strains,all_times)
    return(all_labels)
}

