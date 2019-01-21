model_comparison_analyses <- function(wd, multi_chain=TRUE,
                                      adaptive_period,
                                      parTab, exposureTab,
                                      dat_file, options, typing=TRUE,
                                      times, n=1000, PTchain=FALSE,
                                      individuals=c(3,3,3,3,3),
                                      n1=10000){
    chain <- as.data.frame(load_mcmc_chains(wd,parTab,FALSE,
                                            1,adaptive_period,multi_chain,FALSE,PTchain)[["chain"]])
    cr <- options$cr
    form <- options$form
    
    dat <- read.csv(dat_file)
    dat <- as.matrix(rbind(times, dat))
    rownames(dat) <- NULL
    
    ## Create model solving functions
    f <- create_model_group_func_cpp(parTab,exposureTab,version="model",
                                     form=form,typing = typing,cross_reactivity = cr)
    lik <- create_model_group_func_cpp(parTab,exposureTab,dat=dat,PRIOR_FUNC=NULL,version="posterior",
                                       form=form,typing = typing,cross_reactivity = cr,
                                       individuals=individuals)
    
    ## Generate residuals
    res <- NULL
    mle_res <- generate_mle_residuals(f,chain,times,dat, 3)
    
    ## All calculations
    all_calcs <- NULL
    calcs <- calc_dic_2(lik, chain)
    all_calcs <- c(calcs)
    names(all_calcs) <- names(calcs)
    DIC <- calculate_DIC(chain)
    BIC <- calculate_BIC(chain,parTab,dat[2:nrow(dat),])
    AIC <- calculate_AIC(chain,parTab)
    WAIC <- calculate_WAIC(chain, parTab, dat, f, n)
    loo_estimates <- calculate_loo(chain, parTab, dat, f, n1)
    loo_pareto_k <- loo_estimates[[2]]
    loo_estimates <- loo_estimates[[1]]
    n_pars <- nrow(parTab[parTab$fixed == 0,])

    elpd_loo <- loo_estimates$estimates[1,1]
    elpd_loo_se <- loo_estimates$estimates[1,2]
    p_loo <- loo_estimates$estimates[2,1]
    p_loo_se <- loo_estimates$estimates[2,2]
    looic <- loo_estimates$estimates[3,1]
    looic_se <- loo_estimates$estimates[3,2]
    
    
    return(list(DIC=DIC,BIC=BIC,res=res,calcs=all_calcs,AIC=AIC,n_pars,
                WAIC=WAIC[[1]],pwaic=WAIC[[2]],mle_res=mle_res,
                elpd_loo = elpd_loo, elpd_loo_se = elpd_loo_se,
                p_loo = p_loo, p_loo_se = p_loo_se,
                looic = looic, looic_se = looic_se,
                loo_estimate=loo_estimates,
                pareto_k=loo_pareto_k
                ))
}


generate_mle_residuals <- function(model_func, chain, times,
                                   dat, n_individuals){
    pars <- get_best_pars(chain)
    y <- model_func(pars,times)
    y <- apply(y,2,function(x) rep(x, each=n_individuals))
    y[y > pars[4]] <- pars[4]
    y <- floor(y)

    ## Observed - predicted
    res <- dat[2:nrow(dat),] - y

    ## Collapse predicted vaues to a vector
    y <- c(y)

    ## Collapse residuals to a vector
    res <- c(res)

    ## Return predicted against residuals
    return(cbind(y=y,residual=res))
}
    

generate_residuals <- function(model_func, chain, times,
                               n_samples=1000, n_individuals=3){
    samples <- sample(1:nrow(chain),n_samples,replace=FALSE)
    res <- NULL
    index <- 1
    for(i in samples){
        pars <- zikaInfer::get_index_pars(chain, i)
        y <- f(pars,times)
        y <- apply(y,2,function(x) rep(x, each=n_individuals))
        y[y > pars[4]] <- pars[4]
        y <- floor(y)
        res[index] <- sum((y - dat[2:nrow(dat),])^2)
        index <- index + 1
    }
    return(res)    
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

typing <- TRUE
multi_chain <- FALSE
PTchain <- TRUE
