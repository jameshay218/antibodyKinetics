load_mcmc_chains <- function(location="",parTab,unfixed=TRUE, thin=1, burnin=100000, multi=TRUE, chainNo=FALSE, PTchain=FALSE){
    if(multi){
        chains <- Sys.glob(file.path(location,"*multivariate_chain.csv"))
    } else {
        chains <- Sys.glob(file.path(location,"*univariate_chain.csv"))
    }
    if(PTchain){
        chains <- Sys.glob(file.path(location,"*_chain.csv"))
    }
    print(chains)
    if(length(chains) < 1){
        message("Error - no chains found")
        return(NULL)
    }

    ## Read in the MCMC chains with fread for speed
    read_chains <- lapply(chains,data.table::fread,data.table=FALSE)

    ## Thin and remove burn in
    read_chains <- lapply(read_chains, function(x) x[seq(1,nrow(x),by=thin),])
    read_chains <- lapply(read_chains,function(x) x[x$sampno > burnin,])
    print(lapply(read_chains, nrow))
    
    if(chainNo){
        for(i in 1:length(read_chains)) read_chains[[i]]$chain <- i
    }
    
    ## Get the estimated parameters only
    if(unfixed){
        fixed <- parTab$fixed
        read_chains <- lapply(read_chains, function(x) x[,c(which(fixed==0)+1,ncol(x))])
    }

    ## Try to create an MCMC list. This might not work, which is why we have a try catch
    list_chains <- tryCatch({
        tmp_list <- lapply(read_chains,coda::as.mcmc)
        tmp_list <- coda::as.mcmc.list(tmp_list)
    }, warning = function(w){
        print(w)
        NULL
    }, error = function(e){
        print(e)
        NULL
    },
    finally = {
        tmp_list
    })
    
    chain <- as.mcmc(do.call("rbind",read_chains))
    return(list("list"=list_chains,"chain"=chain))
}

summarise_chain <- function(chain){
  tmp <- summary(as.mcmc(chain))
  return(cbind(tmp$statistics,tmp$quantiles))
}

median.quantile <- function(x){
  out <- quantile(x, probs = c(0.025,0.5,0.975))
  names(out) <- c("ymin","y","ymax")
  return(out)
}
