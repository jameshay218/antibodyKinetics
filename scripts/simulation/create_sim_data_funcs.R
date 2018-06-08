create_data <- function(runName, parTab_file,
                     exposureTab_file, form,typing,cr,priming,
                     ngroup=5,nstrain=5,nindiv=3,
                     times,wd="~/net/home/ferret/inputs/data_normal",
                     group=NA,strain=NA,
                     normal=TRUE){
    parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
    exposureTab <- read.csv(exposureTab_file,stringsAsFactors=FALSE)
    parTab[parTab$names == "MAX_TITRE","values"] <- 13
    parTab[parTab$names == "tp","values"] <- 12
    parTab[parTab$names == "S","values"] <- 1.5
    if(!is.na(group)){
        exposureTab <- exposureTab[exposureTab$group == group,]
        if(!is.na(strain)){
            exposureTab <- exposureTab[exposureTab$strain == strain,]
        }
        
        parTab <- parTab[parTab$id %in% c(NA,"all",unique(exposureTab$id)) | parTab$names == "x",]
    }
    
    ## Simulate data and save
    individuals <- rep(nindiv,ngroup)
    pars <- parTab[parTab$names %in% c("S","EA","MAX_TITRE"),"values"]
    names(pars) <- c("S","EA","MAX_TITRE")
    print(pars)
                                        #pars <- c("S"=0.79,"EA"=0.2,"MAX_TITRE"=)
    f <- create_model_group_func_cpp(parTab,exposureTab,version="model",form=form,typing = typing,cross_reactivity = cr)
    dat <- f(parTab$values, times)
    dat <- floor(dat)
    dat <- apply(dat,2,function(x) rep(x, each=nindiv))
    index <- 1
    difs <- NULL
    for(i in 1:nrow(dat)){
        for(j in 1:ncol(dat)){
            tmp <- dat[i,j]
            dat[i,j] <- add_noise(pars,dat[i,j],normal)
            tmp <- tmp - dat[i,j]
            difs[index] <- tmp
            index <- index + 1
        }
    }
    #dat <- rbind(times, dat)
    rownames(dat) <- NULL

    ## Filenames
    filename <- paste0(wd,"/",runName,"_data.csv")
    if(!is.na(group)){
        filename <- paste0(wd,"/",runName,"_",group,"_data.csv")
    }
    write.csv(dat, filename,row.names=FALSE)
    return(difs)
}

