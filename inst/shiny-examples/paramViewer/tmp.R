makePlot <- !is.null(parameters$exposureTab) && nrow(parameters$exposureTab) > 0
            
            if(makePlot){
                ## Update values for the currently selected exposure
                values <- c(inputs$mu,inputs$tp,inputs$dp,inputs$ts,inputs$m)
                parameters$parTab[parameters$parTab$id == inputs$exposure_select,"values"] <- values
                tmpTab <- parameters$parTab
                exposureTab <- parameters$exposureTab
                overallPars <- data.frame(names=c("lower_bound","S","EA","MAX_TITRE"),id="all",
                                          values=c(inputs$lower_bound,inputs$S,inputs$EA,inputs$max_titre),
                                          type="all",
                                          exposure=NA,
                                          strain=NA,
                                          order=NA,
                                          fixed=1,
                                          steps=0.1,
                                          lower_bound=c(-1000,0,0,0),
                                          upper_bound=c(0,1,1,15))
                cr_pars <- data.frame(names=c("beta","c","sigma","y0_mod"),id="all",
                                      values=c(inputs$beta,inputs$c,-Inf,-20),
                                      type="all",
                                      exposure=NA,
                                      strain=NA,
                                      order=NA,
                                      fixed=1,
                                      steps=0.1,
                                      lower_bound=c(-20,0,-20,-20),
                                      upper_bound=c(2,20,2,2))
                mod_pars <- data.frame(names="mod",id=NA,
                                       values=c(inputs$mod1,inputs$mod2,inputs$mod3,inputs$mod4),
                                       type=NA,
                                       exposure=NA,
                                       strain=NA,
                                       order=c(1,2,3,4),
                                       fixed=1, steps=0.1,lower_bound=0,upper_bound=1)
                x_pars <- data.frame(names="x",id=NA,
                                     values=c(0,50,60,1000,1000,0,10,1000,1000,0,1000,1000,0,200,0),
                                     type=NA, exposure=c(rep("A",5),rep("B",4),rep("C",3),rep("D",2),"E"),
                                     strain=c("A","B","C","D","E","B","C","D","E","C","D","E","D","E","E"),
                                     order=NA,fixed=1,steps=0.1,lower_bound=0,upper_bound=2000)
                parTab <- rbind(overallPars,tmpTab,cr_pars,x_pars,mod_pars)
                parTab[parTab$names=="m","values"] <- exp(parTab[parTab$names=="m","values"])
                print(exposureTab)
                print(lapply(exposureTab,class))
            }
