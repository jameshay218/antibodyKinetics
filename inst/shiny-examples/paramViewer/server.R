library(shiny)
library(ggplot2)

options(shiny.maxRequestSize=1000*1024^2)
tmax=100
## Names of all the allowable strains. Can change this to actual strain names
## if preferred
exposure_strains <- c("A","B","C","D","E",
                      "F","G","H","I","J")
shinyServer(
    function(inputs, output, session){
      
################################
        ## STRAIN RESTRICTIONS
################################
        
        ## Gets the number of allowable strains
        n_exposures <- reactive({
            inputs$n_strains
        })
 ## When choosing which exposure strain to use, only display the number of strains we have
        output$choose_exposure_strain <- renderUI({
            selectInput("exposure_strain","Exposure",as.list(exposure_strains[1:n_exposures()]))
        })

        ## Similarly, only display the present strains when determining which titres an exposure
        ## affects
        output$choose_exposure_affects <- renderUI({
            checkboxGroupInput("exposure_affects","Affects",
                               choices=as.list(exposure_strains[1:n_exposures()]),
                               selected=c(1:as.numeric(n_exposures())),
                               inline=TRUE)
        })
################################

################################
        ## EXPOSURE TYPE RESTRICTIONS
################################
        get_available_exposure_types <- reactive({
            types <- NULL
            if(inputs$typing_flags == 0){
                types <- c("All"="all")
            } else if(inputs$typing_flags == 1){
                types <- c("Infection 1"="infection1","Vaccine 1"="vacc1","Vaccine 2"="vacc2",
                           "Adjuvanted 1"="adj1","Adjuvanted 2"="adj2","Infection 2"="infection2")
            } else {
                types <- c("Infection"="infection","Vaccine"="vaccine","Adjuvanted"="adjuvanted")
            }
            return(types)
        })
 
       get_available_exposure_types_cr <- reactive({
            types <- NULL
            if(inputs$cr_flags == 0){
                types <- c("None"="none")
            } else if(inputs$cr_flags == 1){
                types <- c("All"="all")
            } else {
                if(inputs$typing_flags == 0){
                    types <- c("All"="all")
                } else if(inputs$typing_flags == 1){
                    types <- c("Infection 1"="infection1","Vaccine 1"="vacc1","Vaccine 2"="vacc2",
                               "Adjuvanted 1"="adj1","Adjuvanted 2"="adj2","Infection 2"="infection2")
                } else {
                    types <- c("Infection"="infection","Vaccine"="vaccine","Adjuvanted"="adjuvanted")
                }
            }
            return(types)
        })

        output$choose_exposure_type <- renderUI({
            selectInput("type","Type",
                        get_available_exposure_types(),
                        selected=1)
        })
        
        output$choose_exposure_type_cr<- renderUI({
            selectInput("type_cr","CR Type:",
                        get_available_exposure_types_cr(),
                        selected=1)
        })
        
     

## Create exposure entry
        parameters <- reactiveValues()
        parameters$exposureTab <- NULL
        parameters$parTab <- NULL
        
        
        observeEvent(inputs$add_exposure,{
            primed <- 0
            if(inputs$primed) primed <- 1
            exposureRow <- isolate(data.frame(id=inputs$ID,
                                      values=inputs$ti,
                                      type=inputs$type,
                                      exposure=inputs$exposure_strain,
                                      strain=inputs$exposure_affects,
                                      order=1,
                                      primed=primed,
                                      group=inputs$group,
                                      end=tmax,
                                      next_t=tmax))

            parTabRow <- isolate(data.frame(names=c("mu","tp","dp","ts","m"),
                                            id=inputs$ID,
                                            values=c(inputs$mu,inputs$tp,inputs$dp,inputs$ts,inputs$m),
                                            type=inputs$type,
                                            exposure=inputs$exposure_strain,
                                            strain=rep(inputs$exposure_affects,each=5),
                                            order=1,
                                            fixed=1,steps=0.1,lower_bound=0,upper_bound=c(15,20,1,100,1)))
            parameters$exposureTab <- isolate(rbind(parameters$exposureTab, exposureRow))
            parameters$parTab <- isolate(rbind(parameters$parTab, parTabRow))
        })
        
        observeEvent(inputs$remove_exposure,{
            newDat <- isolate(parameters$exposureTab[parameters$exposureTab$id != inputs$exposure_select,])
            newParTab <- isolate(parameters$parTab[parameters$parTab$id != inputs$exposure_select,])
            parameters$exposureTab <- newDat
            parameters$parTab <- newParTab
            
        })
        
        get_available_exposure_ids <- eventReactive(c(inputs$add_exposure,inputs$remove_exposure),{
            tmp <- isolate(unique(parameters$exposureTab$id))
            return(tmp)
        })
        
        output$exposure_table_ids <- renderUI({
            selectInput("exposure_select","Exposure",
                        get_available_exposure_ids(),
                        selected=1)
        })
        
################################
################################
        ## PLOTS
################################
        
        output$main_plot <- renderPlot({
            ## Check if we have the data to make the plot
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
                f <- create_model_group_func_cpp(parTab,exposureTab,form="isolated",cross_reactivity=FALSE,typing=FALSE)
                pars <- parTab$values
                print(pars)
                times <- seq(0,100,by=10)
                dat <- f(pars,times)
                print(dat)

                ## Plotting time
#                ggplot(dat) +
#                    geom_line(aes(x=x,y=y)) +
 #                   geom_vline(data=parameters$exposureTab,xintercept=parameters$exposureTab$values,col="red")
            }
        })
        
    })
