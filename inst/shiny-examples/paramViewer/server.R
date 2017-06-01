library(shiny)
source("helpers.R",local=TRUE)
source("global_parameters.R",local=TRUE)
library(rhandsontable)
library(ggplot2)

#options(shiny.maxRequestSize=1000*1024^2)
shinyServer(
    function(inputs, output, session){
        parameters <- reactiveValues()
        ## Keep track of the currently selected 
        parameters$selectedID <- 1
        parameters$selectedType <- 1
        parameters$exposureTab <- NULL
        parameters$parTab <- NULL
        parameters$crTab <- data.frame(names=c("None","all",weak_types,strong_types),values=-Inf,stringsAsFactors=FALSE)
        tmp <- as.data.frame(unique(t(apply(expand.grid(exposure_strains,exposure_strains),1,sort))))
        colnames(tmp) <- c("Strain 1","Strain 2")
        parameters$antigenicDistTab <- data.frame(tmp,"Distance"=0,stringsAsFactors=FALSE)
        
        
        source("exposure_management.R",local=TRUE)
        source("cr_table.R",local=TRUE)
        source("cr_management.R",local=TRUE)
################################
        ## EXPOSURE TABLE MANAGEMENT
################################
        ## On initialisation, make the exposure table NULL
        output$exposure_table_id <- renderUI({
            out <- NULL
            if(inputs$typing_flags == 0){
                ids <- get_ids(parameters)
                out <- selectInput("exposure_select","Exposure",
                                   ids)
            } else {
                out <- "none"
                if(!is.null(parameters$parTab)){
                    #types <- get_types(inputs)
                    types <- unique(parameters$parTab$type)
                    out <- selectInput("exposure_select","Exposure",
                                       types)
                }
            }
            out
        })

        
################################
        ## EXPOSURE TYPE RESTRICTIONS
################################
        ## This is for displaying available exposure types when adding exposures.
        ## The exposures available depends on which analysis is being considered
        get_available_exposure_types <- reactive({get_types(inputs)})

        ## As above, but relating to sigma parameters. If not using cross reactivity, then no
        ## sigma. Otherwise, depends on if we are using typed cross reactivity of not.
        get_available_exposure_types_cr <- reactive({
            types <- NULL
            if(inputs$cr_flags == 0){
                types <- c("None"="none")
            } else if(inputs$cr_flags == 1){
                types <- c("all"="all")
            } else {
                types <- get_types(inputs)
            }
            return(types)
        })
        
                                        # Related to the above, display the output select input based on
        ## availability
        output$choose_exposure_type_cr<- renderUI({
            selectInput("type_cr","CR Type:",
                        get_available_exposure_types_cr(),
                        selected=1)
        })
################################
        source("exposure_table.R",local=TRUE)
        source("exposure_table_management.R",local=TRUE)
        source("plots.R",local=TRUE)

        output$download_all <- downloadHandler(
            filename="parTab.csv",
            
            content=function(file){
                top_parTab <- data.frame(names=c("lower_bound","S","EA","MAX_TITRE"), id="all",
                                         values=c(inputs$lower_bound,0.79,0.2,inputs$max_titre),
                                         type="all",
                                         exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                         lower_bound=c(-1000,0,0,0),upper_bound=c(0,1,1,100),stringsAsFactors=FALSE)

                if(inputs$cr_flags != 0){
                    tmpCrTab <- parameters$crTab[parameters$crTab$names %in% get_available_exposure_types_cr(),]
                    cr_values <- tmpCrTab$values
                    cr_names <- tmpCrTab$names
                    bot_parTab <- data.frame(names=c("beta","c",rep("sigma",length(cr_names)),"y0_mod"),id="all",
                                             values=c(inputs$beta,inputs$c,cr_values,inputs$y0_mod),
                                             type=c("all","all",cr_names,"all"),
                                             exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                             lower_bound=c(-20,0,-20,-20),upper_bound=c(2,20,2,2),stringsAsFactors=FALSE)
                } else {
                    bot_parTab <- data.frame(names=c("beta","c","sigma","y0_mod"),id="all",
                                             values=c(inputs$beta,inputs$c,-Inf,inputs$y0_mod),
                                             type=c("all","all","all","all"),
                                             exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                             lower_bound=c(-20,0,-20,-20),upper_bound=c(2,20,2,2),stringsAsFactors=FALSE)
                }

                mod_parTab <- data.frame(names="mod",id=NA,values=c(inputs$mod1,inputs$mod2,inputs$mod3,inputs$mod4),
                                         type="all",exposure=NA,strain=NA,order=NA,fixed=1,steps=0.1,
                                         lower_bound=0,upper_bound=1,stringsAsFactors=FALSE)

                distance_parTab <- data.frame(names="x",id=NA,values=parameters$antigenicDistTab$Distance,
                                              type="all",exposure=parameters$antigenicDistTab$Strain.1,
                                              strain=parameters$antigenicDistTab$Strain.2,
                                              order=NA,fixed=1,steps=0.1,lower_bound=0,upper_bound=10000,
                                              stringsAsFactors=FALSE)
                print(colnames(top_parTab))
                print(colnames(parameters$parTab))
                print(colnames(bot_parTab))
                print(colnames(distance_parTab))
                print(colnames(mod_parTab))
                tmpTab <- parameters$parTab
                tmpTab[tmpTab$names == "m","values"] <- exp(tmpTab[tmpTab$names == "m","values"])
                parTab <- rbind(top_parTab,tmpTab,bot_parTab,distance_parTab,mod_parTab)
                
                write.csv(parTab,file,row.names=FALSE)
            }
        )
    })


