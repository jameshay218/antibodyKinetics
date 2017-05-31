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
        parameters$antigenicDistTab <- reactive({
            return(data.frame(t(combn(exposure_strains[1:inputs$n_strains],2,)),values=0,stringsAsFactors=FALSE))
        })
        
        source("exposure_management.R",local=TRUE)
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

    })
