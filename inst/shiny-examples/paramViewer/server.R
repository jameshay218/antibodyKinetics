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
        
        source("exposure_management.R",local=TRUE)

################################
        ## EXPOSURE TABLE MANAGEMENT
################################
        ## On initialisation, make the exposure table NULL
        parameters$exposureTab <- NULL

      
        output$exposure_table_ids <- renderUI({
            out <- NULL
            if(inputs$typing_flags == "0"){
                out <- selectInput("exposure_select","Exposure",
                            unique(parameters$exposureTab$id),
                            selected=1)
            } else {
                types <- NULL
                if(inputs$typing_flags == 0){
                    types <- c("all"="all")
                } else if(inputs$typing_flags == 1){
                    types <- weak_types
                } else {
                    types <- strong_types
                }
                out <- selectInput("exposure_select","Exposure",
                                   types,
                                   selected=1)
            }
            out
        })



        
################################
        ## EXPOSURE TYPE RESTRICTIONS
################################
        ## This is for displaying available exposure types when adding exposures.
        ## The exposures available depends on which analysis is being considered
        get_available_exposure_types <- reactive({
            types <- NULL
            if(inputs$typing_flags == 0){
                types <- c("all"="all")
            } else if(inputs$typing_flags == 1){
                types <- weak_types
            } else {
                types <- strong_types
            }
            return(types)
        })

        ## As above, but relating to sigma parameters. If not using cross reactivity, then no
        ## sigma. Otherwise, depends on if we are using typed cross reactivity of not.
       get_available_exposure_types_cr <- reactive({
            types <- NULL
            if(inputs$cr_flags == 0){
                types <- c("None"="none")
            } else if(inputs$cr_flags == 1){
                types <- c("all"="all")
            } else {
                if(inputs$typing_flags == 0){
                    types <- c("all"="all")
                } else if(inputs$typing_flags == 1){
                    types <- weak_types
                } else {
                    types <- strong_types
                }
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
