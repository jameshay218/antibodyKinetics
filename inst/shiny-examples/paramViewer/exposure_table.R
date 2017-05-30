##################################
## Exposure Handsontable
##################################
## This bit of code makes the exposureTab an interactive table that we can edit.
observe({
    if(!is.null(inputs$exposure_table)){
        parameters[["previous"]] <- isolate(parameters[["exposureTab"]])
        exposureTab = hot_to_r(inputs$exposure_table)
    } else {
        if (is.null(parameters[["exposureTab"]]))
            exposureTab <- NULL
        else
            exposureTab <- parameters[["exposureTab"]]
    }
    parameters[["exposureTab"]] <- exposureTab
})

output$exposure_table <- renderRHandsontable({
    exposureTab <- parameters$exposureTab
    if(!is.null(exposureTab)){
        rhandsontable(exposureTab, useTypes = FALSE, stretchH="all")
    }})

output$export_exposures <- downloadHandler(
    filename = "exposureTab.csv",
    content=function(file){
        write.csv(parameters$exposureTab,file,row.names=FALSE)
    })

observeEvent(inputs$upload_exposures,{
    newTab <- read.csv(isolate(inputs$exposure_tab_input$datapath),stringsAsFactors=FALSE)
    
    ## If the new parameter table has types in it, need to set the flags etc correctly
    if(any(strong_types%in% newTab$type)){
        updateSelectInput(session,"typing_flags",selected=2)
    } else if(any(weak_types %in% newTab$type)){
        updateSelectInput(session,"typing_flags",selected=1)
    } else {
        updateSelectInput(session,"typing_flags",selected=0)
    }
    updateNumericInput(session,"tmax",value=max(newTab$end))
    updateNumericInput(session,"n_strains",value=length(unique((newTab$strain))))
    parameters$exposureTab <- newTab

    ## THEN NEED TO UPDATE THE OVERALL PARAMETER TABLE
    
})
