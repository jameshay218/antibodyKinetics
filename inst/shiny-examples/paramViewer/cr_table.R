##################################
## CR Handsontable
##################################
## This bit of code makes the antigenicTab an interactive table that we can edit.
observe({
    if(!is.null(inputs$antigenic_table)){
        parameters[["previous_antigenic"]] <- isolate(parameters[["antigenicDistTab"]])
        antigenicTab = hot_to_r(inputs$antigenic_table)
    } else {
        if (is.null(parameters[["antigenicDistTab"]]))
            antigenicTab <- NULL
        else
            antigenicTab <- parameters[["antigenicDistTab"]]
    }
    parameters[["antigenicDistTab"]] <- antigenicTab
})


output$antigenic_table <- renderRHandsontable({
    antigenicTab <- parameters$antigenicDistTab
    if(!is.null(antigenicTab)){
        rhandsontable(antigenicTab, useTypes = FALSE, stretchH="all")
    }})

output$export_cr <- downloadHandler(
    filename = "antigenicDistTab.csv",
    content=function(file){
        write.csv(parameters$antigenicDistTab,file,row.names=FALSE)
    })

observeEvent(inputs$upload_antigenic_distances,{
    newTab <- read.csv(isolate(inputs$antigenic_tab_input$datapath),stringsAsFactors=FALSE)
    if(!is.null(newTab)){
        updateNumericInput(session,"n_strains",value=length(unique(newTab[,1])))
        parameters$antigenicDistTab <- newTab
    }
})
