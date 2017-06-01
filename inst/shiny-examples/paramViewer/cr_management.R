observeEvent(inputs$sigma_value,{
    parameters$crTab[parameters$crTab$names == inputs$type_cr,"values"] <- inputs$sigma_value
})

observeEvent(inputs$type_cr,{
    value <- parameters$crTab[parameters$crTab$names == inputs$type_cr,"values"]
    updateNumericInput(session,"sigma_value",value=value)
})

observeEvent(inputs$n_strains,{
    tmp <- as.data.frame(unique(t(apply(expand.grid(exposure_strains,exposure_strains),1,sort))))
    colnames(tmp) <- c("Strain 1","Strain 2")
    parameters$antigenicDistTab <- data.frame(tmp,"Distance"=0,stringsAsFactors=FALSE)
})
