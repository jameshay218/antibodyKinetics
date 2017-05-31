observeEvent(inputs$sigma_value,{
    parameters$crTab[parameters$crTab$names == inputs$type_cr,"values"] <- inputs$sigma_value
    print(parameters$crTab)
})

observeEvent(inputs$type_cr,{
    value <- parameters$crTab[parameters$crTab$names == inputs$type_cr,"values"]
    updateNumericInput(session,"sigma_value",value=value)
})
