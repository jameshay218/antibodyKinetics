library(shiny)
library(shinyBS)
library(ggplot2)
library(gridExtra)

tmax=100

shinyUI(    
    navbarPage("Antibody Kinetics Model",
               #' Parameter exploration panel
               tabPanel("Parameters",
                        sidebarPanel(
                            h4(strong("Main parameters")),
                            fluidRow(
                                column(3, numericInput("lower_bound","Lower bound",0,min=-1000,max=0)),
                                column(3, numericInput("S","S",0.79,min=0,max=1)),
                                column(3,numericInput("EA","EA",0.2,min=0,max=1)),
                                column(3,numericInput("max_titre","Max Titre",13,min=0,max=25)),
                                column(3,numericInput("n_strains","No. strains",5,min=1,max=10)),
                                column(3,numericInput("tmax","Max time",100,min=10,max=1000)
                                       )),
                            hr(),
                            fluidRow(
                                column(4,selectInput("cr_flags",
                                                     "CR options",
                                                     choices=c(
                                                         "None"=0,
                                                         "Cross reactivity"=1,
                                                         "Typed CR"=2),
                                   selected=0)),
                                column(4, selectInput("typing_flags",
                                                      "Type options",
                                                      choices=c(
                                                          "None"=0,
                                                          "Strong typing"=1,
                                                          "Weak typing"=2),
                                                      selected=0),
                                       bsTooltip("typing_flags","No typing means that each exposure has unique parameters. Weak typing means that the order and type matter. Strong typing means that parameters are exposure type specific only","top",options=list(container="body"))),
                                column(4,selectInput("form","Model form",
                                                     c("Competitive"=1,
                                                       "Isolated"=2),
                                                     selected=1))
                            ),
                            hr(),
                            
                            ## ISOLATE EXPOSURE PARAMETERS
                            h4(strong("Exposures")),
                            fluidRow(
                                column(3, textInput("ID","ID","A1A")),
                                column(3, numericInput("ti","Time",0,min=0,max=1000)),
                                column(4, uiOutput("choose_exposure_type")),
                                column(3, numericInput("group","Group",1,min=0,max=5)),
                                column(3,uiOutput("choose_exposure_strain")),
                                column(5,uiOutput("choose_exposure_affects")),
                                column(3,h5(strong("Primed?")),checkboxInput("primed","Primed",FALSE))
                                ),
                            fluidRow(column(12,align="center",actionButton("add_exposure",strong("Add")))),
                            ##
                            
                            hr(),
                            h4(strong("Priming and CR parameters")),
                            fluidRow(
                                column(4, numericInput("beta","(log) Beta",-20,min=-20,max=2)),
                                column(4, numericInput("c","(log) c",4,min=-20,max=2)),
                                column(4,numericInput("y0_mod","y0 mod",1,min=0,max=1))
                            ),
                            fluidRow(
                                column(8,uiOutput("choose_exposure_type_cr")),
                                column(4,numericInput("sigma_value","Sigma value",-3,min=-20,max=2))
                            ),
                            hr(),
                            h4(strong("Antigenic seniority modifiers")),
                            fluidRow(
                                column(3,numericInput("mod1","1",1,min=0,max=1)),
                                column(3,numericInput("mod2","2",1,min=0,max=1)),
                                column(3,numericInput("mod3","3",1,min=0,max=1)),
                                column(3,numericInput("mod4","4",1,min=0,max=1))
                            )
                        ),
                        mainPanel(
                            fluidRow(
                                column(12,
                                       plotOutput("main_plot")
                                       )
                            ),
                            fluidRow(
                                column(8,uiOutput("exposure_table_ids")),
                                column(4, actionButton("remove_exposure","Remove"))
                            ),
                            fluidRow(sliderInput("mu", "Boost", value=8, min=0, max=20,step=0.1)),
                            fluidRow(sliderInput("tp", "Time to peak", value=12, min=0, max=25,step=0.1)),
                            fluidRow(sliderInput("dp", "Proportional drop", value=0.5, min=0, max=1,step=0.001)),
                            fluidRow(sliderInput("ts", "Time to 2nd waning phase", value=10, min=0, max=25,step=0.1)),
                            fluidRow(sliderInput("m", "Long term waning (log)", value=log(0.003), min=-10, max=0,step=0.01))
                            
                        )
                        )
               )
)





