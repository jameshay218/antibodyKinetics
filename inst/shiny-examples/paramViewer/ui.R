library(shiny)
library(shinyBS)
library(rhandsontable)
library(ggplot2)
library(gridExtra)

tmax=100

shinyUI(    
    navbarPage("Antibody Kinetics Model",
               #' Parameter exploration panel
               tabPanel("Exposures",
                        ## ISOLATE EXPOSURE PARAMETERS
                        sidebarPanel(
                             
                            fluidRow(
                                selectInput("typing_flags",
                                            "Type options",
                                            choices=c(
                                                "None"=0,
                                                "Strong typing"=1,
                                                "Weak typing"=2),
                                            selected=0)),
                            fluidRow(
                                selectInput("form","Model form",
                                            c("Competitive"=1,
                                              "Isolated"=2),
                                            selected=1)),
                            fluidRow(numericInput("n_strains","No. strains",5,min=1,max=10)),
                            h4(strong("Exposures")),
                            fluidRow(uiOutput("choose_exposure_id")),
                            fluidRow(uiOutput("choose_exposure_type")),
                            fluidRow(uiOutput("choose_exposure_ti")),
                            fluidRow(uiOutput("choose_exposure_group")),
                            fluidRow(uiOutput("choose_exposure_strain")),
                            fluidRow(uiOutput("choose_exposure_affects")),
                            fluidRow(uiOutput("choose_is_primed")),
                            hr(),
                            fluidRow(actionButton("add_exposure",strong("Add")),
                                     actionButton("remove_exposure",strong("Remove")),
                                     actionButton("clear_exposures",strong("Clear")),
                                     downloadButton("export_exposures",strong("Download"))
                                     ),
                            br(),
                            fluidRow(
                                fileInput("exposure_tab_input",strong("Exposure table input")),
                                actionButton("upload_exposures","Upload")
                            )
                        ),
                        mainPanel(
                            h4(strong("Protocol")),
                            fluidRow(
                                plotOutput("protocol_plot")
                            ),
                            hr(),
                            h4(strong("Exposure table")),
                            fluidRow(
                                rHandsontableOutput("exposure_table")
                            )
                        )
                        ),
               tabPanel("Trajectories",
                        sidebarPanel(
                            h4(strong("Main parameters")),
                            fluidRow(
                                column(4,numericInput("tmax","Max time", 100, min=10,max=1000)),
                                column(4,numericInput("lower_bound","Lower bound",0,min=-1000,max=0)),
                                column(4,numericInput("max_titre","Max Titre",15,min=0,max=25))
                            ),
                            hr(),
                            fluidRow(uiOutput("exposure_table_id")),
                            fluidRow(selectInput("exposure_exposure_select","Exposure strain","none",selected=1)),
                            fluidRow(selectInput("exposure_strain_select","Affected strain","none",selected=1)),
                            fluidRow(column(8,uiOutput("select_mu")),column(2,actionButton("randomise_mu","Random"))),
                            fluidRow(column(8,uiOutput("select_tp")),column(2,actionButton("randomise_tp","Random"))),
                            fluidRow(column(8,uiOutput("select_dp")),column(2,actionButton("randomise_dp","Random"))),
                            fluidRow(column(8,uiOutput("select_ts")),column(2,actionButton("randomise_ts","Random"))),
                            fluidRow(column(8,uiOutput("select_m")),column(2,actionButton("randomise_m","Random"))),
                            fluidRow(actionButton("randomise_all","Randomise all")),
                            hr(),
                            h4(strong("Antigenic seniority modifiers")),
                            fluidRow(
                                column(3,numericInput("mod1","1",1,min=0,max=1)),
                                column(3,numericInput("mod2","2",1,min=0,max=1)),
                                column(3,numericInput("mod3","3",1,min=0,max=1)),
                                column(3,numericInput("mod4","4",1,min=0,max=1))
                            ),
                            fluidRow(
                                downloadButton("download_all",strong("Download Parameters"))
                            )
                        ),
                        mainPanel(
                            fluidRow(
                                plotOutput("main_plot")
                            )                           
                        )
                        ),
               tabPanel("Cross reactivity",
                        sidebarPanel(
                            fluidRow(
                                selectInput("cr_flags",
                                            "CR options",
                                            choices=c(
                                                "None"=0,
                                                "Cross reactivity"=1,
                                                "Typed CR"=2),
                                            selected=0)),
                            h4(strong("Priming and CR parameters")),
                            fluidRow(
                                column(4, numericInput("beta","(log) Beta",-20,min=-20,max=2)),
                                column(4, numericInput("c","(log) c",4,min=-20,max=2)),
                                column(4, numericInput("y0_mod","y0 mod",-20,min=0,max=1))
                            ),
                            fluidRow(
                                column(8,uiOutput("choose_exposure_type_cr")),
                                column(4,numericInput("sigma_value","Sigma value",-3,min=-20,max=2))
                            ),
                            fluidRow(
                                fileInput("antigenic_tab_input",strong("Antigenic distanceinput")),
                                actionButton("upload_antigenic_distances","Upload")
                            ),
                            fluidRow(
                                downloadButton("export_cr",strong("Download"))
                            ),
                            hr()
                        ),
                        mainPanel(
                            h4(strong("Exposure distances")),
                            fluidRow(
                                rHandsontableOutput("antigenic_table")
                            )
                        )
                        )
               
               )
)





