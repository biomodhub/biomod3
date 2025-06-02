source("functions.R", local = T)
library(shinyjs)
library(shinyWidgets)
library(shinyalert)
library(ggplot2)

# rm(list = ls())

# ## SHINY packages
# list.packages = c("shiny", "shinyFiles", "shinyalert", "shinyDirectoryInput", "shinyjs", "shinyWidgets", "shinyhelper")
# ## DATA MANIPULATION packages
# list.packages = c(list.packages
#                   , "data.table", "foreach", "dplyr", "DT", "reshape2")
# ## GRAPHICS packages
# list.packages = c(list.packages
#                   , "ggplot2", "ggthemes", "ggdendro", "ggrepel", "ggExtra", "gridExtra"
#                   , "viridis", "RColorBrewer", "plotly", "cowplot", "patchwork")
# ## OTHER packages
# list.packages = c(list.packages
#                   , "markdown", "zip", "raster", "rintrojs"
#                   , "SPOT" ## designLHD
#                   , "ade4" ## quasieuclid
#                   , "FD" ## gowdis
# )
# 
# load.shinyDirectoryInput = requireNamespace("shinyDirectoryInput")
# load.devtools = requireNamespace("devtools")
# if (!load.shinyDirectoryInput)
# {
#   if (!load.devtools)
#   {
#     install.packages("devtools")
#   }
#   library(devtools)
#   devtools::install_github('wleepang/shiny-directory-input')
#   # devtools::install_github("thomasp85/shinyFiles")
# }
# check.packages = sapply(list.packages, .loadPackage)
# load.packages = sapply(list.packages, require, character.only = TRUE)

bmformat_potential <- namesObjectsByType(type = 'BIOMOD.formated.data')
caliblines_potential <- c(namesObjectsByType(type = 'matrix'), namesObjectsByType(type = 'data.frame'))


ui <- fluidPage(
  useShinyjs(),
  shinyWidgets::setBackgroundColor("#FAFAFA"),
  titlePanel(
    fluidRow(
      style = HTML("background-color: #9BCD9B ; margin-top: 20px; margin-bottom: 20px;")
      , column(10, headerPanel("ShinyOptions", windowTitle = "ShinyOptions"))
      , column(2, tags$img(#src = "https://github.com/biomodhub/biomod2/blob/master/docs/articles/pictures/LogoBiomod.png"
                            src = "biomod2.png"
                           , align = "right", height = '70px')
      )
    )
  ), 
  
  br(),
  fluidRow(
  column(3,
         selectInput("datatype", label = "Select datatype",
                     choices = c("binary", "abundance", "count", "relative", "ordinal"), selected = "binary",
                     width = "200px"),
         selectInput("base", label = "Select options base",
                     choices = c("bigboss", "default"), selected = "bigboss",
                     width = "200px"),
         selectInput("models", "Select models:", "", multiple = TRUE, width = "200px"),
         br()
         ),
  column(4,
         textInputIcon(inputId = "objectName", label = "Choose a name for your options object (create in your environnement)", value = "myBiomodOptions",
                       icon = list(NULL, icon("circle-pause")), width = "250px"),
         actionButton(inputId = "alldatasets", label = "Create options for all datasets",
                      icon = icon("screwdriver-wrench"), width = "250px")
         ),
  column(5,
         p("If you wish to customise each dataset: "),
         actionButton(inputId = "bmformat_button", label = "Select your biomod2 data",
                      icon = icon("rectangle-list"), width = "200px"),
         br(),
         br(),
         hidden(actionButton(inputId = "caliblines_button", label = "Select your CV table",
                             icon = icon("rectangle-list"), width = "200px")),
         br(),
         br(),
         hidden(strong(id = "question", "For which dataset you want to save these options ?")),
         fluidRow(
           column(6,
                  hidden(selectInput(inputId = "whichPA", label = "PA", choices = c("All the same"), selected = "All the same", width = "150px"))
           ),
           column(6,
                  hidden(selectInput(inputId = "whichCV", label = "CV", choices = c("All the same"), selected = "All the same", width = "150px"))
           )
         ),
         hidden(
           actionButton(inputId = "datasetSave", label = "Validate for this dataset", icon = icon("floppy-disk"), width = "200px"),
           actionButton(inputId = "progress", label = "", icon = icon("hourglass-half"))
         ),

         br(),
         br(),
         hidden(
           actionButton(inputId = "finalOptions", label = "create final options", icon = icon("screwdriver-wrench"), width = "200px")
         )
  )
  ),
  
  br(),
  tags$style(HTML("
                  .tabbable > .nav > li > a                  {font-weight: bold; background-color: #EDEDED;  color:#9BCD9B}
                  .tabbable > .nav > li[class=active]    > a {background-color: #FAFAFA; color:#548B54}
                  ")),
  
  tabsetPanel(id = "inTabset",
              tabPanel("ANN", uiOutput("panelANN")),
              tabPanel("CTA", uiOutput("panelCTA")),
              tabPanel("FDA", uiOutput("panelFDA")),
              tabPanel("GAM", uiOutput("panelGAM")),
              tabPanel("GBM", uiOutput("panelGBM")),
              tabPanel("GLM", uiOutput("panelGLM")),
              tabPanel("MARS", uiOutput("panelMARS")),
              tabPanel("MAXENT", uiOutput("panelMAXENT")),
              tabPanel("MAXNET", uiOutput("panelMAXNET")),
              tabPanel("RF", uiOutput("panelRF")),
              tabPanel("RFd", uiOutput("panelRFd")),
              tabPanel("SRE", uiOutput("panelSRE")),
              tabPanel("XGBOOST", uiOutput("panelXGBOOST"))
  )

)

server <- function(input, output, session){
  
  ### Models tabs in function of data.type and options base =====================
  hide_several_tabs("inTabset", c("ANN","CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXENT", "MAXNET", "RF",
                                  "RFd", "SRE", "XGBOOST"))
  models_possible <- reactive({
    if (input$datatype == "binary"){
      models_possible <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'RFd', 'SRE', 'XGBOOST')
    } else if (input$datatype == "ordinal"){
      models_possible <- c('CTA', 'FDA', 'GAM', 'GLM', 'MARS', 'RF', 'XGBOOST')
    } else {
      models_possible <- c('CTA', 'GAM', 'GBM', 'GLM', 'MARS', 'RF', 'XGBOOST')
    }
    models_possible
  })
  
  options_base <- reactive({
    ddpcr::quiet(options <- biomod2::bm_ModelingOptions(data.type = input$datatype,
                                           strategy = input$base,
                                           models = models_possible()),
                   )
    options
  })
  
  
  observe({
    updateSelectInput(session, "models", choices = models_possible()
    )})
  
  observeEvent(ignoreInit = TRUE, list(input$datatype, input$models), {
    hide_several_tabs("inTabset", setdiff(c("ANN","CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXENT", "MAXNET", "RF",
                                    "RFd", "SRE", "XGBOOST"), input$models))
    for (i in 1:length(input$models)){
      showTab(inputId = "inTabset", target = input$models[i])
    }
  })
  
  
  observeEvent(input$base, {
               output$panelANN <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "ANN")))})
               output$panelCTA <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "CTA")))})
               output$panelFDA <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "FDA")))})
               output$panelGAM <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "GAM")))})
               output$panelGBM <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "GBM")))})
               output$panelGLM <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "GLM")))})
               output$panelMARS <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "MARS")))})
               output$panelMAXNET <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "MAXNET")))})
               output$panelMAXENT <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "MAXENT")))})
               output$panelRF <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "RF")))})
               output$panelRFd <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "RFd")))})
               output$panelSRE <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "SRE")))})
               output$panelXGBOOST <- renderUI({eval(parse(text = ui_tabOptions(options_base(), "XGBOOST")))})
               }
               )
  
  #=============================================================================
  
  
  ### Basic: for all datasets ==================================================
  observeEvent(input$alldatasets, {
    eval(parse(text = get_inputs(options_base(), input$models)))
    options_created <- create_options(list_inputs, options_base(), "for_all_datasets",
                                      input$models, input$datatype, input$base, input$objectName)
    assign(input$objectName, options_created, envir = .GlobalEnv)
    updateTextInputIcon(session, "objectName", icon  = list(NULL, icon("circle-check")))
  })
  
  observeEvent(input$objectName, {
    updateTextInputIcon(session, "objectName", icon = list(NULL, icon("circle-pause")))
  })
  #=============================================================================
  
  
  ### Get bm.format and cvtable ================================================
  observeEvent(input$bmformat_button,
    shinyalert(html = TRUE, text = tagList(
      selectInput("bmformat_name_in_env", "Your BIOMOD.formated.data object", choices = bmformat_potential)
    ))
  )
  
  # bm.format <- NULL  
  # cvtable <- NULL
  RV = reactiveValues(progress_table = data.frame(),
                      bm.format = NULL, 
                      calib.lines = NULL, 
                      PA_names = c(),
                      CV_names = c(),
                      user.val = list())
  
  
  observeEvent(input$bmformat_name_in_env, {
    if (!is.null(input$bmformat_name_in_env) && input$bmformat_name_in_env != "No BIOMOD.formated.data object available"){
      RV$bm.format <- get(input$bmformat_name_in_env, envir = .GlobalEnv)
      RV$PA_names <- get_PA_datasets(RV$bm.format)
    }
    if(!is.null(RV$bm.format)){
      shinyjs::show("caliblines_button", anim = T)
      shinyjs::show("question", anim = T)
      shinyjs::show("whichPA", anim = T)
      shinyjs::show("whichCV", anim = T)
      updateSelectInput(session, "whichPA", choices = RV$PA_names)
      shinyjs::show("datasetSave", anim = T)
      shinyjs::show("progress", anim = T)
      #shinyjs::show("finalOptions", anim = T)
      RV$progress_table <- init_progress(RV$PA_names, c("All the same"))
    }
  })
  
  observeEvent(input$caliblines_button,
               shinyalert(html = TRUE, text = tagList(
                 selectInput("caliblines_name_in_env", "Your CV table object", choices = caliblines_potential)
               ))
  )
  
  observeEvent(input$caliblines_name_in_env, {
    if (!is.null(input$caliblines_name_in_env) && input$caliblines_name_in_env != "No matrix object available" && input$caliblines_name_in_env != "No data.frame object available"){
      RV$calib.lines <- get(input$caliblines_name_in_env, envir = .GlobalEnv)
    }
    if(!is.null(RV$calib.lines)){
      updateSelectInput(session, "whichCV", choices = get_CV_datasets(RV$bm.format, RV$calib.lines))
      RV$CV_names <- get_CV_datasets(RV$bm.format, RV$calib.lines)
    }
    RV$progress_table <- init_progress(RV$PA_names, get_CV_datasets(RV$bm.format, RV$calib.lines))
  })
  #=============================================================================
  

  ### Progress table ===========================================================
  
  observeEvent(input$progress, {
    shinyalert(
      title = "Datasets you have already saved", html = T,
      text = tagList(plotOutput("progressPlot")),
      closeOnClickOutside = TRUE, showConfirmButton = FALSE,
    )
  } 
  )
  
  output$progressPlot <- renderPlot({progress_plot(RV$progress_table)})
  #=============================================================================
  
  
  ### Validate Options for one dataset =========================================
  
  observeEvent(input$datasetSave, {
    ## Update user.val
    eval(parse(text = get_inputs(options_base(), input$models)))
    RV$user.val <- update_userval(RV$user.val, input$whichPA, input$whichCV,
                                  list_inputs, RV$PA_names, RV$CV_names, input$models, options_base())
    ## Updata progress_table 
    RV$progress_table <- update_progress(RV$progress_table, input$whichPA, input$whichCV)
    
    if(all(RV$progress_table$progress)){
      shinyjs::show("finalOptions", anim = T)
    }
  })
  #=============================================================================
  
  
  ### Create final options =====================================================
  observeEvent(input$finalOptions, {
    options_created <- biomod2::bm_ModelingOptions(data.type = input$datatype,
                                                   bm.format = RV$bm.format,
                                                   calib.lines = RV$calib.lines,
                                                   strategy = "user.defined",
                                                   models = input$models,
                                                   user.val = RV$user.val,
                                                   user.base = input$base
                                                   )
    assign(input$objectName, options_created, envir = .GlobalEnv)
    updateTextInputIcon(session, "objectName", icon  = list(NULL, icon("circle-check")))
  })
  
  #=============================================================================
  ## Promis Maya je vais nettoyer un peu tout ça! 
}

###################################################################################################################################
# Create a Shiny app object
shinyApp(ui = ui, server = server)

