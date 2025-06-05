source("functions.R", local = T)
library(shiny)
library(shinyjs)

all <- "all"
metrics <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
             , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC'
             , 'BOYCE', 'MPA', 'Accuracy', 'Recall', 'Precision', 'F1',
             'RMSE','MSE','MAE','Rsquared','Rsquared_aj','Max_error')

names_list <- namesObjectsByType('list')
names_df <- namesObjectsByType('data.frame')
names_char <- namesObjectsByType('character')
names_Bmo <- namesObjectsByType('BIOMOD.models.options')


#=======================================================================

ui <- fluidPage(
  shinyWidgets::setBackgroundColor("#fdf2ed"),
  titlePanel(
    fluidRow(
      style = HTML("background-color: #eecab8 ; margin-top: 5px; margin-bottom: 5px; color: #c05217")
      , column(10, headerPanel("help.Params", windowTitle = "help.Params"))
      , column(2, tags$img(#src = "https://github.com/biomodhub/biomod2/blob/master/docs/articles/pictures/LogoBiomod.png"
        src = "logo.png"
        , align = "right", height = '70px')
      )
    )
  ), 
  br(),
  
  #tags$img(src="aa.png", align = "right"),
  tags$style(HTML("
                  .tabbable > .nav > li > a                  {font-weight: bold; background-color: #fdf2ed ;  color:#c05217 ; border: 1px solid #ADABA9}
                  .tabbable > .nav > li[class=active]    > a {background-color: #eecab8; color:#c05217}
                  ")),
  #tags$figure(tags$img(src="logo.png", align = "right", height = '70px')),
  tabsetPanel(
    
    source("UI_panelPA.R", local = T)$value,
    source("UI_panelCV.R", local = T)$value,
    source("UI_panelOPT.R", local = T)$value,
    source("UI_panelEM.R", local = T)$value
    
  )
)

#=======================================================================

server <- function(input, output){
  #bs_themer()
  observeEvent(input$createPA, {
    print(paste0("params.PA <- list('PA.strategy' = '" ,input$PAstrategy, "', 'PA.nb.rep' = ", input$PA.nb.rep, 
                 ", 'PA.nb.absences' =", input$PA.nb.absences, ", 'PA.dist.min' =", input$PA.dist.min,
                 ", 'PA.dist.max' =", input$PA.dist.max, ", 'PA.fact.aggr' = ", input$PA.fact.aggr,
                 ", 'PA.sre.quant' =", input$PA.sre.quant, ")"))
  })
  
  observeEvent(input$createCV, {
    print(paste0("params.CV <- list( 'CV.nb.rep' = ", input$CV.nb.rep, ", 'CV.strategy' = '", input$CVstrategy,
                 "', 'CV.perc' = ", input$CV.perc, ", 'CV.k' = ", input$CV.k, ", 'CV.balance' = '", input$CV.balance,
                 "', 'CV.env.var' = '", input$CV.env.var, "', 'CV.strat' = '", input$CV.strat, "' )"))
  })
  
  observeEvent(input$createOPT, {
    print(paste0("params.OPT <- list( 'OPT.strategy' = '", input$OPTstrategy, "' , 'OPT.user.val' = ", input$OPT.user.val,
                 ", 'OPT.user.base' = '", input$OPT.user.base, "', 'OPT.user' = ", input$OPT.user, ")"))
  })
  
  observeEvent(input$createEM, {
    print(paste0("params.EM <- list( 'models.chosen' =", input$models.chosen, ", 'em.by' = '", input$em.by, 
                 "', 'metric.select' = '", input$metric.select, "', 'metric.select.thresh' =", input$metric.select.thresh,
                 #", 'metric.select.table' = ", input$metric.select.table,
                 ", 'metric.select.dataset' = '", input$metric.select.dataset, 
                 "', 'EMci.alpha' =", input$EMci.alpha, ", 'EMwmean.decay' = ", input$EMwmean.decay, ")"))
  })
}


###################################################################################################################################
# Create a Shiny app object
shinyApp(ui = ui, server = server)