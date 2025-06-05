##========================================= EM panel ===========================================
tabPanel(title = "Ensemble Models",
         br(),
         fluidPage(
           column(9, style = "background: #eecab8; border: 1px solid #333c34; border-radius: 8px;",
                  tags$h2("params.EM", style = "color: #c05217"),
                  p(icon("location-crosshairs") ,
                    strong("params.EM"), " can be use in", 
                    em("BIOMOD_Wrap"), " and ", em("MS_EnsembleModeling"), " of ", strong("biomod3.")),
                  p(icon("book"),
                    "These parameters are pass to",
                    a("BIOMOD_EnsembleModeling", href = "https://biomodhub.github.io/biomod2/reference/BIOMOD_EnsembleModeling.html"),
                    "Check the documentation for definition."),
                  br()
           )),
         br(),
         fluidRow(
           column(6,
                  selectInput(inputId = "models.chosen", label = "models.chosen", choices = c(all, names_char), selected = "all"),
                  selectInput(inputId = "emby", label = "em.by", choices = c('all', 'algo', 'PA', 'PA+algo', 'PA+run'), selected = "all"),
                  selectInput(inputId = "metric.select", label = "metric.select", choices = metrics, selected = NULL, multiple = T),
                  numericInput(inputId = "metric.select.thresh", label = "metric.select.thresh", value = 0.7, min = 0, max = 1),
                  #selectInput(inputId = "metric.select.table", label = "metric.select.table", choices = c("", names_df), selected = ""),
                  selectInput(inputId = "metric.select.dataset", label = "metric.select.dataset", choices = c('calibration', 'validation'), selected = "validation"),
                  numericInput(inputId = "EMci.alpha", label = em("EMci.alpha"), value = 0.05, min = 0, max = 1),
                  textInput("EMwmean.decay", label = em("EMwmean.decay"), value = "'proportional'"),
                  br(),
                  actionButton(inputId = "createEM", label = "Create params.EM", icon = icon("gears"))),
           column(6,
                  br(),
                  br(),
                  conditionalPanel (
                    condition = "input.emby == 'all'",
                    tags$img(src = "EMbyall.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.emby == 'algo'",
                    tags$img(src = "EMbyalgo.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.emby == 'PA'",
                    tags$img(src = "EMbyPA.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.emby == 'PA+algo'",
                    tags$img(src = "EMbyPAalgo.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.emby == 'PA+run'",
                    tags$img(src = "EMbyPArun.png", align = "center", height = "250px"))
           )
         )
)