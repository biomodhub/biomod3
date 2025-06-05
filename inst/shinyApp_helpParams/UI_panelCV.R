##=========================================CV panel ===========================================
tabPanel(title = "CrossValidation",
         br(),
         fluidPage(
           column(8, style = "background: #eecab8; border: 1px solid #333c34; border-radius: 8px;",
                  tags$h2("params.CV", style = "color: #c05217"),
                   p(icon("location-crosshairs"),
                     strong("params.CV"), " can be use in",
                     em("BIOMOD_Wrap"), " and ", em("MS_Modelling"), " of ", strong("biomod3.")),
                   p(icon("book"),
                     "These parameters are pass to",
                     a("BIOMOD_Modeling", href = "https://biomodhub.github.io/biomod2/reference/BIOMOD_Modeling.html"),
                    "and", 
                    a("bm_CrossValidation.", href = "https://biomodhub.github.io/biomod2/reference/bm_CrossValidation.html"), "Check the documentation for definition."),
                  br()
           )),
         br(),
         fluidRow(
           column(6, 
                  selectInput("CVstrategy", label = "CV.strategy", choices = c("random", "kfold", "block", "strat", "env"), selected = "random"),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'random' | input.CVstrategy == 'kfold'",
                    numericInput(inputId = "CV.nb.rep", label = "Cv.nb.rep", value = 2, min = 0)),
                  conditionalPanel (
                    condition = "input.CVstrategy  == 'random'",
                    numericInput(inputId = "CV.perc", label = "CV.perc", value = 0.7, min = 0, max = 1)),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'env' | input.CVstrategy == 'kfold' | input.CVstrategy == 'strat'",
                    numericInput(inputId = "CV.k", label = "CV.k", value = NULL, min = 2)),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'strat' | input.CVstrategy == 'env'",
                    selectInput(inputId = "CV.balance", label = "CV.balance", choices = c("presences", "absences"), selected = "presences")),
                  conditionalPanel (
                    condition = "input.CVstrategy  == 'strat'",
                    selectInput(inputId = "CV.strat", label = "CV.strat", choices = c("x","y", "both"), selected = "both")),
                  br(),
                  actionButton(inputId = "createCV", label = "Create params.CV", icon = icon("gears"))),
           column(6,
                  br(),
                  br(),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'random'",
                    tags$img(src = "CV_random.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'kfold'",
                    tags$img(src = "CV_kfold.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'env'",
                    tags$img(src = "CV_env.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'block'",
                    tags$img(src = "CV_block.png", align = "center", height = "250px")),
                  conditionalPanel (
                    condition = "input.CVstrategy == 'strat'",
                    tags$img(src = "CV_strat.png", align = "center", height = "250px"))
           )
         )
         
)