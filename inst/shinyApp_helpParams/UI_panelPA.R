##=========================================PA panel ===========================================
tabPanel(title = "Pseudo Absences",
         br(),
         fluidPage(
           column(8, style = "background: #eecab8; border: 1px solid #333c34; border-radius: 8px;",
                  tags$h2("params.PA", style = "color: #c05217"),
                  p(icon("location-crosshairs") ,
                    strong("params.PA"), " can be use in", 
                    em("BIOMOD_Wrap"), " and ", em("MS_FormatingData"), " of ", strong("biomod3.")),
                  p(icon("book"),
                    "These parameters are pass to",
                    a("BIOMOD_FormatingData", href = "https://biomodhub.github.io/biomod2/reference/BIOMOD_FormatingData.html"),
                    "and",
                    a("bm_PseudoAbsences.", href = "https://biomodhub.github.io/biomod2/reference/bm_PseudoAbsences.html"),
                    "Check the documentation for definition."),
                  br()
           )),
         br(),
         fluidRow(
           #column(1),
           column(5, style = "background: #fdf2ed; border-radius: 8px;",
                  selectInput("PAstrategy", label = "PA.strategy", choices = c("random", "disk", "sre"), selected = "random"),
                  numericInput(inputId = "PA.nb.rep", label = "PA.nb.rep", value = 3, min = 0),
                  numericInput(inputId = "PA.nb.absences", label = "PA.nb.absences", value = 1000, min = 0),
                  
                  conditionalPanel (
                    condition = "input.PAstrategy == 'disk'",
                    numericInput(inputId = "PA.dist.min", label = "PA.dist.min", value = NULL, min = 0),
                    numericInput(inputId = "PA.dist.max", label = "PA.dist.max", value = NULL, min = 0),
                    numericInput(inputId = "PA.fact.aggr", label = "PA.fact.aggr", value = NULL, min = 0)),
                  
                  conditionalPanel (
                    condition = "input.PAstrategy == 'sre'",
                    numericInput(inputId = "PA.sre.quant", label = "PA.sre.quant", value = 0.025, min = 0)),
                  br(),
                  actionButton(inputId = "createPA", label = "Create params.PA", icon = icon("gears"))),
           column(1),
           column(5, align = "center",
                  br(),
                  br(),
                  conditionalPanel (
                    condition = "input.PAstrategy == 'random'",
                    tags$img(src = "PA_random.png", height = "250px"),
                    em("Random selection of pseudo absences")),
                  conditionalPanel (
                    condition = "input.PAstrategy == 'disk'",
                    tags$img(src = "PA_disk.png", height = "250px"),
                    em("Selection of pseudo absences with a minimal distance disk")),
                  conditionalPanel (
                    condition = "input.PAstrategy == 'sre'",
                    tags$img(src = "PA_sre.png", height = "250px"),
                    em("Selection of pseudo absences taking in account the environnemental conditions"))
           )
         )
)