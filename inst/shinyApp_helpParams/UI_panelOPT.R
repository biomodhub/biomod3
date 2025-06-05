##========================================= Options panel ===========================================
tabPanel(title = "Options",
         br(),
         fluidPage(
           column(8, style = "background: #eecab8; border: 1px solid #333c34; border-radius: 8px;",
                  tags$h2("params.OPT", style = "color: #c05217"),
                  p(icon("location-crosshairs") ,
                    strong("params.OPT"), " can be use in", 
                    em("BIOMOD_Wrap"), " and ", em("MS_Modeling"), " of ", strong("biomod3.")),
                  p(icon("book"),
                    "These parameters are pass to",
                    a("BIOMOD_Modeling", href = "https://biomodhub.github.io/biomod2/reference/BIOMOD_Modeling.html"),
                    "and",
                    a("bm_ModelingOptions.", href = "https://biomodhub.github.io/biomod2/reference/bm_ModelingOptions.html"),
                    "Check the documentation for definition."),
                  br()
           )),
         br(),
         fluidRow(
           column(6, 
                  selectInput("OPTstrategy", label = "OPT.strategy", choices = c("default", "bigboss", "tuning", "user.defined"), selected = "default"),
                  conditionalPanel (
                    condition = "input.OPTstrategy  == 'user.defined'",
                    selectInput(inputId = "OPT.user", label = "OPT.user", choices = c("", names_Bmo), selected = ""),
                    selectInput(inputId = "OPT.user.base", label = "OPT.user.base", choices = c("default", "bigboss"), selected = ""),
                    selectInput(inputId = "OPT.user.val", label = "OPT.user.val", choices = c("", names_list), selected = "")),
                  br(),
                  actionButton(inputId = "createOPT", label = "Create params.OPT", icon = icon("gears"))),
           
           column(6,
                  br(),
                  conditionalPanel (
                    condition = "input.OPTstrategy  == 'default'",
                    strong("Only default parameter values of default parameters of the single models functions are retrieved.")),
                  conditionalPanel (
                    condition = "input.OPTstrategy  == 'bigboss'",
                    strong("Parameters pre-defined by biomod2 team and that are available in the dataset OptionsBigboss. ")),
                  conditionalPanel (
                    condition = "input.OPTstrategy  == 'tuning'",
                    strong("Calling the bm_Tuning function to try and optimize some default values.")),
                  conditionalPanel (
                    condition = "input.OPTstrategy  == 'user.defined'",
                    strong("Updates default or bigboss parameters with some parameters values defined by the user."),
                    br(),
                    br(),
                    br(),
                    p("OPT.user must be an BIOMOD.models.options objects", br(),"(obtain with", a("bm_ModelingOptions", href = "https://biomodhub.github.io/biomod2/reference/bm_ModelingOptions.html"), ")."),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    p("OPT.user.val must be an list with the structure", br(),"list(name_model = list(namedataset = list(options)))"))
           )
         )
)