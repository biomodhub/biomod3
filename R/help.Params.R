# help.Params ---------------------------------------------------------
##' @name help.Params
##' @author Helene Blancheteau
##' 
##' @title Shiny application to help set the parameter for the parameters in 
##' BIOMOD_Wrap or MS functions.
##' 
##' @description ...
##'              
##'  
##' @export
##' 
##' 

help.Params = function()
{
  if (!isNamespaceLoaded("shiny")) { 
    if (!requireNamespace('shiny', quietly = TRUE)) stop("Package 'shiny' not found")
  }
  
  if (!isNamespaceLoaded("biomod2")) { 
    if (!requireNamespace('biomod2', quietly = TRUE)) stop("Package 'biomod2' not found")
  }
  
  library(shiny)
  library(shinyWidgets)
  
  params.PA <- list()
  params.CV <- list()
  params.OPT <- list()
  params.EM <- list()
  
  appDir <- system.file("shinyApp_helpParams", package = "biomod3")
  if (appDir == "") {
    stop("Could not find shinyApp directory. Try re-installing `biomod3`.", call. = FALSE)
  }
  #shiny::runApp(appDir, display.mode = "normal")
  
  output <- capture.output(shiny::runApp(appDir, display.mode = "normal"))
  output <- unlist(strsplit(output, '"'))
  if(!is.null(output)){
    output <- output[seq(2,length(output), by = 2)]
  }
  eval(parse(text = output))
  
  return(list('params.PA' = params.PA, 
              'params.CV' = params.CV, 
              'params.OPT' = params.OPT, 
              'params.EM' = params.EM))
}
