# bm_shinyOptions ---------------------------------------------------------
##' @name bm_shinyOptions
##' @author Helene Blancheteau
##' 
##' @title Shiny application to help set the options for the algorithms
##' in biomod2
##' 
##' @description ...
##'              
##'  
##' @export
##' 
##' 



bm_shinyOptions = function()
{
  if (!isNamespaceLoaded("shiny")) { 
    if (!requireNamespace('shiny', quietly = TRUE)) stop("Package 'shiny' not found")
  }
  
  if (!isNamespaceLoaded("biomod2")) { 
    if (!requireNamespace('biomod2', quietly = TRUE)) stop("Package 'biomod2' not found")
  }
  
  library(shiny)
  library(shinyWidgets)
  
  appDir <- system.file("shinyApp", package = "biomod3")
  #appDir = "~/biomod_stacked/ShinyParameters/inst/shinyApp"
  if (appDir == "") {
    stop("Could not find shinyApp directory. Try re-installing `biomod3`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
