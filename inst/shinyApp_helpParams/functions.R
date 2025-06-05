namesObjectsByType <- function(type){
  names_obj <- objects(envir = globalenv())
  text  <- paste0(" inherits(", names_obj, ", '", type, "' )", sep = "")
  yes <- sapply(text , FUN = function(x){eval(parse(text = x), envir = globalenv())})
  names_type <- names_obj[yes]
  if(rlang::is_empty(names_type)) {names_type <- paste0("No ", type, " object available")}
  names_type
}

# print_params_PA <- function(PA.strategy, PA.nb.rep, PA.nb.absences, PA.dist.min, PA.dist.max, PA.fact.aggr, PA.sre.quant)
# 
# 
# observeEvent(input$createPA, {
#   print(paste0("params.PA <- list('PA.strategy' = '" ,input$PAstrategy, "', 'PA.nb.rep' = ", input$PA.nb.rep, 
#                ", 'PA.nb.absences' =", input$PA.nb.absences, ", 'PA.dist.min' =", input$PA.dist.min,
#                ", 'PA.dist.max' =", input$PA.dist.max, ", 'PA.fact.aggr' = ", input$PA.fact.aggr,
#                ", 'PA.sre.quant' =", input$PA.sre.quant, ")"))
# })
# 
# observeEvent(input$createCV, {
#   print(paste0("params.CV <- list( 'CV.nb.rep' = ", input$CV.nb.rep, ", 'CV.strategy' = '", input$CVstrategy,
#                "', 'CV.perc' = ", input$CV.perc, ", 'CV.k' = ", input$CV.k, ", 'CV.balance' = '", input$CV.balance,
#                "', 'CV.env.var' = '", input$CV.env.var, "', 'CV.strat' = '", input$CV.strat, "' )"))
# })
# 
# observeEvent(input$createOPT, {
#   print(paste0("params.OPT <- list( 'OPT.strategy' = '", input$OPTstrategy, "' , 'OPT.user.val' = ", input$OPT.user.val,
#                ", 'OPT.user.base' = '", input$OPT.user.base, "', 'OPT.user' = ", input$OPT.user, ")"))
# })
# 
# observeEvent(input$createEM, {
#   print(paste0("params.EM <- list( 'models.chosen' =", input$models.chosen, ", 'em.by' = '", input$em.by, 
#                "', 'metric.select' = '", input$metric.select, "', 'metric.select.thresh' =", input$metric.select.thresh,
#                #", 'metric.select.table' = ", input$metric.select.table,
#                ", 'metric.select.dataset' = '", input$metric.select.dataset, 
#                "', 'EMci.alpha' =", input$EMci.alpha, ", 'EMwmean.decay' = ", input$EMwmean.decay, ")"))
# })