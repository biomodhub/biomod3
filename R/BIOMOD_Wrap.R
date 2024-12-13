# BIOMOD_Wrap ---------------------------------------------------------
##' @name BIOMOD_Wrap
##' @author Helene Blancheteau
##' 
##' @title Run a range of species distribution models
##' 
##' @description This function allows to run a workflow for one or several species using biomod2. 
##' The parameters could be change inside the arguments 'params. 's. 
##' 
##' 
##' @param ms.project.name a \code{character} corresponding to the name of your project for your multispecies modelling
##' @param dir.name (\emph{optional, default} \code{.}) \cr
##' A \code{character} corresponding to the modeling folder
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' 
##' @param resp.name a \code{character vector} corresponding to the species name
##' 
##' @param resp.var a \code{vector} or a \code{\link[terra:vect]{SpatVector}} containing your response variable (See Details).
##' @param data.type a \code{character}, corresponding to the response data type to be used, must be either 
##' \code{binary}, \code{count}, \code{ordinal}, \code{relative}, or \code{abundance}. If data.type is not provided,
##' \code{biomod2} will try to guess.
##' 
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' 
##' @param resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be
##' used to build the species distribution model(s)
##' 
##' @param eval.resp.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector}, a \code{\link[terra:vect]{SpatVector}}
##' without associated data (\emph{if presence-only}), 
##' or a \code{\link[terra:vect]{SpatVector}} object containing binary data  
##' (\code{0} : absence, \code{1} : presence, \code{NA} : indeterminate) 
##'  for a single species that will be used to 
##' evaluate the species distribution model(s) with independent data
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##'  \code{SpatialPoints}  (if presence-only) or \code{SpatialPointsDataFrame}
##'  object containing binary data.}
##' @param eval.expl.var.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables (in 
##' columns or layers) that will be used to evaluate the species distribution model(s) with 
##' independent data.
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' @param eval.resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to evaluate 
##' the species distribution model(s) with independent data
##' 
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' @param models.pa (\emph{optional, default} \code{NULL}) \cr 
##' A \code{list} containing for each model a \code{vector} defining which pseudo-absence datasets 
##' are to be used, must be among \code{colnames(bm.format@PA.table)}
##' @param weights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to observation weights (one per 
##' observation, see Details)
##' @param prevalence (\emph{optional, default} \code{NULL}) \cr 
##' A \code{numeric} between \code{0} and \code{1} corresponding to the species prevalence to 
##' build '\emph{weighted response weights}' (see Details)
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{OR}, 
##' \code{ORSS}, \code{BOYCE}, \code{MPA}, \code{RMSE}, \code{MAE}, \code{MSE}, \code{Rsquared}, \code{Rsquared_aj},
##' \code{Max_error}, \code{Accuracy}, \code{Recall}, \code{Precision}, \code{F1}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param scale.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions should be scaled with a 
##' binomial GLM or not
##' 
##' @param em.algo a \code{vector} corresponding to the ensemble models that will be computed, 
##' must be among \code{'EMmean'}, \code{'EMmedian'}, \code{'EMcv'}, \code{'EMci'}, 
##' \code{'EMca'}, \code{'EMwmean'}
##' 
##' @param params.PA a \code{list} with the species names associated to the parameters of PA
##' @param params.CV a \code{list} with the species names associated to the parameters of Cross-Validation. See BIOMOD_Modeling
##' @param params.OPT a \code{list} with the species names associated to the options of the algorithms. See BIOMOD_Modeling
##' 
##' @param filter.raster (\emph{optional, default} \code{FALSE}) \cr 
##' If \code{expl.var} is of raster type, a \code{logical} value defining whether \code{resp.var} 
##' is to be filtered when several points occur in the same raster cell
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' 
##' 
##' @return
##' 
##' A \code{\link{MS.models.out}} object acting as a proxi for the created \code{BIOMOD.models.out}.
##' 
##' @importFrom foreach foreach %do%
##' @importFrom biomod2 BIOMOD_Modeling
##' 
##' @export
##' 
##' 
###################################################################################################

setGeneric("BIOMOD_Wrap", function(ms.project.name,
                                   dir.name = ".",
                                   modeling.id = as.character(format(Sys.time(), "%s")),
                                   data.type = "binary",
                                   resp.name,
                                   resp.var,
                                   resp.xy = NULL,
                                   expl.var,
                                   eval.resp.var = NULL,
                                   eval.resp.xy = NULL,
                                   eval.expl.var = NULL,
                                   filter.raster = FALSE,
                                   params.PA,
                                   models, ### mettre un défaut ? 
                                   models.pa = NULL,
                                   metric.eval = c("KAPPA", "TSS", "ROC"),
                                   weights = NULL,
                                   prevalence = NULL,
                                   scale.models = FALSE,
                                   var.import = 0,
                                   params.CV,
                                   params.OPT,
                                   em.algo,
                                   params.EM,
                                   seed.val = NULL,
                                   nb.cpu = 1){
  standardGeneric("BIOMOD_Wrap")
})

setMethod(f= "BIOMOD_Wrap", signature(ms.project.name = "missing"), function(ms.project.name,
                                                                             dir.name = ".",
                                                                             modeling.id = as.character(format(Sys.time(), "%s")),
                                                                             data.type = "binary",
                                                                             resp.name,
                                                                             resp.var,
                                                                             resp.xy = NULL,
                                                                             expl.var,
                                                                             eval.resp.var = NULL,
                                                                             eval.resp.xy = NULL,
                                                                             eval.expl.var = NULL,
                                                                             filter.raster = FALSE,
                                                                             params.PA,
                                                                             models, ### mettre un défaut ? 
                                                                             models.pa = NULL,
                                                                             metric.eval = c("KAPPA", "TSS", "ROC"),
                                                                             weights = NULL,
                                                                             prevalence = NULL,
                                                                             scale.models = FALSE,
                                                                             var.import = 0,
                                                                             params.CV,
                                                                             params.OPT,
                                                                             em.algo,
                                                                             params.EM,
                                                                             seed.val = NULL,
                                                                             nb.cpu = 1){
  
  .bm_cat(paste0("Workflow for ", resp.name))

  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_Wrap.check.args(dir.name = dir.name,
                                  modeling.id = modeling.id,
                                  data.type = data.type,
                                  resp.name = resp.name,
                                  resp.var = resp.var,
                                  resp.xy = resp.xy,
                                  expl.var = expl.var,
                                  eval.resp.var = eval.resp.var,
                                  eval.resp.xy = eval.resp.xy,
                                  eval.expl.var = eval.expl.var,
                                  filter.raster = filter.raster,
                                  params.PA = params.PA,
                                  models = models,
                                  models.pa = models.pa,
                                  metric.eval = metric.eval,
                                  weights = weights,
                                  prevalence = prevalence,
                                  var.import = var.import,
                                  params.CV = params.CV,
                                  params.OPT = params.OPT,
                                  em.algo = em.algo,
                                  params.EM = params.EM,
                                  seed.val = seed.val,
                                  nb.cpu = nb.cpu
  )
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  cat("\n\t > Formating Data")
  output <- capture.output(formated.data <- BIOMOD_FormatingData(dir.name = dir.name,
                                                                 resp.var = resp.var,
                                                                 expl.var = expl.var,
                                                                 resp.xy = resp.xy,
                                                                 resp.name = resp.name,
                                                                 eval.resp.var = eval.resp.var,
                                                                 eval.expl.var = eval.expl.var,
                                                                 eval.resp.xy = eval.resp.xy,
                                                                 PA.nb.rep = params.PA$PA.nb.rep,
                                                                 PA.nb.absences = params.PA$PA.nb.absences,
                                                                 PA.strategy = params.PA$PA.strategy,
                                                                 PA.dist.min = params.PA$PA.dist.min,
                                                                 PA.dist.max = params.PA$PA.dist.max,
                                                                 PA.sre.quant = params.PA$PA.sre.quant,
                                                                 PA.fact.aggr = params.PA$PA.fact.aggr,
                                                                 PA.user.table = params.PA$PA.user.table,
                                                                 na.rm = T,
                                                                 seed.val = NULL,
                                                                 filter.raster = filter.raster))
  
  cat("\n\t > Single Models")
  output <- capture.output(single.models <- BIOMOD_Modeling(formated.data,
                                                         modeling.id = modeling.id,
                                                         models = models,
                                                         CV.strategy = params.CV$CV.strategy,
                                                         CV.nb.rep = params.CV$CV.nb.rep,
                                                         CV.perc = params.CV$CV.perc,
                                                         CV.k = params.CV$CV.k,
                                                         CV.balance = params.CV$CV.balance,
                                                         CV.env.var = params.CV$CV.env.var,
                                                         CV.strat = params.CV$CV.strat,
                                                         CV.user.table = params.CV$CV.user.table,
                                                         CV.do.full.models = TRUE,
                                                         OPT.strategy = params.OPT$OPT.strategy,
                                                         OPT.user.val = params.OPT$OPT.user.val,
                                                         OPT.user.base = params.OPT$OPT.user.base,
                                                         OPT.user = params.OPT$OPT.user,
                                                         weights = weights,
                                                         prevalence = prevalence,
                                                         metric.eval = metric.eval,
                                                         var.import = var.import,
                                                         scale.models = scale.models,
                                                         nb.cpu = nb.cpu,
                                                         seed.val = NULL))
  
  cat("\n\t > Ensemble Models")
  output <- capture.output(em.models <- BIOMOD_EnsembleModeling(single.models,
                                                                models.chosen = params.EM$models.chosen,
                                                                em.by = params.EM$em.by,
                                                                em.algo = em.algo,
                                                                metric.select = params.EM$metric.select,
                                                                metric.select.thresh = params.EM$metric.select.thresh,
                                                                metric.select.table = params.EM$metric.select.table,
                                                                metric.select.dataset = params.EM$metric.select.dataset,
                                                                metric.eval = metric.eval,
                                                                var.import = var.import,
                                                                EMci.alpha = params.EM$EMci.alpha,
                                                                EMwmean.decay = params.EM$EMwmean.decay,
                                                                nb.cpu = nb.cpu,
                                                                seed.val = NULL,
                                                                do.progress = FALSE))
  
  ## laisser les messages ?!?! Rajouter quelques messages ? 
  
  .bm_cat("Done")
  return(list("formated.data" = formated.data,
              "single.models" = single.models,
              "ensemble.models" = em.models))
} )

#======================================================================================================

## BIOMOD_Wrap for ms 

setMethod("BIOMOD_Wrap", signature(ms.project.name = "character"), function(ms.project.name,
                                                                            dir.name = ".",
                                                                            modeling.id = as.character(format(Sys.time(), "%s")),
                                                                            data.type = "binary",
                                                                            resp.name,
                                                                            resp.var,
                                                                            resp.xy = NULL,
                                                                            expl.var,
                                                                            eval.resp.var = NULL,
                                                                            eval.resp.xy = NULL,
                                                                            eval.expl.var = NULL,
                                                                            filter.raster = FALSE,
                                                                            params.PA,
                                                                            models, 
                                                                            models.pa = NULL,
                                                                            metric.eval = c("KAPPA", "TSS", "ROC"),
                                                                            weights = NULL,
                                                                            prevalence = NULL,
                                                                            scale.models = FALSE,
                                                                            var.import = 0,
                                                                            params.CV,
                                                                            params.OPT,
                                                                            em.algo,
                                                                            params.EM,
                                                                            seed.val = NULL,
                                                                            nb.cpu = 1){
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_Wrap.check.args(dir.name = dir.name,
                                  modeling.id = modeling.id,
                                  data.type = data.type,
                                  resp.name = resp.name,
                                  resp.var = resp.var,
                                  resp.xy = resp.xy,
                                  expl.var = expl.var,
                                  eval.resp.var = eval.resp.var,
                                  eval.resp.xy = eval.resp.xy,
                                  eval.expl.var = eval.expl.var,
                                  filter.raster = filter.raster,
                                  params.PA = params.PA,
                                  models = models,
                                  models.pa = models.pa,
                                  metric.eval = metric.eval,
                                  weights = weights,
                                  prevalence = prevalence,
                                  var.import = var.import,
                                  params.CV = params.CV,
                                  params.OPT = params.OPT,
                                  em.algo = em.algo,
                                  params.EM = params.EM,
                                  seed.val = seed.val,
                                  nb.cpu = nb.cpu
  )
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  formated.data <- MS_FormatingData(ms.project.name = ms.project.name,
                                     dir.name = dir.name,
                                     resp.name = resp.name,
                                     resp.var = resp.var,
                                     expl.var = expl.var,
                                     data.type = data.type,
                                     resp.xy = resp.xy,
                                     eval.resp.var = eval.resp.var,
                                     eval.expl.var = eval.expl.var,
                                     eval.resp.xy = eval.resp.xy,
                                     params = params.PA, 
                                     #single.formated.data = NULL,
                                     #ms.formated.data = NULL,
                                     filter.raster = filter.raster,
                                     seed.val = NULL)

  
  single.models <- MS_Modeling(formated.data,
                               modeling.id = modeling.id,
                               models = models,
                               params.CV = params.CV,
                               params.OPT = params.OPT,
                               weights = weights,
                               prevalence = prevalence,
                               metric.eval = metric.eval,
                               var.import = var.import,
                               scale.models = FALSE,
                               nb.cpu = nb.cpu,
                               seed.val = NULL)
  
  ### petit gros problème 
  params.EM <- .inversion_species_params(params.EM)
  em.models <- MS_EnsembleModeling(single.models,
                                   models.chosen = params.EM$models.chosen,
                                   em.by = params.EM$em.by,
                                   em.algo = em.algo,
                                   metric.select = params.EM$metric.select,
                                   metric.select.thresh = params.EM$metric.select.thresh,
                                   metric.select.table = params.EM$metric.select.table,
                                   metric.select.dataset = params.EM$metric.select.dataset,
                                   metric.eval = metric.eval,
                                   var.import = var.import,
                                   EMci.alpha = params.EM$EMci.alpha,
                                   EMwmean.decay = params.EM$EMwmean.decay,
                                   nb.cpu = nb.cpu,
                                   seed.val = NULL)
  
  return(list("formated.data" = formated.data,
              "single.models" = single.models,
              "ensemble.models" = em.models))
})


# ---------------------------------------------------------------------------- #

.BIOMOD_Wrap.check.args <- function(dir.name,
                                    modeling.id,
                                    data.type,
                                    resp.name,
                                    resp.var,
                                    resp.xy,
                                    expl.var,
                                    eval.resp.var,
                                    eval.resp.xy,
                                    eval.expl.var,
                                    filter.raster,
                                    params.PA,
                                    models,
                                    models.pa,
                                    metric.eval,
                                    weights,
                                    prevalence,
                                    var.import,
                                    params.CV,
                                    params.OPT,
                                    em.algo,
                                    params.EM,
                                    seed.val,
                                    nb.cpu)
{
  ## 1. Check dir.name and modeling.id
  if (!dir.exists(dir.name)) {
    stop(paste0("Modeling folder '", dir.name, "' does not exist"))
  }
  if (!is.character(modeling.id) || length(modeling.id) > 1) { stop("modeling.id must be a 'character' of length 1") }
  
  ## 2. Check params.PA
  if (length(resp.name) == 1){
    params.PA <- check.params.PA(params.PA)
  } else {
    if (missing(params.PA)){
      params.PA <- list()
    }
    for (sp in resp.name){
      params.PA[[sp]] <- check.params.PA(params.PA[[sp]])
    }
  }
  
  ## 3. Check modeling parameters
  
  ## 3.1 models
  models <- unique(models)
  models.switch.off <- NULL
  
  ## check if model is supported
  if (data.type == "binary"){
    avail.models.list <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF','RFd', 'SRE', 'XGBOOST')
  } else if (data.type == "ordinal") {
    avail.models.list <- c('CTA', 'FDA', 'GAM', 'GLM', 'MARS', 'RF', 'XGBOOST')
  } else {
    avail.models.list <- c('CTA', 'GAM', 'GBM', 'GLM', 'MARS', 'RF', 'XGBOOST')
  }
  .fun_testIfIn(TRUE, paste0("models with ", data.type, " data type"), models, avail.models.list)
  
  
  ## 3.2 models.pa argument
  if (!is.null(models.pa)) {
    if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
      .fun_testIfInherits(TRUE, "models.pa", models.pa, "list")
      .fun_testIfIn(TRUE, "unlist(models.pa)", unlist(models.pa), colnames(bm.format@PA.table))
      .fun_testIfIn(TRUE, "names(models.pa)", names(models.pa), models)
      if (length(models.pa) != length(models)) {
        mod.miss = models[-which(models %in% names(models.pa))]
        list.miss = rep(list(colnames(bm.format@PA.table)), length(mod.miss))
        names(list.miss) = mod.miss
        models.pa = c(models.pa, list.miss)
        warning(paste0(paste0(mod.miss, collapse = ", ")
                       , " have been assigned to all PA datasets as no information was given")
                , immediate. = TRUE)
      }
    } else {
      warning("models.pa has been disabled because no PA datasets have been given", immediate. = TRUE)
      models.pa = NULL
    }
  }
  
  ## 3.3 prevalence arguments --------------------------------------------
  if (!is.null(prevalence)) {
    .fun_testIf01(TRUE, "prevalence", prevalence)
  } else {
    prevalence = 0.5
  }
  
  
  ## 3.5 metric.eval arguments -------------------------------------------
  metric.eval <- unique(metric.eval)
  if (data.type == "binary") {
    avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                              , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC'
                              , 'BOYCE', 'MPA')
  } else if (data.type == "ordinal") {
    avail.eval.meth.list <- c("Accuracy", "Recall", "Precision", "F1")
  } else {
    avail.eval.meth.list <- c('RMSE','MSE',"MAE","Rsquared","Rsquared_aj","Max_error")
  }
  .fun_testIfIn(TRUE, paste0("metric.eval with ", data.type, " data type"), metric.eval, avail.eval.meth.list)
  
  ## 3.6 var.import argument
  if (is.null(var.import)) { var.import = 0 }
  
  
  ## 4. Check params.CV
  if (length(resp.name) == 1){
    params.CV <- check.params.CV(params.CV)
  } else {
    if (missing(params.CV)){
      params.CV <- list()
    }
    for (sp in resp.name){
      params.CV[[sp]] <- check.params.CV(params.CV[[sp]])
    }
  }
  
  ## 5. Check params.OPT
  if (length(resp.name) == 1){
    params.OPT <- check.params.OPT(params.OPT)
  } else {
    if (missing(params.OPT)){
      params.OPT <- list()
    }
    for (sp in resp.name){
      params.OPT[[sp]] <- check.params.OPT(params.OPT[[sp]])
    }
  }
  
  ## 6. Check params.EM
  if (length(resp.name) == 1){
    params.EM <- check.params.EM(params.EM)
  } else {
    if (missing(params.EM)){
      params.EM <- list()
    }
    for (sp in resp.name){
      params.EM[[sp]] <- check.params.EM(params.EM[[sp]])
    }
  }
  
  
  ## 7. Check general parameters
  if (!is.null(seed.val)) { set.seed(seed.val) }
  
  return(list(data.type = data.type,
              filter.raster = filter.raster,
              params.PA = params.PA,
              models = models,
              models.pa = models.pa,
              metric.eval = metric.eval,
              var.import = var.import,
              params.CV = params.CV,
              params.OPT = params.OPT,
              em.algo = em.algo,
              params.EM = params.EM,
              weigths = weights,
              prevalence = prevalence))
}



check.params.PA <- function(params.PA){
  names.params.PA <- c("PA.nb.rep", "PA.nb.absences", "PA.strategy", "PA.dist.min", "PA.dist.max",
                       "PA.sre.quant", "PA.fact.aggr", "PA.user.table")
  if (missing(params.PA) || is.null(params.PA)){
    params.PA <- list(PA.nb.rep = 0,
                      PA.nb.absences = 1000,
                      PA.strategy = NULL,
                      PA.dist.min = 0,
                      PA.dist.max = NULL,
                      PA.sre.quant = 0.025,
                      PA.fact.aggr = NULL,
                      PA.user.table = NULL)
  } else {
    .fun_testIfIn(TRUE, "names of params.PA", names(params.PA), names.params.PA)
    #other check ? They will be done a few moment later anyway ? 
  }
  return(params.PA)
}


check.params.CV <- function(params.CV){
  names.params.CV <- c("CV.nb.rep", "CV.strategy", "CV.perc", "CV.k", "CV.balance",
                       "CV.env.var", "CV.strat", "CV.user.table", "CV.do.full.models")
  if (missing(params.CV) || is.null(params.CV)){
    params.CV <- list(CV.strategy = 'random',
                      CV.nb.rep = 1,
                      CV.perc = 0.7,
                      CV.k = NULL,
                      CV.balance = NULL,
                      CV.env.var = NULL,
                      CV.strat = NULL,
                      CV.user.table = NULL,
                      CV.do.full.models = TRUE)
  } else {
    .fun_testIfIn(TRUE, "names of params.CV", names(params.CV), names.params.CV)
    
    ## 1. Check strategy argument -------------------------------------
    .fun_testIfIn(TRUE, "strategy", params.CV$CV.strategy, c("random", "kfold", "block", "strat", "env", "user.defined"))
    
    ## 2.a Check nb.rep / perc argument -------------------------------
    if (params.CV$CV.strategy %in% c("random", "kfold")) {
      .fun_testIfPosInt(TRUE, "nb.rep", params.CV$CV.nb.rep)
      if (params.CV$CV.nb.rep < 1) { stop("nb.rep must be an integer >= 1") }
      
      if (params.CV$CV.strategy == "random") {
        if (is.null(params.CV$CV.perc)) {
          stop("perc (or CV.perc) is required when strategy = 'random'")
        }
        .fun_testIf01(TRUE, "perc", params.CV$CV.perc)
        if (params.CV$CV.perc < 0.5) {
          warning("You chose to allocate more data to validation than to calibration of your model
                (perc<0.5)\nMake sure you really wanted to do that. \n", immediate. = TRUE)
        } else if (params.CV$CV.perc == 1) {
          params.CV$CV.nb.rep <- 0
          warning(paste0("The models will be evaluated on the calibration data only "
                         , "(nb.rep=0 and no independent data) \n\t "
                         , "It could lead to over-optimistic predictive performances.\n")
                  , immediate. = TRUE)
        }
      }
    }
    
    ## 2.b Check k argument -------------------------------------------
    if (params.CV$CV.strategy %in% c("kfold", "strat", "env")) {
      .fun_testIfPosInt(TRUE, "k", params.CV$CV.k)
      if (params.CV$CV.k < 2) { stop("k must be an integer >= 2") }
    }
    
    ## 2.c Check env.var argument -------------------------------------------
    if (params.CV$CV.strategy %in% c("env")) {
      if (is.null(env.var)) {
        env.var <- names(expl.var)
      } else {
        .fun_testIfIn(TRUE, "env.var", params.CV$env.var, names(expl.var))
      }
    }
    ## 3. Check balance / strat argument ------------------------------
    if (params.CV$CV.strategy %in% c("strat", "env")) {
      .fun_testIfIn(TRUE, "balance", params.CV$CV.balance, c("presences","absences"))
      
      if (params.CV$CV.strategy == "strat") {
        .fun_testIfIn(TRUE, "strat", params.CV$CV.strat, c("x", "y", "both"))
      }
    }
    
    ## 4. Check user.table argument -----------------------------------
    if (params.CV$CV.strategy == "user.defined") {
      if (is.null(params.CV$CV.user.table)) {
        stop("user.table must be a matrix or a data.frame") 
      } else {
        .fun_testIfInherits(TRUE, "user.table", params.CV$CV.user.table, c("matrix", "data.frame"))
        if (inherits(params.CV$CV.user.table, 'data.frame')) {
          params.CV$CV.user.table <- as.matrix(params.CV$CV.user.table)
        }
      }
    }
  }
  return(params.CV)
}


check.params.OPT <- function(params.OPT ){
  names.params.OPT <- c("OPT.strategy", "OPT.user.val", "OPT.user.base", "OPT.user")
  if (missing(params.OPT) || is.null(params.OPT)){
    params.OPT <- list(OPT.strategy = 'bigboss',
                       OPT.user.val = NULL,
                       OPT.user.base = 'bigboss',
                       OPT.user = NULL)
  } else {
    .fun_testIfIn(TRUE, "names of params.OPT", names(params.OPT), names.params.OPT)
    
    .fun_testIfIn(TRUE, "OPT.strategy", params.OPT$OPT.strategy, c("bigboss", "default", "tuned", "user.defined"))
    
    if (!is.null(params.OPT$OPT.user)) {
      .fun_testIfInherits(TRUE, "OPT.user", params.OPT$OPT.user, "BIOMOD.models.options")
    }
  } 
  return(params.OPT)
}

check.params.EM <- function(params.EM){
  names.params.EM <- c("models.chosen", "em.by", "metric.select", "metric.select.thresh", "metric.select.table",
                       "metric.select.dataset", "EMci.alpha", "EMwmean.decay")
  if (missing(params.EM) || is.null(params.EM)){
    params.EM <- list(models.chosen = 'all',
                      em.by = 'PA+run',
                      metric.select = 'all',
                      metric.select.thresh = NULL,
                      metric.select.table = NULL,
                      metric.select.dataset = NULL,
                      EMci.alpha = 0.05,
                      EMwmean.decay = 'proportional')
  } else {
    .fun_testIfIn(TRUE, "names of params.EM", names(params.EM), names.params.EM)
  }
  return(params.EM)
}


.inversion_species_params <- function(liste){
  names_liste <- names(liste)
  names_arguments <- names(unlist(liste))
  for (n in names_liste){
    names_arguments <- sub(paste0(n,"."), "",names_arguments)
  }
  names_arguments <- unique(names_arguments)
  new <- list()
  for (a in names_arguments){
    list_a <- list()
    for (l in names_liste){
      list_a[[l]] <- liste[[l]][[a]]
    }
    new[[a]] <- list_a
  }
  return(new)
}
