
## -------------------------------------------------------------------------- #
## 1. MS.models.out -----------------------------------------------------
## -------------------------------------------------------------------------- #

##' @name MS.models.out
##' @author Hélène Blancheteau
##' 
##' @title \code{MS_Modeling()} output object class
##' 
##' @export

setClass("MS.models.out",
         representation(ms.project = 'character',
                        modeling.id = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        data.type = 'character',
                        expl.var.names = 'character',
                        models.computed = 'list',
                        models.failed = 'list',
                        has.evaluation.data = 'logical',
                        scale.models = 'logical',
                        link = 'character'),
         prototype(modeling.id = as.character(format(Sys.time(), "%s")),
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   has.evaluation.data = FALSE,
                   scale.models = TRUE),
         validity = function(object){ return(TRUE) } )

setClass("BIOMOD.stored.ms.models.out",
         representation(val = 'MS.models.out'),
         prototype(val = NULL),
         validity = function(object) { return(TRUE) })


##' @name getters.out
##' @aliases get_species_data
##' @aliases get_eval_data
##' @aliases get_options
##' @aliases get_calib_lines
##' @aliases get_formal_data
##' @aliases get_projected_models
##' @aliases free
##' @aliases get_predictions
##' @aliases get_kept_models
##' @aliases get_built_models
##' @aliases get_evaluations
##' @aliases get_variables_importance
##' 
##' @title Functions to extract informations from \code{\link{MS.models.out}}, 
##' \code{\link{MS.projection.out}} or \code{\link{MS.ensemble.models.out}} objects
##' 
##' @param obj a \code{\link{MS.formated.data}}, \code{\link{MS.models.out}}, \code{\link{MS.projection.out}} or 
##' \code{\link{MS.ensemble.models.out}} object
##' @param sp a \code{character} name of a specific species 
##' @param as.data.frame a \code{logical} defining whether output should be returned as 
##' \code{data.frame} or \code{array} object
##' @param subinfo a \code{character} corresponding to the information to be extracted, must be 
##' among \code{NULL}, \code{expl.var.names}, \code{resp.var}, \code{expl.var}, \code{MinMax}, 
##' \code{eval.resp.var}, \code{eval.expl.var} (see Details)
##' @param evaluation a \code{logical} defining whether evaluation data should be used or not
##' @param full.name (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing model names to be kept, must be either \code{all} or a 
##' sub-selection of model names that can be obtained with the \code{\link{get_built_models}} 
##' function
##' 
##' @param PA (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing pseudo-absence set to be loaded, must be among \code{PA1}, 
##' \code{PA2}, \code{...}, \code{allData}
##' @param run (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing repetition set to be loaded, must be among \code{RUN1}, 
##' \code{RUN2}, \code{...}, \code{allRun}
##' @param algo (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} containing algorithm to be loaded, must be either 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}
##' 
##' @param merged.by.PA (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing merged pseudo-absence set to be loaded, must be among \code{PA1}, 
##' \code{PA2}, \code{...}, \code{mergedData}
##' @param merged.by.run (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing merged repetition set to be loaded, must be among \code{RUN1}, 
##' \code{RUN2}, \code{...}, \code{mergedRun}
##' @param merged.by.algo (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} containing merged algorithm to be loaded, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}, \code{mergedAlgo}
##' @param filtered.by (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric selected to filter single models to build the 
##' ensemble models, must be among \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, 
##' \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, 
##' \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' 
##' @param metric.eval (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric to be kept, must be among \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, 
##' \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param expl.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing explanatory variables to be kept, that can be obtained with the 
##' \code{\link{get_formal_data}(obj, subinfo = 'expl.var.names')} function
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric selected to transform predictions into binary 
##' values, must be among \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
##' \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, 
##' \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric to filter predictions, must be among \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, 
##' \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' 
##' @param model.as.col (\emph{optional, default} \code{FALSE}) \cr
##' A \code{boolean} given to \code{\link{get_predictions}}. If \code{TRUE} 
##' prediction are returned as a wide \code{data.frame} with each column containing
##' predictions for a single model and corresponding to the old output given by
##' \pkg{biomod2} in version < 4.2-2. If \code{FALSE} predictions are returned 
##' as a long \code{data.frame} with many additional informations readily 
##' available.
##' @param ... additional parameters

setGeneric("get_species_data", function(obj, ...) { standardGeneric("get_species_data") }) ## 012
setGeneric("get_eval_data", function(obj, ...) { standardGeneric("get_eval_data") }) ## 012

setGeneric("get_options", function(obj, ...) { standardGeneric("get_options") }) ## A
setGeneric("get_calib_lines", function(obj, ...) { standardGeneric("get_calib_lines") }) ## A

setGeneric("get_projected_models", function(obj, ...) { standardGeneric("get_projected_models") }) ## B
setGeneric("free", function(obj, ...) { standardGeneric("free") }) ## B

setGeneric("get_predictions", function(obj, ...) { standardGeneric("get_predictions") }) ## ABC

setGeneric("get_kept_models", function(obj, ...) { standardGeneric("get_kept_models") }) ## C

setGeneric("get_formal_data", function(obj, ...) { standardGeneric("get_formal_data") }) ## AC
setGeneric("get_built_models", function(obj, ...) { standardGeneric("get_built_models") }) ## AC
setGeneric("get_evaluations", function(obj, ...) { standardGeneric("get_evaluations") }) ## AC
setGeneric("get_variables_importance", function(obj, ...) { standardGeneric("get_variables_importance") }) ## AC


## show.MS.models.out ---------------------------------------------------
##' 
##' @rdname MS.models.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('MS.models.out'), function(object) {
  .bm_cat("MS.models.out")
  cat("\nModeling folder :", object@dir.name, fill = .Options$width)
  cat("\nSpecies modeled :", object@sp.name, fill = .Options$width)
  cat("\nModeling id :", object@modeling.id, fill = .Options$width)
  cat("\nConsidered variables :", object@expl.var.names, fill = .Options$width)
  cat("\n\nComputed Models : ", object@models.computed, fill = .Options$width)
  .bm_cat()
})

## get_options.MS.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_options", "MS.models.out", function(obj, sp) {
  nameFolder <- file.path(obj@dir.name, sp)
  model <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".models.out"))))
  model_options <- load_stored_object(model@models.options) 
  return(model_options)
})

## get_calib_lines.MS.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##'

setMethod("get_calib_lines", "MS.models.out",
          function(obj, as.data.frame = FALSE, PA = NULL, run = NULL)
          {
            out_total <- foreach (sp = obj@sp.name, .combine = rbind) %do% {
              nameFolder <- file.path(obj@dir.name, sp)
              model <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".models.out"))))
              out <- load_stored_object(model@calib.lines)
              
              if (!is.null(out) && as.data.frame == TRUE) {
                tmp <- melt(out, varnames = c("points", "PA_run"))
                tmp$PA = strsplit(sub("^_", "", tmp$PA_run), "_")[[1]][1]
                tmp$run = strsplit(sub("^_", "", tmp$PA_run), "_")[[1]][2]
                out <- tmp[, c("PA", "run", "points", "value")]
                colnames(out)[4] = "calib.lines"
                
                keep_lines <- .filter_outputs.df(out, subset.list = list(PA = PA, run = run))
                out <- out[keep_lines, ]
              }
              names(out) <- paste(names(out), sp, sep = "_")
              return(out)
            }
            return(out_total)
          }
)

## get_formal_data.MS.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_formal_data", "MS.models.out",
          function(obj, sp, subinfo = NULL)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            model <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".models.out"))))
            
            if (is.null(subinfo)) {
              return(load_stored_object(model@formated.input.data))
            } else if (subinfo == 'MinMax') {
              env = as.data.frame(get_formal_data(model)@data.env.var)
              MinMax = foreach(i = 1:ncol(env)) %do%
                {
                  x = env[, i]
                  if (is.numeric(x)) {
                    return(list(min = min(x, na.rm = TRUE)
                                , max = max(x, na.rm = TRUE)))
                  } else if (is.factor(x)) {
                    return(list(levels = levels(x)))
                  }
                }
              names(MinMax) = colnames(env)
              return(MinMax)
            } else if (subinfo == 'expl.var') {
              return(as.data.frame(get_formal_data(model)@data.env.var))
            } else if (subinfo == 'expl.var.names') {
              return(model@expl.var.names)
            } else if (subinfo == 'resp.var') {
              return(get_formal_data(model)@data.species)
            } else if (subinfo == 'eval.resp.var') {
              return(as.numeric(get_formal_data(model)@eval.data.species))
            } else if (subinfo == 'eval.expl.var') {
              return(as.data.frame(get_formal_data(model)@eval.data.env.var))
            } else { stop("Unknown subinfo tag")}
          }
)



## get_predictions.MS.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "MS.models.out",
          function(obj, sp, evaluation = FALSE
                   , full.name = NULL, PA = NULL, run = NULL, algo = NULL
                   , model.as.col = FALSE)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            model <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".models.out"))))
            
            if (evaluation && (!model@has.evaluation.data)) {
              warning("!   Calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if (evaluation) {
              out <- load_stored_object(mod@models.prediction.eval)
            } else { 
              out <- load_stored_object(mod@models.prediction)
            }
            
            # subselection of models_selected
            keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name, PA = PA
                                                                     , run = run, algo = algo))
            out <- out[keep_lines, ]
            if (model.as.col) {
              out <- .transform_model.as.col(out)
            }
            return(out)
          }
)

## get_built_models.MS.models.out ---------------------------------------------------
##' @rdname getters.out
##' @export
##' 

setMethod("get_built_models", "MS.models.out",
          function(obj, full.name = NULL, PA = NULL, run = NULL, algo = NULL)
          { 
            out <- obj@models.computed
            for (sp in obj@sp.name){
              vec <- out[[sp]]
              keep_ind <- .filter_outputs.vec(vec, obj.type = "mod", subset.list = list(full.name = full.name, PA = PA
                                                                                      , run = run, algo = algo))
              out[[sp]] <- vec[keep_ind]
            }
            return(out)
          }
)

## get_evaluations.MS.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_evaluations", "MS.models.out",
          function(obj, sp, full.name = NULL, PA = NULL, run = NULL, algo = NULL, metric.eval = NULL)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            model <- get(load(file.path(nameFolder, paste0(sp, ".", obj@modeling.id, ".models.out"))))
            out <- load_stored_object(model@models.evaluation)
            if (nrow(out) == 0) {
              cat("\n! models have no evaluations\n")
              return(invisible(NULL))
            } else {
              keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name, PA = PA
                                                                       , run = run, algo = algo
                                                                       , metric.eval = metric.eval))
              out <- out[keep_lines, ]
              return(out)
            }
          }
)

## get_variables_importance.MS.models.out ---------------------------------------------------
##' @rdname getters.out
##' @export
##' 

setMethod("get_variables_importance", "MS.models.out",
          function(obj, sp, full.name = NULL, PA = NULL, run = NULL, algo = NULL, expl.var = NULL)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            model <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".models.out"))))
            out <- load_stored_object(model@variables.importance)
            if (mod@variables.importance@link == '') {
              cat("\n! models have no variables importance\n")
              return(invisible(NULL))
            } else {
              keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name, PA = PA
                                                                       , run = run, algo = algo
                                                                       , expl.var = expl.var))
              out <- out[keep_lines, ]
              return(out)
            }
          }
)
