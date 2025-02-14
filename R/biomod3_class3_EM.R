## --------------------------------------------------------------------------- #
## MS.ensemble.models.out ---------------------------------------------
## --------------------------------------------------------------------------- #

##' @name MS.ensemble.models.out
##' @aliases MS.ensemble.models.out-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_EnsembleModeling()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_EnsembleModeling}}, and used by 
##' \code{\link{BIOMOD_LoadModels}}, \code{\link{BIOMOD_PresenceOnly}} and 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' 
##' @slot ms.project a \code{character} corresponding to the name of the multispecies project
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the
##'   simulation set
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory
##'   variables
##' @slot data.type a \code{character} corresponding to the data type
##' @slot models.out a \code{\link{BIOMOD.stored.models.out-class}} object
##'   containing informations from \code{\link{BIOMOD_Modeling}} object
##' @slot em.by a \code{character} corresponding to the way kept models have
##'   been combined to build the ensemble models, must be among
##'   \code{PA+run}, \code{PA+algo}, \code{PA},
##'   \code{algo}, \code{all}
##' @slot em.computed a \code{vector} containing names of ensemble models
##' @slot em.failed a \code{vector} containing names of failed ensemble models
# ##' @slot em.models_needed a \code{list} containing single models for each ensemble model
##' @slot em.models_kept a \code{list} containing single models for each ensemble model
##' @slot link a \code{character} containing the file name of the saved object
##' 
##' @param object a \code{\link{MS.ensemble.models.out}} object
##' 
##' 
##' @rdname MS.ensemble.models.out
##' @export
##' 

#  Class Definition ---------------------------------------------------------

setClass("MS.ensemble.models.out",
         representation(ms.project = 'character',
                        modeling.id = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        data.type = 'character',
                        models.out = 'BIOMOD.stored.ms.models.out',
                        em.by = 'character',
                        em.computed = 'list',
                        em.failed = 'list',
                        em.models_kept = 'list',
                        link = 'character'),
         prototype(modeling.id = '.',
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   data.type = '',
                   models.out = new('BIOMOD.stored.ms.models.out'),
                   em.by = character(),
                   em.computed = NULL,
                   em.failed = NULL,
                   em.models_kept = NULL),
         validity = function(object){ return(TRUE) })


#  Other functions ----------------------------------------------------------
## show.MS.ensemble.models.out ---------------------------------------------
##' 
##' @rdname MS.ensemble.models.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('MS.ensemble.models.out'), function(object) {
  .bm_cat("MS.ensemble.models.out")
  cat("\nsp.name :", object@sp.name, fill = .Options$width)
  cat("\nexpl.var.names :", object@expl.var.names, fill = .Options$width)
  cat("\n")
  cat("\nmodels computed:", fill = .Options$width)
  show(object@em.computed)
  cat("\nmodels failed:", fill = .Options$width)
  show(object@em.failed)
  .bm_cat()
})

## get_formal_data.MS.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_formal_data", "MS.ensemble.models.out",
          function(obj, sp, subinfo = NULL)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            em <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".ensemble.models.out"))))
            if (is.null(subinfo)) {
              return(load_stored_object(em@models.out))
            } else {
              bm_form = get_formal_data(em)
              return(get_formal_data(bm_form, subinfo = subinfo))
            }
          }
)


## get_built_models.MS.ensemble.models.out ---------------------------------
##'
##' @rdname getters.out
##' @export
##' 

setMethod("get_built_models", "MS.ensemble.models.out",
          function(obj, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL)
          {
            out <- obj@em.computed
            keep_ind <- .filter_outputs.vec(out, obj.type = "em", subset.list = list(full.name = full.name
                                                                                     , merged.by.PA = merged.by.PA
                                                                                     , merged.by.run = merged.by.run
                                                                                     , merged.by.algo = merged.by.algo
                                                                                     , filtered.by = filtered.by
                                                                                     , algo = algo))
            out <- out[keep_ind]
            return(out)
          }
)


## get_kept_models.MS.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_kept_models", "MS.ensemble.models.out", function(obj) { return(obj@em.models_kept) })


## get_predictions.MS.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "MS.ensemble.models.out",
          function(obj, sp, evaluation = FALSE, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL,
                   model.as.col = FALSE)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            em <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".ensemble.models.out"))))
             # check evaluation data availability
            if (evaluation && (!get_formal_data(em)@has.evaluation.data)) {
              warning("!   Calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if (evaluation) { 
              out <- load_stored_object(em@models.prediction.eval)
            } else { 
              out <- load_stored_object(em@models.prediction)
            }
            
            # subselection of models_selected
            keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name
                                                                     , merged.by.algo = merged.by.algo
                                                                     , merged.by.run = merged.by.run
                                                                     , merged.by.PA = merged.by.PA
                                                                     , filtered.by = filtered.by
                                                                     , algo = algo))
            out <- out[keep_lines, ]
            if (model.as.col) {
              out <- .transform_model.as.col(out)
            }
            return(out)
          }
)


## get_evaluations.MS.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_evaluations", "MS.ensemble.models.out",
          function(obj, sp, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL, metric.eval = NULL)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            em <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".ensemble.models.out"))))
            out <- load_stored_object(em@models.evaluation)
            if (nrow(out) == 0) {
              cat("\n! models have no evaluations\n")
              return(invisible(NULL))
            } else {
              keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name
                                                                       , merged.by.algo = merged.by.algo
                                                                       , merged.by.run = merged.by.run
                                                                       , merged.by.PA = merged.by.PA
                                                                       , filtered.by = filtered.by
                                                                       , algo = algo
                                                                       , metric.eval = metric.eval))
              out <- out[keep_lines, ]
              return(out)
            }
          }
)

## get_variables_importance.MS.ensemble.models.out -------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_variables_importance", "MS.ensemble.models.out",
          function(obj, sp, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL, expl.var = NULL)
          {
            nameFolder <- file.path(obj@dir.name, sp)
            em <- get(load(file.path(nameFolder, paste0(sp,".", obj@modeling.id, ".ensemble.models.out"))))
            out <- load_stored_object(em@variables.importance)
            keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name
                                                                     , merged.by.algo = merged.by.algo
                                                                     , merged.by.run = merged.by.run
                                                                     , merged.by.PA = merged.by.PA
                                                                     , filtered.by = filtered.by
                                                                     , algo = algo
                                                                     , expl.var = expl.var))
            out <- out[keep_lines, ]
            return(out)
          }
)


