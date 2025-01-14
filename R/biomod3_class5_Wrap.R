## --------------------------------------------------------------------------- #
## BIOMOD.wrap.out ---------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.wrap.out
##' @aliases BIOMOD.wrap.out-class
##' @author Helene Blancheteau
##' 
##' @title \code{BIOMOD_Wrap()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_Wrap}}, and used by 
##' \code{\link{BIOMOD_ProjectionWrap}}
##' 
##' 
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the
##'   simulation set
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory
##'   variables
##' @slot data.type a \code{character} corresponding to the data type
##' @slot models.out a \code{\link{BIOMOD.stored.models.out-class}} object
##'   containing informations from \code{\link{BIOMOD_Modeling}} object
##' 
##' 
##' @rdname BIOMOD.wrap.out
##' @export
##' 

#  Class Definition ---------------------------------------------------------

setClass("BIOMOD.wrap.out",
         representation(formated.data = "BIOMOD.formated.data",
                        single.models = "BIOMOD.models.out",
                        ensemble.models = "BIOMOD.ensemble.models.out",
                        output = "character"),
         validity = function(object){ return(TRUE) })


#  Other functions ----------------------------------------------------------
## show.BIOMOD.wrap.out ---------------------------------------------
##' 
##' @rdname BIOMOD.wrap.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.wrap.out'), function(object) {
  .bm_cat("BIOMOD.wrap.out")
  cat("\nsp.name :", object@formated.data@sp.name, fill = .Options$width)
  cat("\nexpl.var.names :", object@single.models@expl.var.names, fill = .Options$width)
  cat("\n")
  cat("\nSingle models computed:", toString(object@single.models@models.computed), fill = .Options$width)
  cat("\nEnsemble models computed:", toString(object@ensemble.models@em.computed), fill = .Options$width)
  .bm_cat()
})

