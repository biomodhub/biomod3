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
##' @slot formated.data a \code{BIOMOD.formated.data} object
##' @slot single.models a \code{BIOMOD.models.out} object
##' @slot ensemble.models a \code{BIOMOD.ensemble.models.out} object
##' @slot output a \code{character} link to the output file
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

