## --------------------------------------------------------------------------- #
## 1. MS.formated.data ---------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name MS.formated.data
##' @aliases MS.formated.data-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_FormatingData()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
##' \code{\link{bm_Tuning}}, \code{\link{bm_CrossValidation}} and 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @param object a \code{\link{MS.formated.data}} object
##' 
##' @slot data.type a \code{character} corresponding to the data type
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates
##' @slot data.species a \code{vector} containing the species observations (\code{0}, \code{1} or 
##' \code{NA})
##' @slot data.env.var a \code{data.frame} containing explanatory variables
##' @slot data.mask a \code{\link[terra:rast]{SpatRaster}} object containing the mask of the 
##' studied area
##' @slot has.data.eval a \code{logical} value defining whether evaluation data is given
##' @slot has.filter.raster a \code{logical} value defining whether filtering have been done or not
##' @slot eval.coord (\emph{optional, default} \code{NULL}) \cr 
##' A 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates for evaluation data
##' @slot eval.data.species (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the species observations (\code{0}, \code{1} or \code{NA}) for 
##' evaluation data
##' @slot eval.data.env.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{data.frame} containing explanatory variables for evaluation data
##' @slot biomod2.version a \code{character} corresponding to the biomod2 version
##' 
##' 
##' @name MS.formated.data-class
##' @rdname MS.formated.data
##' @importFrom terra rast app is.factor subset extract cellFromXY `add<-` 
##' classify rasterize values
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("MS.formated.data",
         representation(data.type = 'character',
                        dir.name = 'character',
                        ms.project = 'character',
                        sp.name = 'character',
                        coord = "data.frame",
                        data.species = "ANY",
                        data.env.var = "data.frame",
                        data.mask = "list",
                        PA.strategy = "character",
                        PA.table = "list",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame",
                        has.filter.raster = "logical",
                        biomod3.version = "character"),
         validity = function(object) { 
           check.data.mask <- suppressWarnings(
             all(sapply(object@data.mask, function(x) inherits(x, "PackedSpatRaster")))
           )
           if (check.data.mask) {
             return(TRUE)
           } else {
             return(FALSE)
           }
         })


setMethod('plot', signature(x = 'MS.formated.data', y = "missing"),
          function(x,
                   calib.lines = NULL,
                   plot.type,
                   plot.output, 
                   PA,
                   run,
                   plot.eval,
                   point.size = 1.5,
                   do.plot = TRUE){
          #Plot adapté
          })


setMethod('show', signature('MS.formated.data'),
          function(object)
          {
            .bm_cat("MS.formated.data")
            cat("\ndir.name = ", object@dir.name, fill = .Options$width)
            cat("\nname of the project = ", object@ms.project, fill = .Options$width)
            cat("\nnumber of species = ", length(object@sp.name), fill = .Options$width)
            
            ## Summary of data.species ?? Peut être très lourd
            cat("\n\n\t",
                ncol(object@data.env.var),
                'explanatory variables\n',
                fill = .Options$width)
            print(summary(object@data.env.var))
            
            if (object@has.data.eval) {
              cat("\n\nEvaluation data :", fill = .Options$width)
              # cat("\n\t",
              #     sum(object@eval.data.species, na.rm = TRUE),
              #     'presences, ',
              #     sum(object@eval.data.species == 0, na.rm = TRUE),
              #     'true absences and ',
              #     sum(is.na(object@eval.data.species), na.rm = TRUE),
              #     'undefined points in dataset',
              #     fill = .Options$width)
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            .bm_cat()
          }
)


setMethod('summary', signature(object = 'MS.formated.data'),
          function(object, calib.lines = NULL) {
            nameFolder <- file.path(object@dir.name, object@ms.project, ".BIOMOD_DATA", "single.formated.data")
            table <- foreach(sp = object@sp.name, .combine = rbind) %do% {
              sfd <- get(load(file.path(nameFolder, paste0(sp,".sfd"))))
              return(data.frame("species" = sp , summary(sfd, calib.lines = calib.lines)))
            }
            table
          }
)

