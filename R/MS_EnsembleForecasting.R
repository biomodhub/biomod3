###################################################################################################
##' @name MS_EnsembleForecasting
##' @author Helene Blancheteau
##' 
##' @title Project a range of calibrated species distribution models onto new environment
##' 
##' @description This function allows to project a range of models built with the 
##' \code{\link{BIOMOD_Modeling}} function onto new environmental data (\emph{which can 
##' represent new areas, resolution or time scales for example}).
##' 
##' 
##' @param ms.mod a \code{\link{BIOMOD.models.out}} object returned by the 
##' \code{\link{BIOMOD_Modeling}} function
##' @param proj.name a \code{character} corresponding to the name (ID) of the projection set 
##' (\emph{a new folder will be created within the simulation folder with this name})
##' @param new.env A \code{matrix}, \code{data.frame} or
##'  \code{\link[terra:rast]{SpatRaster}} object containing the new 
##' explanatory variables (in columns or layers, with names matching the
##'  variables names given to the \code{\link{BIOMOD_FormatingData}} function to build 
##' \code{ms.mod}) that will be used to project the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} are still supported such as 
##' \code{RasterStack} objects. }
##' 
##' @param new.env.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{new.env} is a \code{matrix} or a \code{data.frame}, a 2-columns \code{matrix} or 
##' \code{data.frame} containing the corresponding \code{X} and \code{Y} coordinates that will be 
##' used to project the species distribution model(s)
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into binary values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{ms.mod}) or \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
##' \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, 
##' \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into filtered values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{ms.mod}) or \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
##' \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, 
##' \code{ETS}, \code{BOYCE}, \code{MPA}
##' 
##' @param compress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} or a \code{character} value defining whether and how objects should be 
##' compressed when saved on hard drive. Must be either \code{TRUE}, \code{FALSE}, \code{xz} or 
##' \code{gzip} (see Details)
##' @param digits (\emph{optional, default} \code{0}) \cr 
##' A \code{integer} value defining the number of digits of the predictions. 
##' @param build.clamping.mask (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether a clamping mask should be built and saved on hard 
##' drive or not (see Details)
##' 
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' 
##' @param \ldots (\emph{optional, see Details)}) 
##' 
##' 
##' @return
##' 
##' A \code{\link{BIOMOD.projection.out}} object containing models projections, or links to saved 
##' outputs. \cr Models projections are stored out of \R (for memory storage reasons) in 
##' \code{proj.name} folder created in the current working directory :
##' \enumerate{
##'   \item the output is a \code{data.frame} if \code{new.env} is a \code{matrix} or a 
##'   \code{data.frame}
##'   \item it is a \code{\link[terra:rast]{SpatRaster}} if \code{new.env} is a 
##'   \code{\link[terra:rast]{SpatRaster}} (or several \code{\link[terra:rast]{SpatRaster}} 
##'   objects, if \code{new.env} is too large)
##'   \item raw projections, as well as binary and filtered projections (if asked), are saved in 
##'   the \code{proj.name} folder
##' }
##' 
##' 
##' @details 
##' 
##' If \code{models.chosen = 'all'}, projections are done for all calibration and pseudo absences 
##' runs if applicable. \cr These projections may be used later by the 
##' \code{\link{BIOMOD_EnsembleForecasting}} function. \cr \cr 
##' 
##' If \code{build.clamping.mask = TRUE}, a raster file will be saved within the projection folder. 
##' This mask values will correspond to the number of variables in each pixel that are out of their 
##' calibration / validation range, identifying locations where predictions are uncertain. \cr \cr
##' 
##' \code{...} can take the following values :
##' \itemize{
##   \item \code{clamping.level} : a \code{logical} value defining whether \code{clamping.mask} 
##   cells with at least one variable out of its calibration range are to be removed from the 
##   projections or not
##'   \item \code{omit.na} : a \code{logical} value defining whether all not fully referenced 
##'   environmental points will get \code{NA} as predictions or not
##'   \item \code{on_0_1000} : a \code{logical} value defining whether \code{0 - 1} probabilities 
##'   are to be converted to \code{0 - 1000} scale to save memory on backup
##'   \item \code{do.stack} : a \code{logical} value defining whether all projections are to be 
##'   saved as one \code{\link[terra:rast]{SpatRaster}} object or several 
##'   \code{\link[terra:rast]{SpatRaster}} files (\emph{the default if projections are too heavy to 
##'   be all loaded at once in memory})
##'   \item \code{keep.in.memory} : a \code{logical} value defining whether all projections are 
##'   to be kept loaded at once in memory, or only links pointing to hard drive are to be returned
##'   \item \code{output.format} : a \code{character} value corresponding to the projections 
##'   saving format on hard drive, must be either \code{.grd}, \code{.img}, \code{.tif} or \code{.RData} (the 
##'   default if \code{new.env} is given as \code{matrix} or \code{data.frame})
##'   \item \code{overwrite} : a \code{logical}. FALSE by default. If some projections with the same project ID 
##'   and modeling ID are already create, does biomod2 crush it or not ? 
##' }
##' 
##' 
##' @keywords models projection
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_RangeSize}}
##' @family Main functions
##' 
##' 
##' @importFrom foreach foreach %dopar% 
##' @importFrom terra rast subset nlyr writeRaster terraOptions wrap unwrap
##'  mem_info app is.factor mask
##' @importFrom utils capture.output
##' @importFrom abind asub
##' 
##' @export
##' 
##' 
###################################################################################################

MS_EnsembleForecasting <- function(ms.mod,
                          proj.name,
                          new.env,
                          new.env.xy = NULL,
                          models.chosen = 'all',
                          metric.binary = NULL,
                          metric.filter = NULL,
                          compress = TRUE,
                          digits = 0,
                          build.clamping.mask = TRUE,
                          nb.cpu = 1,
                          seed.val = NULL,
                          ...) {
  .bm_cat("Do Ensemble Models Projection")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .MS_EnsembleForecasting.check.args(ms.mod, proj.name, new.env, new.env.xy
                                    , models.chosen, metric.binary, metric.filter, compress, seed.val, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Create output object ----------------------------------------------------------------------
  proj_out <- new('MS.projection.out',
                  proj.name = proj.name,
                  dir.name = ms.mod@dir.name,
                  sp.name =  ms.mod@sp.name,
                  expl.var.names = ms.mod@expl.var.names,
                  models.projected = models.chosen,
                  coord = new.env.xy,
                  data.type = ms.mod@data.type,
                  modeling.id = ms.mod@modeling.id,
                  type = "em")
  
  
  cat("\n")
  workflow <- foreach(sp = ms.mod@sp.name) %do% {
    cat("\n\t Projection of", sp)
    # 1. Récupération ms.mod 
    bm.em <- get(load(file.path(ms.mod@dir.name, sp, paste0(sp, ".", ms.mod@modeling.id,".ensemble.models.out"))))
    
    models.chosen.sp <- models.chosen[[sp]]
    # 2. Run MS_EnsembleModeling
    output <- capture.output(proj_sp <- BIOMOD_EnsembleForecasting(bm.em,
                                                          proj.name = proj.name,
                                                          new.env = new.env,
                                                          new.env.xy = new.env.xy,
                                                          models.chosen = models.chosen.sp,
                                                          metric.binary = metric.binary,
                                                          metric.filter = metric.filter,
                                                          compress = TRUE,
                                                          build.clamping.mask = TRUE,
                                                          nb.cpu = 1,
                                                          seed.val = NULL))
    models.chosen[[sp]] <- proj_sp@models.projected
    
  }
  proj_out@models.projected <- models.chosen
  
  .bm_cat("Done")
  return(proj_out)
}

# .MS_EnsembleForecasting.check.args---------------------------------------------

.MS_EnsembleForecasting.check.args <- function(ms.mod, proj.name, new.env, new.env.xy,
                                      models.chosen, metric.binary, metric.filter, compress, seed.val, ...)
{
  args <- list(...)
  
  ## 1. Check ms.mod ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "ms.mod", ms.mod, "MS.ensemble.models.out")
  
  ## 2. Check proj.name -------------------------------------------------------
  if (is.null(proj.name)) {
    stop("\nYou must define a name for Projection Outputs")
  } 
  
  ## 3. Check new.env ---------------------------------------------------------
  .fun_testIfInherits(TRUE, "new.env", new.env, c('matrix', 'data.frame', 'SpatRaster','Raster'))
  
  if (inherits(new.env, 'matrix')) {
    if (any(sapply(get_formal_data(ms.mod, sp = ms.mod@sp.name[1], "expl.var"), is.factor))) {
      stop("new.env cannot be given as matrix when model involves categorical variables")
    }
    new.env <- data.frame(new.env)
  } else if (inherits(new.env, 'data.frame')) {
    # ensure that data.table are coerced into classic data.frame
    new.env <- as.data.frame(new.env) 
  }
  
  if (inherits(new.env, 'Raster')) {
    # conversion into SpatRaster
    if (any(raster::is.factor(new.env))) {
      new.env <- .categorical_stack_to_terra(raster::stack(new.env),
                                             expected_levels = head(get_formal_data(ms.mod, sp = ms.mod@sp.name[1], subinfo = "expl.var"))
      )
    } else {
      new.env <- rast(new.env)
    }
  }
  
  if (inherits(new.env, 'SpatRaster')) {
    .fun_testIfIn(TRUE, "names(new.env)", names(new.env), ms.mod@expl.var.names, exact = TRUE)
    new.env <- new.env[[ms.mod@expl.var.names]]
    new.env.mask <- .get_data_mask(new.env, value.out = 1)
    new.env <- mask(new.env, new.env.mask)
  } else {
    .fun_testIfIn(TRUE, "colnames(new.env)", colnames(new.env), ms.mod@expl.var.names, exact = TRUE)
    new.env <- new.env[ , ms.mod@expl.var.names, drop = FALSE]
  }
  
  which.factor <- which(sapply(new.env, is.factor))
  if (length(which.factor) > 0) {
    new.env <- .check_env_levels(new.env, 
                                 expected_levels = head(get_formal_data(ms.mod, sp = ms.mod@sp.name[1], subinfo = "expl.var")))
  }
  
  ## 4. Check new.env.xy ------------------------------------------------------
  if (!is.null(new.env.xy)  & !inherits(new.env, 'SpatRaster')) {
    new.env.xy = as.data.frame(new.env.xy)
    if (ncol(new.env.xy) != 2 || nrow(new.env.xy) != nrow(new.env)) {
      stop("invalid xy coordinates argument given -- dimensions mismatch !")
    }
  } else {
    new.env.xy = data.frame()
  }
  ## 5. Check models.chosen ---------------------------------------------------
  if ( is.null(models.chosen) | (length(models.chosen) == 1 && models.chosen[1] == 'all')) {
    models.chosen <- as.list(rep("all", length(ms.mod@sp.name)))
    names(models.chosen) <- ms.mod@sp.name
  } else {
    .fun_testIfInherits(TRUE, "models.chosen", models.chosen, "list")
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  
  ## 6. Check metric.binary & metric.filter -----------------------------------
  if (!is.null(metric.binary) | !is.null(metric.filter)) {
    if ( ms.mod@data.type != "binary"){
      cat ("No metric.binary or metric.filter are needed with",ms.mod@data.type, "data")
      metric.binary <- NULL
      metric.filter <- NULL
    } else {
      models.evaluation <- get_evaluations(ms.mod, sp = ms.mod@sp.name[1])
      if (is.null(models.evaluation)) {
        warning("Binary and/or Filtered transformations of projection not ran because of models evaluation information missing")
      } else {
        available.evaluation <- unique(models.evaluation$metric.eval)
        if (!is.null(metric.binary) && metric.binary[1] == 'all') {
          metric.binary <- available.evaluation
        } else if (!is.null(metric.binary) && sum(!(metric.binary %in% available.evaluation)) > 0) {
          warning(paste0(toString(metric.binary[!(metric.binary %in% available.evaluation)]),
                         " Binary Transformation were switched off because no corresponding evaluation method found"))
          metric.binary <- metric.binary[metric.binary %in% available.evaluation]
        }
        
        if (!is.null(metric.filter) && metric.filter[1] == 'all') {
          metric.filter <- available.evaluation
        } else if (!is.null(metric.filter) && sum(!(metric.filter %in% available.evaluation)) > 0) {
          warning(paste0(toString(metric.filter[!(metric.filter %in% available.evaluation)]),
                         " Filtered Transformation were switched off because no corresponding evaluation method found"))
          metric.filter <- metric.filter[metric.filter %in% available.evaluation]
        }
      }
    }
  }
  
  ## 7. Check compress --------------------------------------------------------
  if (compress == 'xz') {
    compress <- ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }
  
  ## 9. Check output.format ---------------------------------------------------
  output.format <- args$output.format # raster output format
  if (!is.null(output.format)) {
    if (!output.format %in% c(".img", ".grd", ".tif", ".RData")) {
      stop(paste0("output.format argument should be one of '.img','.grd', '.tif' or '.RData'\n"
                  , "Note : '.img','.grd', '.tif' are only available if you give environmental condition as a SpatRaster object"))
    }
    if (output.format %in% c(".img", ".grd", ".tif") && !inherits(new.env, "SpatRaster")) {
      warning("output.format was automatically set to '.RData' because environmental conditions are not given as a raster object")
    }
  } else {
    output.format <- ifelse(!inherits(new.env, "SpatRaster"), ".RData", ".tif")
  }
  
  ## 9. Check do.stack --------------------------------------------------------
  do.stack <- ifelse(is.null(args$do.stack), TRUE, args$do.stack)
  if (!inherits(new.env, 'SpatRaster')) {
    if (!do.stack) { 
      cat("\n\t\t! 'do.stack' arg is always set as TRUE for data.frame/matrix dataset") 
    }
    do.stack <- TRUE
  } else if (do.stack) { 
    # test if there is enough memory to work with multilayer SpatRaster
    capture.output({
      test <-
        mem_info(
          subset(new.env, 1),
          n = 2 * length(models.chosen[1]) + nlyr(new.env)
        )
    })
    if (test["needed"] >= test["available"]) { 
      terraOptions(todisk = TRUE) 
    }
  } else if (output.format == "RData"){
    cat("\n\t\t! 'do.stack' arg is always set as TRUE for .RData output format") 
    do.stack <- TRUE
  }
  
  ## 10.on_0_1000 --------------------------------
  
  on_0_1000 <- ifelse(is.null(args$on_0_1000), TRUE, args$on_0_1000)
  if (ms.mod@data.type %in% c("count","abundance","ordinal")) {on_0_1000 <- FALSE}
  
  ## 11.Check overwrite
  overwrite <- ifelse(is.null(args$overwrite), ifelse(do.stack, TRUE, FALSE), args$overwrite)
  if (!overwrite){
    cat("\n\t\t! 'overwrite' arg is set as FALSE. Projections that have already been saved will not be redone. 
        Please be carfeul if you have changed the models in the meantime. ")
  }
  
  return(list(proj.name = proj.name,
              new.env = new.env,
              new.env.xy = new.env.xy,
              models.chosen = models.chosen,
              metric.binary = metric.binary,
              metric.filter = metric.filter,
              compress = compress,
              do.stack = do.stack,
              output.format = output.format,
              omit.na = ifelse(is.null(args$omit.na), TRUE, args$omit.na),
              do.stack = do.stack,
              keep.in.memory = ifelse(is.null(args$keep.in.memory), TRUE, args$keep.in.memory),
              on_0_1000 = on_0_1000,
              seed.val = seed.val,
              overwrite = overwrite))
}
