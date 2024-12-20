###################################################################################################
##' @name BIOMOD_ProjectionWrap
##' @author Helene Blancheteau
##' 
##' @title Project a range of calibrated species distribution models onto new environment
##' 
##' @description This function allows to project a range of models built with the 
##' \code{\link{BIOMOD_Wrap}} function onto new environmental data (\emph{which can 
##' represent new areas, resolution or time scales for example}).
##' 
##' 
##' @param bm.wrap a \code{\link{BIOMOD.wrap.out}} object returned by the 
##' \code{\link{BIOMOD_Modeling}} function
##' @param proj.name a \code{character} corresponding to the name (ID) of the projection set 
##' (\emph{a new folder will be created within the simulation folder with this name})
##' @param new.env A \code{matrix}, \code{data.frame} or
##'  \code{\link[terra:rast]{SpatRaster}} object containing the new 
##' explanatory variables (in columns or layers, with names matching the
##'  variables names given to the \code{\link{BIOMOD_FormatingData}} function to build 
##' \code{bm.wrap}) that will be used to project the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} are still supported such as 
##' \code{RasterStack} objects. }
##' 
##' @param new.env.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{new.env} is a \code{matrix} or a \code{data.frame}, a 2-columns \code{matrix} or 
##' \code{data.frame} containing the corresponding \code{X} and \code{Y} coordinates that will be 
##' used to project the species distribution model(s)
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all}, \code{single}, \code{ensemble} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into binary values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{bm.wrap}) or \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
##' \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, 
##' \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into filtered values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{bm.wrap}) or \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
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
##' 
##' 
##' @keywords models projection
##' 
##' @export
##' 
##' 
###################################################################################################

BIOMOD_ProjectionWrap <- function(bm.wrap,
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
  .bm_cat("Do Single Models Projection")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_ProjectionWrap.check.args(bm.wrap, proj.name, new.env, new.env.xy
                                        , models.chosen, metric.binary, metric.filter, compress, seed.val, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Create output object ----------------------------------------------------------------------
  proj_out <- new('BIOMOD.projection.out',
                  proj.name = proj.name,
                  dir.name = bm.wrap@formated.data@dir.name,
                  sp.name =  bm.wrap@formated.data@sp.name,
                  expl.var.names = bm.wrap@formated.data@expl.var.names,
                  #models.projected = models.chosen,
                  scale.models = bm.wrap@single.models@scale.models,
                  coord = new.env.xy,
                  data.type = bm.wrap@single.models@data.type,
                  modeling.id = bm.wrap@single.models@modeling.id,
                  type = "wrap")
  
  
  cat("\n")
  cat("\n\t Projection of single models")

  output <- capture.output(proj_single <- BIOMOD_Projection(bm.wrap@single.models,
                                                        proj.name = proj.name,
                                                        new.env = new.env,
                                                        new.env.xy = new.env.xy,
                                                        models.chosen = models.chosen.single,
                                                        metric.binary = metric.binary,
                                                        metric.filter = metric.filter,
                                                        compress = compress,
                                                        build.clamping.mask = build.clamping.mask,
                                                        nb.cpu = nb.cpu,
                                                        digits = digits))
  
  cat("\n\t Projection of ensemble models")
  
  output <- capture.output(proj_ens <- BIOMOD_Projection(bm.wrap@ensemble.models,
                                                            proj.name = proj.name,
                                                            new.env = new.env,
                                                            new.env.xy = new.env.xy,
                                                            models.chosen = models.chosen.ens,
                                                            metric.binary = metric.binary,
                                                            metric.filter = metric.filter,
                                                            compress = compress,
                                                            build.clamping.mask = build.clamping.mask,
                                                            nb.cpu = nb.cpu,
                                                            digits = digits))
  
  
    
  
  proj_out@models.projected <- c(proj_single@models.projected, proj_ens@models.projected)
  # proj_out@models.out <- fusion !!!!
  # proj_out@proj.out <- fusion !!!
  
  .bm_cat("Done")
  return(proj_out)
}

# .BIOMOD_ProjectionWrap.check.args---------------------------------------------

.BIOMOD_ProjectionWrap.check.args <- function(bm.wrap, proj.name, new.env, new.env.xy,
                                          models.chosen, metric.binary, metric.filter, compress, seed.val, ...)
{
  args <- list(...)
  
  ## 1. Check bm.wrap ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.wrap", bm.wrap, "BIOMOD.wrap.out")
  
  ## 2. Check proj.name -------------------------------------------------------
  if (is.null(proj.name)) {
    stop("\nYou must define a name for Projection Outputs")
  } 
  
  ## 3. Check new.env ---------------------------------------------------------
  .fun_testIfInherits(TRUE, "new.env", new.env, c('matrix', 'data.frame', 'SpatRaster','Raster'))
  
  if (inherits(new.env, 'matrix')) {
    if (any(sapply(get_formal_data(bm.wrap, sp = bm.wrap@sp.name[1], "expl.var"), is.factor))) {
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
                                             expected_levels = head(get_formal_data(bm.wrap, sp = bm.wrap@sp.name[1], subinfo = "expl.var"))
      )
    } else {
      new.env <- rast(new.env)
    }
  }
  
  if (inherits(new.env, 'SpatRaster')) {
    .fun_testIfIn(TRUE, "names(new.env)", names(new.env), bm.wrap@expl.var.names, exact = TRUE)
    new.env <- new.env[[bm.wrap@expl.var.names]]
    new.env.mask <- .get_data_mask(new.env, value.out = 1)
    new.env <- mask(new.env, new.env.mask)
  } else {
    .fun_testIfIn(TRUE, "colnames(new.env)", colnames(new.env), bm.wrap@expl.var.names, exact = TRUE)
    new.env <- new.env[ , bm.wrap@expl.var.names, drop = FALSE]
  }
  
  which.factor <- which(sapply(new.env, is.factor))
  if (length(which.factor) > 0) {
    new.env <- .check_env_levels(new.env, 
                                 expected_levels = head(get_formal_data(bm.wrap, sp = bm.wrap@sp.name[1], subinfo = "expl.var")))
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
    models.chosen <- as.list(rep("all", length(bm.wrap@sp.name)))
    names(models.chosen) <- bm.wrap@sp.name
  } else {
    .fun_testIfInherits(TRUE, "models.chosen", models.chosen, "list")
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  
  ## 6. Check metric.binary & metric.filter -----------------------------------
  if (!is.null(metric.binary) | !is.null(metric.filter)) {
    if ( bm.wrap@data.type != "binary"){
      cat ("No metric.binary or metric.filter are needed with",bm.wrap@data.type, "data")
      metric.binary <- NULL
      metric.filter <- NULL
    } else {
      models.evaluation <- get_evaluations(bm.wrap, sp = bm.wrap@sp.name[1])
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
  if (bm.wrap@data.type %in% c("count","abundance","ordinal")) {on_0_1000 <- FALSE}
  
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
