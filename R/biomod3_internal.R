
# .getOS <- function()
# {
#   sysinf = Sys.info()
#   if (!is.null(sysinf))
#   {
#     os = sysinf['sysname']
#     if (os == 'Darwin') os = "mac"
#   } else
#   {
#     os = .Platform$OS.type
#     if (grepl("^darwin", R.version$os)) os = "mac"
#     if (grepl("linux-gnu", R.version$os)) os = "linux"
#   }
#   return(tolower(os))
# }


## BIOMOD SPECIFIC CAT ----------------------------------------------------------------------------

.bm_cat <- function(x = NULL, ...)
{
  if (is.null(x))
  {
    cat("\n")
    cat(paste(rep("-=", round(.Options$width / 2)), collapse = ""))
    cat("\n")
  } else
  {
    x.length = nchar(x) + 2
    y.length = (.Options$width - x.length) / 2
    cat("\n")
    cat(paste(rep("-=", round(y.length / 2)), collapse = "")
        , x
        , paste(rep("-=", round(y.length / 2)), collapse = ""))
    cat("\n")
  }
}


## TEST PARAMETERS --------------------------------------------------------------------------------

.fun_testIfInherits <- function(test, objName, objValue, values)
{
  if (!inherits(objValue, values)) {
    stop(paste0("\n", paste0(objName, " must be a '", paste0(values[1:(length(values) -1)], collapse = "', '")
                             , ifelse(length(values) > 1, paste0("' or '", values[length(values)]), "")
                             , "' object")))
    test <- FALSE
  }
  return(test)
}

.fun_testIfIn <- function(test, objName, objValue, values, exact = FALSE)
{
  if (any(! objValue %in% values) ||
      (exact == TRUE && any(! values %in% objValue))) {
    stop(paste0("\n", objName, " must be '", 
                ifelse(length(values) > 1, 
                       paste0(paste0(values[1:(length(values) -1)], collapse = "', '"),
                              "' or '", values[length(values)])
                       , paste0(values,"'"))))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosNum <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    stop(paste0("\n", objName, " must be a numeric"))
    test <- FALSE
  } else if (objValue < 0) {
    stop(paste0("\n", objName, " must be a positive numeric"))
    test <- FALSE
  }
  return(test)
}

.fun_testIf01 <- function(test, objName, objValue)
{
  test <- .fun_testIfPosNum(test, objName, objValue)
  if (test && objValue > 1) {
    stop(paste0("\n", objName, " must be a 0 to 1 numeric"))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosInt <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    cat(paste0("\n", objName, " must be a integer"))
    test <- FALSE
  } else if (any(objValue < 0) || any(objValue %% 1 != 0)) {
    cat(paste0("\n", objName, " must be a positive integer"))
    test <- FALSE
  }
  return(test)
}

.fun_testMetric <- function(test, objName, objValue, values)
{
  if (!is.null(objValue)) {
    test <- .fun_testIfInherits(test, objName, objValue, "character")
    if(length(objValue) != 1) {
      stop(paste0("`",objName,"` must only have one element"))
    }
    test <- .fun_testIfIn(test, objName, objValue, values)
  }
  return(test)
}




## TOOLS for SpatRaster ---------------------------------------------------------------------------
## used in biomod2_classes_1

##' @name rast.has.values
##' @title Check whether SpatRaster is an empty \code{rast()} 
##' @description Check whether SpatRaster is an empty \code{rast()} 
##' @param x SpatRaster to be checked
##' @return a boolean
##' @keywords internal
##' @importFrom terra values

rast.has.values <- function(x)
{
  stopifnot(inherits(x, 'SpatRaster'))
  suppressWarnings(any(!is.na(values(x))))
}



## used in bm_PseudoAbsences, BIOMOD_Projection, BIOMOD_EnsembleForecasting
.get_data_mask <- function(expl.var, value.out = 1)
{
  stopifnot(inherits(expl.var, 'SpatRaster'))
  classify(as.numeric(!any(is.na(expl.var))),
           matrix(c(0, NA, 
                    1, value.out), 
                  ncol = 2, byrow = TRUE)) 
}



## EXTRACT model names according to specific infos ------------------------------------------------
## used in biomod2_classes_4, BIOMOD_EnsembleModeling

.extract_modelNamesInfo <- function(model.names, obj.type, info, as.unique = TRUE)
{
  ## 0. CHECK parameters ----------------------------------------------------------------
  if (!is.character(model.names)) { stop("model.names must be a character vector") }
  .fun_testIfIn(TRUE, "obj.type", obj.type, c("mod", "em"))
  if (obj.type == "mod") {
    .fun_testIfIn(TRUE, "info", info, c("species", "PA", "run", "algo"))
  } else if (obj.type == "em") {
    .fun_testIfIn(TRUE, "info", info, c("species", "merged.by.PA", "merged.by.run", "merged.by.algo", "filtered.by", "algo"))
  }
  
  ## 1. SPLIT model.names ---------------------------------------------------------------
  info.tmp <- strsplit(model.names, "_")
  
  if (obj.type == "mod") {
    res <- switch(info
                  , species = sapply(info.tmp, function(x) x[1])
                  , PA = sapply(info.tmp, function(x) x[2])
                  , run = sapply(info.tmp, function(x) x[3])
                  , algo = sapply(info.tmp, function(x) x[4]))
  } else if (obj.type == "em") {
    res <- switch(info
                  , species = sapply(info.tmp, function(x) x[1])
                  , merged.by.PA = sapply(info.tmp, function(x) x[3])
                  , merged.by.run = sapply(info.tmp, function(x) x[4])
                  , merged.by.algo = sapply(info.tmp, function(x) x[5])
                  , filtered.by = sapply(info.tmp, function(x) sub(".*By", "", x[2]))
                  , algo = sapply(info.tmp, function(x) sub("By.*", "", x[2])))
  }
  
  ## 2. RETURN either the full vector, or only the possible values ----------------------
  if (as.unique == TRUE) {
    return(unique(res))
  } else {
    return(res)
  }
}



## FORMAT .format_proj.df output from long to wide with models as columns -------------------------
## used in biomod2_classes_3

##' @name .transform_model.as.col
##' @author Remi Patin
##' 
##' @title Transform predictions data.frame from long to wide with models as columns
##' 
##' @description This function is used internally in \code{\link{get_predictions}}
##' to ensure back-compatibility with former output of \code{\link{get_predictions}}
##' (i.e. for \pkg{biomod2} version < 4.2-2). It transform a long \code{data.frame}
##' into a wide \code{data.frame} with each column corresponding to a single model.
##' 
##' \emph{Note that the function is intended for internal use but have been 
##' made available for compatibility with \pkg{ecospat}} 
##'
##' @param df a long \code{data.frame}, generally a standard output of 
##' \code{\link{get_predictions}}
##' 
##' @return a wide \code{data.frame}
##' @importFrom stats reshape
##' @export
##' @keywords internal

.transform_model.as.col <- function(df)
{
  df <- reshape(df[, c("full.name", "points", "pred")], idvar = "points"
                , timevar = "full.name", direction = "wide")[ , -1, drop = FALSE] 
  colnames(df) <- substring(colnames(df), 6)
  df
}


## EXTRACT models outputs according to specific infos ---------------------------------------------
## used in biomod2_classes_3.R

.filter_outputs.df <- function(out, subset.list)
{
  keep_subset <- foreach(sub.i = names(subset.list)) %do%
    {
      keep_lines <- 1:nrow(out)
      if (!is.null(subset.list[[sub.i]])) {
        .fun_testIfIn(TRUE, sub.i, subset.list[[sub.i]], unique(out[, sub.i]))
        keep_lines <- which(out[, sub.i] %in% subset.list[[sub.i]])
      }
      return(keep_lines)
    }
  keep_lines <- Reduce(intersect, keep_subset)
  if (length(keep_lines) == 0) {
    warning(paste0("No information corresponding to the given filter(s) ("
                   , paste0(names(subset.list), collapse = ', '), ")"))
  }
  return(keep_lines)
}

.filter_outputs.vec <- function(out_names, obj.type, subset.list)
{
  keep_layers <- out_names
  if ("full.name" %in% names(subset.list) && !is.null(subset.list[["full.name"]])) {
    .fun_testIfIn(TRUE, "full.name", subset.list[["full.name"]], out_names)
    keep_layers <- which(out_names %in% subset.list[["full.name"]])
  } else {
    keep_subset <- foreach(sub.i = names(subset.list)) %do%
      {
        keep_tmp <- 1:length(out_names)
        if (!is.null(subset.list[[sub.i]])) {
          .fun_testIfIn(TRUE, sub.i, subset.list[[sub.i]], .extract_modelNamesInfo(out_names, obj.type = obj.type, info = sub.i, as.unique = TRUE))
          keep_tmp <- grep(paste0(subset.list[[sub.i]], ifelse(sub.i %in% c("PA", "run"), "_", ""), collapse = "|"), out_names)
        }
        return(keep_tmp)
      }
    keep_layers <- Reduce(intersect, keep_subset)
  }
  if (length(keep_layers) == 0) {
    warning(paste0("No information corresponding to the given filter(s) ("
                   , paste0(names(subset.list), collapse = ', '), ")"))
  }
  return(keep_layers)
}


## EXTRACT projections info and specific layers ---------------------------------------------------
## used in biomod2_classes_3.R

## Extract info (out/binary/filtered + metric) from BIOMOD.projection.out
.extract_projlinkInfo <- function(obj)
{
  stopifnot(inherits(obj, 'BIOMOD.projection.out'))
  out.names <- obj@proj.out@link
  sp.name <- obj@sp.name
  sub(pattern     = paste0(".*?_",sp.name), 
      replacement = "",
      x           = out.names)
  out.binary <- grepl(pattern = "bin\\.",
                      x = out.names)
  out.filt <- grepl(pattern = "filt\\.",
                    x = out.names)
  out.type <- ifelse(out.binary,
                     "bin",
                     ifelse(out.filt,
                            "filt", 
                            "out"))
  # To avoid special case when the species name finish by "bin" (or "filt") ??
  # out.format <- tools::file_ext(out.names)
  # if(grepl(pattern = paste0(sp.name,"/.",out.format), x = out.names)){ 
  #   out.type <- "out"
  # }
  out.metric <- sapply(seq_len(length(out.names)), 
                       function(i){
                         if (out.type[i] == "out") {
                           return("none")
                         } else {
                           begin.cut <-  gsub(".*_", "", out.names[i]) 
                           return(  sub(paste0(out.type[i],".*"), "", begin.cut) )
                         }
                       })
  data.frame("link"   = out.names,
             "type"   = out.type,
             "metric" = out.metric)
}

## Return relevant layers indices from BIOMOD.projection.out given metric.binary or metric.select
.extract_selected.layers <- function(obj, metric.binary, metric.filter)
{
  ## argument  check ----------------------------------------------------------
  stopifnot(inherits(obj, "BIOMOD.projection.out"))
  
  if (!is.null(metric.binary) & !is.null(metric.filter)) {
    stop("cannot return both binary and filtered projection, please provide either `metric.binary` or `metric.filter` but not both.")
  }
  
  df.info <- .extract_projlinkInfo(obj)
  
  available.metrics.binary <- unique(df.info$metric[which(df.info$type == "bin")])
  available.metrics.filter <- unique(df.info$metric[which(df.info$type == "filt")])
  .fun_testMetric(TRUE, "metric.binary", metric.binary, available.metrics.binary)
  .fun_testMetric(TRUE, "metric.filter", metric.filter, available.metrics.filter)
  
  
  ## select layers  -----------------------------------------------------------
  if (!is.null(metric.binary)) {
    selected.layers <- which(df.info$type == "bin" & df.info$metric == metric.binary)  
  } else if ( !is.null(metric.filter)) {
    selected.layers <- which(df.info$type == "filt" & df.info$metric == metric.filter)  
  } else {
    selected.layers <- which(df.info$type == "out")  
  }
  selected.layers
}



## FILL BIOMOD.models.out elements in bm.mod or bm.em objects -------------------------------------
## used in BIOMOD_Modeling, BIOMOD_EnsembleModeling

.fill_BIOMOD.models.out <- function(objName, objValue, mod.out, inMemory = FALSE, nameFolder = ".")
{
  # backup case uses existing @val slot
  if(is.null(objValue)){
    eval(parse(text = paste0("objValue <- mod.out@", objName, "@val")))
  }
  
  save(objValue, file = file.path(nameFolder, objName), compress = TRUE)
  if (inMemory) {
    eval(parse(text = paste0("mod.out@", objName, "@val <- objValue")))
  }
  eval(parse(text = paste0("mod.out@", objName, "@inMemory <- ", inMemory)))
  eval(parse(text = paste0("mod.out@", objName, "@link <- file.path(nameFolder, objName)")))
  return(mod.out)
}





## CATEGORICAL variables tools --------------------------------------------------------------------

## used in biomod2_classes_4, BIOMOD_Modeling, bm_RunModelsLoop

##' @name .get_categorical_names
##' @title Get categorical variable names
##' @description Internal function to get categorical variables name from a data.frame.
##' @param df data.frame to be checked
##' @return a vector with the name of categorical variables
##' @keywords internal

.get_categorical_names <- function(df){
  unlist(sapply(names(df), function(x) {
    if (is.factor(df[, x])) { return(x) } else { return(NULL) }
  }))
}


## used in bm_FindOptimStat, bm_PlotResponseCurves, BIOMOD_Projection, BIOMOD_EnsembleForecasting
.categorical_stack_to_terra <- function(myraster, expected_levels = NULL)
{
  myTerra <- rast(
    sapply(1:raster::nlayers(myraster), 
           function(thislayer){
             rast(myraster[[thislayer]])
           })
  )
  which.factor <- which(raster::is.factor(myraster))
  for (this.layer in which.factor) {
    this.levels <- raster::levels(myraster)[[this.layer]][[1]]
    ind.levels <- ifelse(ncol(this.levels) > 1, 2, 1)
    if (any(duplicated(this.levels[ , ind.levels]))) {
      stop("duplicated levels in environmental raster")
    }
    if (is.null(expected_levels)) {
      # no check to do, just formatting
      this.levels.df <- data.frame(ID = this.levels[,ind.levels],
                                   value = paste0(this.levels[,ind.levels]))
    } else {
      new.levels <- paste0(this.levels[,ind.levels])
      new.index <- this.levels[,1]
      this.layer.name <- names(myraster)[this.layer]
      fit.levels <- levels(expected_levels[,this.layer.name])
      if (!all(new.levels %in% fit.levels)) {
        cat("\n",
            "!! Levels for layer", colnames(expected_levels)[this.layer],
            " do not match.", new.levels[which(!new.levels %in% fit.levels)], "not found in fit data",
            "\n Fit data levels: ",paste0(fit.levels, collapse = " ; "),
            "\n Projection data levels: ",paste0(new.levels, collapse = " ; "))
        stop(paste0("Levels for ", colnames(expected_levels)[this.layer],
                    " do not match."))
      }
      levels.to.add <- which(!fit.levels %in% new.levels)
      if (length(levels.to.add) > 0) {
        max.index <- max(new.index)
        this.levels.df <- data.frame(ID = c(new.index,
                                            max.index +
                                              seq_along(levels.to.add)),
                                     value = c(new.levels, fit.levels[levels.to.add]))
      } else {
        this.levels.df <- data.frame(ID = new.index,
                                     value = new.levels)
      }
    }
    colnames(this.levels.df) <- c("ID", names(myTerra)[this.layer])
    try({myTerra <- categories(myTerra, layer = this.layer, this.levels.df)}, silent = TRUE)
  }
  return(myTerra)
}

## used in BIOMOD_Projection, BIOMOD_EnsembleForecasting
.check_env_levels <-  function(new.env, expected_levels)
{
  which.factor <- which(sapply(new.env, is.factor))
  if (inherits(new.env, 'SpatRaster')) {
    for (this.layer in which.factor) {
      this.layer.name <- names(new.env)[this.layer]
      this.levels <- cats(new.env)[[this.layer]]
      ind.levels <- ifelse(ncol(this.levels) > 1, 2, 1)
      new.levels <- paste0(this.levels[,ind.levels])
      if( any(duplicated(new.levels)) ) {
        stop("duplicated levels in `new.env`")
      }
      new.index <- this.levels[, 1]
      fit.levels <- levels(expected_levels[,this.layer.name])
      if (!all(new.levels %in% fit.levels)) {
        cat("\n",
            "!! Levels for layer", colnames(expected_levels)[this.layer],
            " do not match.", new.levels[which(!new.levels %in% fit.levels)], "not found in fit data",
            "\n Fit data levels: ",paste0(fit.levels, collapse = " ; "),
            "\n Projection data levels: ",paste0(new.levels, collapse = " ; "))
        stop(paste0("Levels for ", colnames(expected_levels)[this.layer],
                    " do not match."))
      }
      levels.to.add <- which(!fit.levels %in% new.levels)
      if (length(levels.to.add) > 0) {
        max.index <- max(new.index)
        this.levels.df <- data.frame(ID = c(new.index,
                                            max.index +
                                              seq_along(levels.to.add)),
                                     value = c(new.levels, fit.levels[levels.to.add]))
        
        colnames(this.levels.df) <- c("ID", names(new.env)[this.layer])
        new.env <- categories(new.env, layer = this.layer, 
                              this.levels.df)
      } 
    }
  } else {
    for (this.layer in which.factor) {
      this.layer.name <- names(new.env)[this.layer]
      this.levels <- levels(new.env[,this.layer.name])
      new.levels <- paste0(this.levels)
      fit.levels <- levels(expected_levels[,this.layer.name])
      if (!all(new.levels %in% fit.levels)) {
        cat("\n",
            "!! Levels for layer", colnames(expected_levels)[this.layer],
            " do not match.", new.levels[which(!new.levels %in% fit.levels)], "not found in fit data",
            "\n Fit data levels: ",paste0(fit.levels, collapse = " ; "),
            "\n Projection data levels: ",paste0(new.levels, collapse = " ; "))
        stop(paste0("Levels for ", colnames(expected_levels)[this.layer],
                    " do not match."))
      }
      levels.to.add <- which(!fit.levels %in% new.levels)
      if (length(levels.to.add) > 0) {
        new.env[ , this.layer.name] <- factor(new.env[ , this.layer.name], 
                                              levels = c(new.levels,
                                                         fit.levels[levels.to.add]))
      }
    }
  }
  new.env 
}


## PREDICTIONS TOOLS ------------------------------------------------------------------------------

## used in biomod2_classes_4

##' 
##' @importFrom terra rast classify subset
##' @importFrom biomod2 predict
##' 

.run_pred <- function(object, Prev = 0.5, dat, mod.name = NULL)
{
  if (is.finite(object$deviance) && 
      is.finite(object$null.deviance) && 
      object$deviance != object$null.deviance)
  {
    if (inherits(dat, 'SpatRaster')) {
      pred <- predict(object = dat, model = object, type = "response", wopt = list(names = mod.name))
    } else {
      pred <- predict(object, dat, type = "response")
    }
  }
  
  if (!exists('pred')) {
    if (inherits(dat, 'SpatRaster')) {
      pred <- subset(dat, 1)
      if (Prev < 0.5) {
        pred <- classify(x = pred, rcl = c(-Inf, Inf, 0))
      } else {
        pred <- classify(x = pred, rcl = c(-Inf, Inf, 1))
      }
    } else {
      if (Prev < 0.5) {
        pred <- rep(0, nrow(dat))
      } else {
        pred <- rep(1, nrow(dat))
      }
    }
  }
  
  return(pred)
}

## Predict on SpatRaster with XGBOOST
## used in biomod2_classes_4

.xgb_pred <- function(model, data, ...) {
  predict(model, newdata = as.matrix(data), ...)
}


## Get formal predictions for ensemble models
## used in biomod2_classes_5

##' 
##' @importFrom biomod2 predict
##' 
.get_formal_predictions <- function(object, newdata, on_0_1000, seedval)
{
  # make prediction of all models required
  formal_predictions <- sapply(object@model,
                               function(mod.name, dir_name, resp_name, modeling.id)
                               {
                                 ## check if model is loaded on memory
                                 if (is.character(mod.name)) {
                                   mod <- get(load(file.path(dir_name, resp_name, "models", modeling.id, mod.name)))
                                 }
                                 temp_workdir = NULL
                                 if (length(grep("MAXENT$", mod.name)) == 1) {
                                   temp_workdir = mod@model_output_dir
                                 }
                                 return(predict(mod, newdata = newdata, on_0_1000 = on_0_1000
                                                , temp_workdir = temp_workdir, seedval = seedval))
                               }, dir_name = object@dir_name, resp_name = object@resp_name, modeling.id = object@modeling.id)
  
  
  return(formal_predictions)
}


## ABUNDANCE tools --------------------------------------------------------------------------------

## used in biomod2_classes_1
.which.data.type <- function(resp.var)
{
  if (is.factor(resp.var)) {
    data.type <- ifelse(nlevels(resp.var) <= 2, "binary", "ordinal")
  } else {
    element <- sort(unique(resp.var))
    if (identical(as.numeric(element), c(0, 1)) | identical(as.numeric(element), 1)) {
      data.type <- "binary"
    } else {
      if (identical(as.numeric(as.integer(element)), element)) {
        data.type <- "count"
      } else {
        data.type <- ifelse(!any(resp.var > 1), "relative", "abundance")
      }
    }
  } 
  return(data.type)
}

## used in bm_RunModelsLoop
.numeric2factor <- function(pred, obs)
{
  nblevels <- length(levels(obs))
  pred <- round(pred)
  pred[pred > nblevels] <- nblevels
  pred[pred < 1] <- 1
  newlevels <- levels(obs)[sort(unique(pred))]
  pred <- factor(pred, labels = newlevels, ordered = TRUE)
  return(pred)
}

## CHOOSE the optimized cutoff between each ordinal variables
## used in bm_RunModelsLoop, bm_Tuning
# .threshold_ordinal <- function(obs, fit, metric.eval)
# {
#   nlevels <- length(levels(obs))
#   
#   #manage extremum
#   fit[fit < 1] <- 1
#   fit[fit > nlevels] <- nlevels
#   
#   #Define all the possibility
#   seq <- list(seq(0.2, 0.8, by = 0.1))
#   possibility <- rep(seq, nlevels-1)
#   names(possibility) <- paste0("Threshold_", 1:(nlevels-1))
#   grid <- expand.grid(possibility)
#   
#   test_seq <- foreach(i = 1:nrow(grid), .combine = rbind) %do%
#     {
#       fit_factor <- fit
#       for (j in 1:(nlevels - 1)) {
#         fit_factor[fit_factor <= j + grid[i, j] & fit_factor >= j ] <- j
#         fit_factor[fit_factor > j + grid[i, j] & fit_factor < j+1 ] <- j+1
#       }
#       stat <- bm_CalculateStatAbun(metric.eval, obs, fit_factor)
#       return(data.frame(grid[i, ], "stat" = stat))
#     }
#   
#   max <- test_seq[test_seq$stat == max(test_seq$stat), ]
#   limits <- colMeans(max) #several maximum !?!
#   
#   fit_factor <- fit
#   for (j in 1:(nlevels - 1)) {
#     fit_factor[fit_factor <= j + limits[j] & fit_factor >= j ] <- j
#     fit_factor[fit_factor > j + limits[j] & fit_factor < j + 1 ] <- j + 1
#   }
#   
#   levels <- 1:length(levels(obs))
#   names(levels) <- levels(obs)
#   fit_factor <- factor(fit_factor, levels = levels, ordered = TRUE)
#   
#   return(list(limits = limits, fit_factor = fit_factor))
# }
