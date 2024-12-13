# MS_EnsembleModeling ---------------------------------------------------------
##' @name MS_EnsembleModeling
##' @author Helene Blancheteau
##' 
##' @title Create and evaluate an ensemble set of models and predictions
##' 
##' @description This function allows to combine a range of models built with the 
##' \code{\link{MS_Modeling}} function in one (or several) ensemble model. Modeling 
##' uncertainty can be assessed as well as variables importance, ensemble predictions can be 
##' evaluated against original data, and created ensemble models can be projected over new 
##' conditions (see Details).
##' 
##' 
##' @param bm.mod a \code{\link{MS.models.out}} object returned by the 
##' \code{\link{MS_Modeling}} function
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' @param em.by a \code{character} corresponding to the way kept models will be combined to build 
##' the ensemble models, must be among \code{all}, \code{algo}, \code{PA}, \code{PA+algo}, 
##' \code{PA+run}
##' @param em.algo a \code{vector} corresponding to the ensemble models that will be computed, 
##' must be among \code{'EMmean'}, \code{'EMmedian'}, \code{'EMcv'}, \code{'EMci'}, 
##' \code{'EMca'}, \code{'EMwmean'}
##' @param metric.select a \code{vector} containing evaluation metric names to be used together with 
##' \code{metric.select.thresh} to exclude single models based on their evaluation scores 
##' (for ensemble methods like probability weighted mean or committee averaging). Must be among  
##' \code{all} (same evaluation metrics than those of \code{bm.mod}), \code{user.defined} 
##' (and defined through \code{metric.select.table}) or \code{POD}, \code{FAR}, \code{POFD}, 
##' \code{SR}, \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, 
##' \code{ORSS}, \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param metric.select.thresh (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to the minimum scores (one for each 
##' \code{metric.select}) below which single models will be excluded from the ensemble model 
##' building
##' @param metric.select.table (\emph{optional, default} \code{NULL}) \cr If 
##' \code{metric.select = 'user.defined'}, a \code{data.frame} containing evaluation scores 
##' calculated for each single models and that will be compared to \code{metric.select.thresh} 
##' values to exclude some of them from the ensemble model building, with \code{metric.select} 
##' rownames, and \code{models.chosen} colnames
##' @param metric.select.dataset (\emph{optional, default} \code{'validation'} 
##' \emph{if possible}). A character determining which dataset should be used to filter and/or 
##' weigh the ensemble models should be among 'evaluation', 'validation' or 'calibration'.
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among  \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, 
##' \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, 
##' \code{BOYCE}, \code{MPA}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param EMci.alpha (\emph{optional, default} \code{0.05}) \cr 
##' A \code{numeric} value corresponding to the significance level to estimate confidence interval
##' @param EMwmean.decay (\emph{optional, default} \code{proportional}) \cr A
##'   value defining the relative importance of the weights (if \code{'EMwmean'}
##'   was given to argument \code{em.algo}). A high value will strongly
##'   discriminate \emph{good} models from the \emph{bad} ones (see Details),
##'   while \code{proportional} will attribute weights proportionally to the
##'   models evaluation scores
##' 
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models predictions and the ensemble models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' @param do.progress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the progress bar is to be rendered or not
##'
##'@return
##' A \code{\link{MS.ensemble.models.out}} object acting as a proxi for the created \code{BIOMOD.ensemble.models.out}.
##' 
##' 
##' 
##' @keywords models ensemble weights
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{bm_ModelingOptions}}, 
##' \code{\link{bm_CrossValidation}}, \code{\link{bm_VariablesImportance}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleForecasting}},
##' \code{\link{bm_PlotEvalMean}}, \code{\link{bm_PlotEvalBoxplot}}, 
##' \code{\link{bm_PlotVarImpBoxplot}}, \code{\link{bm_PlotResponseCurves}}
##' @family Main functions
##' 
##' @importFrom foreach foreach %do%
##' @importFrom biomod2 BIOMOD_Modeling
##' 
##' @export
##' 
##' 
###################################################################################################

MS_EnsembleModeling <- function(ms.mod,
                                models.chosen = 'all',
                                em.by = 'PA+run',
                                em.algo,
                                metric.select = 'all',
                                metric.select.thresh = NULL,
                                metric.select.table = NULL,
                                metric.select.dataset = NULL,
                                metric.eval = c('KAPPA', 'TSS', 'ROC'),
                                var.import = 0,
                                EMci.alpha = 0.05,
                                EMwmean.decay = 'proportional',
                                nb.cpu = 1,
                                seed.val = NULL)
{
  .bm_cat("Build Ensemble Models")

  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .MS_EnsembleModeling.check.args(
    ms.mod = ms.mod, 
    models.chosen = models.chosen,
    em.by = em.by,
    em.algo = em.algo,
    metric.select = metric.select,
    metric.select.thresh = metric.select.thresh,
    metric.select.table = metric.select.table,
    metric.select.dataset = metric.select.dataset,
    metric.eval = metric.eval,
    EMci.alpha = EMci.alpha,
    EMwmean.decay = EMwmean.decay
  )
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)

  
  MSEM<- new(
    'MS.ensemble.models.out',
    ms.project = ms.mod@ms.project,
    modeling.id = ms.mod@modeling.id,
    dir.name = ms.mod@dir.name,
    sp.name = ms.mod@sp.name,
    data.type = ms.mod@data.type,
    expl.var.names = ms.mod@expl.var.names,
    em.computed = list(),
    em.failed = list(),
    em.models_kept = list()
  )
  
  cat("\n")
  workflow <- foreach(sp = ms.mod@sp.name) %do% {
    cat("\n\t Ensemble Modeling of", sp)
    # 1. Récupération ms.mod 
    bm.mod <- get(load(file.path(ms.mod@dir.name, sp, paste0(sp, ".", ms.mod@modeling.id,".models.out"))))
    
    em.by.sp <- em.by[[sp]]
    models.chosen.sp <- models.chosen[[sp]]
    metric.select.sp <- metric.select[[sp]]
    metric.select.thresh.sp <- metric.select.thresh[[sp]]
    metric.select.table.sp <- metric.select.table[[sp]]
    metric.select.dataset.sp <- metric.select.dataset[[sp]]

    
    # 2. Run MS_EnsembleModeling
    output <- capture.output(em_models <- BIOMOD_EnsembleModeling(bm.mod,
                                                                  models.chosen = models.chosen.sp,
                                                                  em.by = em.by.sp,
                                                                  em.algo = em.algo,
                                                                  metric.select = metric.select.sp,
                                                                  metric.select.thresh = metric.select.thresh.sp,
                                                                  metric.select.table = metric.select.table.sp,
                                                                  metric.select.dataset = metric.select.dataset.sp,
                                                                  metric.eval = metric.eval,
                                                                  var.import = var.import,
                                                                  EMci.alpha = EMci.alpha,
                                                                  EMwmean.decay = EMwmean.decay,
                                                                  nb.cpu = 1,
                                                                  seed.val = NULL,
                                                                  do.progress = FALSE))
    
    # 3.Stockage
    MSEM@em.computed[[sp]] <- em_models@em.computed
    MSEM@em.failed[[sp]] <- em_models@em.failed
    MSEM@em.models_kept[[sp]] <- em_models@em.models_kept
  }
  
  ## SAVE MODEL OBJECT ON HARD DRIVE ----------------------------
  name.OUT = paste0(ms.mod@ms.project, '.', ms.mod@modeling.id, '.MS.ensemble.models.out')
  MSEM@link <- file.path(ms.mod@dir.name, name.OUT)
  assign(x = name.OUT, value = MSEM)
  save(list = name.OUT, file = MSEM@link)
  
  .bm_cat("Done")
  return(MSEM)
}


# ---------------------------------------------------------------------------- #

# Argument check function  -----------------------------------------------------

.MS_EnsembleModeling.check.args <- function(ms.mod,
                                            models.chosen,
                                            em.by,
                                            em.algo,
                                            metric.select,
                                            metric.select.thresh,
                                            metric.select.table,
                                            metric.select.dataset,
                                            metric.eval,
                                            EMci.alpha,
                                            EMwmean.decay) { 
  ## 1. Check bm.mod ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "ms.mod", ms.mod, "MS.models.out")
  
  ## 2. Check models.chosen ---------------------------------------------------
  if ( is.null(models.chosen) | (length(models.chosen) == 1 && models.chosen[1] == 'all')) {
    models.chosen <- as.list(rep("all", length(ms.mod@sp.name)))
    names(models.chosen) <- ms.mod@sp.name
  } else {
    .fun_testIfInherits(TRUE, "models.chosen", models.chosen, "list")
    .fun_testIfIn(TRUE, "names(models.chosen)", names(models.chosen), ms.mod@sp.name)
  }
  
  # 3. check argument em.algo ----------------------------------------------
  em.avail.check <- c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean')
  em.avail <- c('EMmean', 'EMcv', 'EMciInf', 'EMciSup', 'EMmedian', 'EMca', 'EMwmean')
  if (missing(em.algo)) {
    em.algo <- 'EMmean'
    cat("\n! setting em.algo to its default value c('EMmean')")
  } else {
    .fun_testIfIn(TRUE, "em.algo", em.algo, em.avail.check)
    em.algo <- unique(em.algo)
    # testCI <- grepl(pattern = "EMci", x = em.algo)
    # if(any(testCI)){
    #   em.algo <- em.algo[-which(testCI)]
    #   em.algo <- c(em.algo, 'EMciInf', 'EMciSup')
    # }
    if (ms.mod@data.type != "binary" & 'EMca' %in% em.algo){
      cat ("\n\t EMca is not available with",ms.mod@data.type, "data")
      em.algo <- em.algo[-which(em.algo == "EMca")]
    }
  }

  
  ## 4. Check metric.select ---------------------------------------------------
  metric.select.user = FALSE
  if (!is.list(metric.select)){
    initial.metric.select <- metric.select
    if (!is.null(metric.select)) {
      if (!is.character(metric.select)) {
        stop("metric.select must be a character vector or NULL")
      }
      if ('user.defined' %in% metric.select) {
        metric.select.user = TRUE
        if (!is.null(metric.select.table)) {
          .fun_testIfIn(TRUE, "models.chosen", models.chosen, colnames(metric.select.table))
          metric.select.table <- metric.select.table[, models.chosen, drop = FALSE]
          metric.select <- rownames(metric.select.table)
        } else {
          stop("metric.select.table must be a data.frame or NULL")
        }
      } 
    initial.metric.select <- metric.select
    metric.select <- as.list(rep(list(metric.select), length(ms.mod@sp.name)))
    names(metric.select) <- ms.mod@sp.name
    }
  } else {
    initial.metric.select <- metric.select[[1]]
    #.fun_testIfInherits(TRUE, "models.chosen", models.chosen, "list")
    .fun_testIfIn(TRUE, "names(metric.select)", names(metric.select), ms.mod@sp.name)
  }
  
  ## 5. metric.select.dataset -------------------------------------------------
  has.validation.data <- any(!is.na((get_evaluations(ms.mod, sp = ms.mod@sp.name[1]))$validation))
  has.evaluation.data <- ms.mod@has.evaluation.data
  
  metric.select.dataset.available <- c("calibration")
  if (has.validation.data) {
    metric.select.dataset.available <- 
      append(metric.select.dataset.available, "validation")
  }
  if (has.evaluation.data) {
    metric.select.dataset.available <- 
      append(metric.select.dataset.available, "evaluation")
  }
  
  if(!is.list(metric.select.dataset)){
    if (is.null(metric.select.dataset)) {
      if (has.validation.data) {
        metric.select.dataset <- "validation"
        cat("\n  ! Ensemble Models will be filtered and/or weighted using validation dataset (if possible). Please use `metric.select.dataset` for alternative options.")
      } else {
        metric.select.dataset <- "calibration"
        cat("\n  ! Ensemble Models will be filtered and/or weighted using calibration dataset. Please use `metric.select.dataset` for alternative options.")
      }
    } else {
      .fun_testIfIn(TRUE, "metric.select.dataset",
                    metric.select.dataset, metric.select.dataset.available)
    }
    metric.select.dataset <- as.list(rep(metric.select.dataset, length(ms.mod@sp.name)))
    names(metric.select.dataset) <- ms.mod@sp.name
  } else {
    .fun_testIfIn(TRUE, "names(metric.select.dataset)", names(metric.select.dataset), ms.mod@sp.name)
  }
  

  ## 6. Check metric.select.thresh --------------------------------------------
  print(metric.select.thresh)
  if (!is.null(initial.metric.select)) {
    if(!is.list(metric.select.thresh)){
      if (!is.null(metric.select.thresh)) {
        if (!is.numeric(metric.select.thresh)) {
          stop("metric.select.thresh must be NULL or a numeric vector")
        }
        if (length(initial.metric.select) != length(metric.select.thresh)) {
          stop("you must specify as many metric.select.thresh as metric.select (if you specify some)")
        }
        cat("\n   > Evaluation & Weighting methods summary :\n")
        if (any(c("RMSE", "MSE", "MAE", "Max_error") %in% metric.select)){
          metric.select.over <- metric.select[-which(metric.select %in% c("RMSE", "MSE", "MAE", "Max_error"))]
          metric.select.thresh.over <- metric.select.thresh[-which(metric.select %in% c("RMSE", "MSE", "MAE", "Max_error"))]
          cat(paste(metric.select.over, metric.select.thresh.over, sep = " over ", collapse = "\n      ")
              , fill = TRUE, labels = "     ")
          
          metric.select.under <- metric.select[which(metric.select %in% c("RMSE", "MSE", "MAE", "Max_error"))]
          metric.select.thresh.under <- metric.select.thresh[which(metric.select %in% c("RMSE", "MSE", "MAE", "Max_error"))]
          cat(paste(metric.select.under, metric.select.thresh.under, sep = " under the best + ", collapse = "\n      ")
              , fill = TRUE, labels = "     ")
        } else {
          cat(paste(metric.select, metric.select.thresh, sep = " over ", collapse = "\n      ")
              , fill = TRUE, labels = "     ")
        }
      } else {
        cat("\n   ! No metric.select.thresh -> All models will be kept for Ensemble Modeling")
        #metric.select.thresh <- ifelse(metric.select %in% c("RMSE", "MSE", "MAE", "Max_error"), 100, 0)
      }
      metric.select.thresh <- as.list(rep(list(metric.select.thresh), length(ms.mod@sp.name)))
      names(metric.select.thresh) <- ms.mod@sp.name
    } else {
      .fun_testIfIn(TRUE, "names(metric.select.thresh)", names(metric.select.thresh), ms.mod@sp.name)
    }
  }
  
  
  
  ## 7. Check metric.eval -----------------------------------------------------
  metric.eval <- unique(metric.eval)
  if (ms.mod@data.type == "binary"){
    avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                              , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC'
                              , 'BOYCE', 'MPA')
  } else if (ms.mod@data.type == "ordinal"){
    avail.eval.meth.list <- c("Accuracy", "Recall", "Precision", "F1")
  } else {
    avail.eval.meth.list <- c('RMSE','MSE',"MAE","Rsquared","Rsquared_aj","Max_error")
  }
  .fun_testIfIn(TRUE, paste0("metric.eval with ", ms.mod@data.type, " data type"), metric.eval, avail.eval.meth.list)
  
  
  ## 8. Check selected EM algo ------------------------------------------------
  if (is.null(metric.select) && 
      any(c("committee.averaging", "prob.mean.weight") %in% em.algo)) {
    stop("You must choose metric.select if you want to compute Committee Averaging or Probability Weighted Mean algorithms")
  }

  
  ## 8.1 Check alpha for Confident interval
  if ("EMci" %in% em.algo) {
    .fun_testIfPosNum(TRUE, "EMci.alpha", EMci.alpha)
    if (EMci.alpha <= 0 | EMci.alpha >= 0.5) {
      stop("EMci.alpha must be a numeric between 0 and 0.5")
    }
  }
  # prob.mean.weight.decay
  ## 8.2 Check decay for wmean
  if ("EMwmean" %in% em.algo) {
    if ((!is.numeric(EMwmean.decay) &&
         !is.character(EMwmean.decay) &&
         !is.function(EMwmean.decay)) ||
        (is.numeric(EMwmean.decay) && EMwmean.decay < 0) ||
        (is.character(EMwmean.decay) && EMwmean.decay != 'proportional')) {
      stop("'EMwmean.decay' should be either 'proportional', a numeric value > 0 or a function")
    }
  }
  
  ## 9. Check em.by -----------------------------------------------------------
  
  if (!is.list(em.by)){
    if(length(em.by) != 1){
      stop("\nem.by should be of length 1")
    }
    em.by.avail.old <- c("PA_dataset"       = "PA",
                         "PA_dataset+repet" = "PA+run",
                         "PA_dataset+algo"  = "PA+algo")
    em.by.avail <- c('PA', 'algo', 'all', 'PA+run', 'PA+algo')
    
    if(missing(em.by)){
      em.by <- "all"
      cat("\n! `em.by` automatically set to 'all'")
    }
    
    .fun_testIfIn(TRUE, "em.by", em.by, em.by.avail)
    em.by <- as.list(rep(em.by, length(ms.mod@sp.name)))
    names(em.by) <- ms.mod@sp.name
  } else {
    .fun_testIfIn(TRUE, "names(em.by)", names(em.by), ms.mod@sp.name)
  }
  
  
  
  # # check that repetition are note merged with full models
  # if(any(grepl(pattern = "RUN",  x = models.chosen)) &&
  #    any(grepl(pattern = "allRun", x = models.chosen)) &&
  #    em.by != 'PA+run') {
  #   cat("\n!!! Removed models using the Full dataset as ensemble models cannot merge repetition dataset (RUN1, RUN2, ...) with Full dataset unless em.by = 'PA+run'.")
  #   models.chosen <- models.chosen[!grepl(pattern = "allRun", x = models.chosen)]
  # }
  
  ## 10. Check that ensemble model have > 1 model to run -------------
  ## make a list of models names that will be combined together according to em.by argument
  # em.mod.assemb <- .get_models_assembling(models.chosen, em.by)
  # ### Check that all EM have > 1 model selected ----------------------------
  # out.check <- foreach(assemb = names(em.mod.assemb), .combine = 'rbind') %do% {
  #   out <- .get_kept_models(
  #     bm.mod, em.mod.assemb[[assemb]], 
  #     metric.select, metric.select.thresh,
  #     metric.select.user, metric.select.table,
  #     metric.select.dataset
  #   )$models.kept
  #   data.frame(models.kept = sapply(out, length),
  #              metric.select = names(out),
  #              assemb = assemb,
  #              row.names = NULL)
  # }
  # 
  # for (thismetric in metric.select)  {
  #   out.check.sub <- out.check[which(out.check$metric.select == thismetric),]
  #   assemb.1 <- out.check.sub[which(out.check.sub$models.kept == 1), "assemb"]
  #   assemb.0 <- out.check.sub[which(out.check.sub$models.kept == 0), "assemb"]
  #   
  #   if(length(assemb.0) > 0 || length(assemb.1) > 0){
  #     cat("\n")
  #     if(length(assemb.0) > 0){
  #       cat("\n     !! Ensemble Model", assemb.0, "selected by", thismetric, "have no model selected and will fail.")
  #     }   
  #     if(length(assemb.1) > 0){
  #       cat("\n     !! Ensemble Model", assemb.1, "selected by", thismetric, "have only one single model selected.")
  #     }
  #     cat("\n     !! Please make sure this is intended or review your selection metrics and threshold.")
  #   }
  # }
  # 
  
  
  
  ## Return Args ------------------------------------------------
  return(list(ms.mod = ms.mod,
              models.chosen = models.chosen,
              em.algo = em.algo,
              metric.select = metric.select,
              metric.select.thresh = metric.select.thresh,
              metric.select.user = metric.select.user,
              metric.select.table = metric.select.table,
              metric.select.dataset = metric.select.dataset,
              metric.eval = metric.eval,
              EMci.alpha = EMci.alpha,
              EMwmean.decay = EMwmean.decay,
              em.by = em.by))
              #em.mod.assemb = em.mod.assemb))
  
}

