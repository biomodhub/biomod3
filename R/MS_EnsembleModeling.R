# MS_EnsembleModeling ---------------------------------------------------------
##' @name MS_EnsembleModeling
##' @author Helene Blancheteau
##' 
##' @title Run a range of species distribution models
##' 
##' @description This function allows to calibrate and evaluate a range of modeling techniques 
##' for a given species distribution. The dataset can be split up in calibration/validation parts,
##' and the predictive power of the different models can be estimated using a range of evaluation 
##' metrics (see Details).
##' 
##' 
##' @param ms.mod a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' @param models.pa (\emph{optional, default} \code{NULL}) \cr 
##' A \code{list} containing for each model a \code{vector} defining which pseudo-absence datasets 
##' are to be used, must be among \code{colnames(ms.mod@PA.table)}
##' 
##' @param CV.strategy a \code{character} corresponding to the cross-validation selection strategy, 
##' must be among \code{random}, \code{kfold}, \code{block}, \code{strat}, \code{env} or 
##' \code{user.defined}
##' @param CV.nb.rep (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'random'} or \code{strategy = 'kfold'}, an \code{integer} corresponding 
##' to the number of sets (repetitions) of cross-validation points that will be drawn
##' @param CV.perc (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'random'}, a \code{numeric} between \code{0} and \code{1} defining the 
##' percentage of data that will be kept for calibration
##' @param CV.k (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'kfold'} or \code{strategy = 'strat'} or \code{strategy = 'env'}, an 
##' \code{integer} corresponding to the number of partitions 
##' @param CV.balance (\emph{optional, default} \code{'presences'}) \cr
##' If \code{strategy = 'strat'} or \code{strategy = 'env'}, a \code{character} corresponding 
##' to how data will be balanced between partitions, must be either \code{presences} or
##' \code{absences} 
##' @param CV.env.var (\emph{optional}) \cr If \code{strategy = 'env'}, a \code{character} 
##' corresponding to the environmental variables used to build the partition. \code{k} partitions 
##' will be built for each environmental variables. By default the function uses all 
##' environmental variables available.
##' @param CV.strat (\emph{optional, default} \code{'both'}) \cr
##' If \code{strategy = 'env'}, a \code{character} corresponding to how data will partitioned 
##' along gradient, must be among \code{x}, \code{y}, \code{both}
##' @param CV.user.table (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'user.defined'}, a \code{matrix} or \code{data.frame} defining for each 
##' repetition (in columns) which observation lines should be used for models calibration 
##' (\code{TRUE}) and validation (\code{FALSE})
##' @param CV.do.full.models (\emph{optional, default} \code{TRUE}) \cr  
##' A \code{logical} value defining whether models should be also calibrated and validated over 
##' the whole dataset (and pseudo-absence datasets) or not
##' 
##' @param OPT.strategy a \code{character} corresponding to the method to select models' 
##' parameters values, must be either \code{default}, \code{bigboss}, \code{user.defined}, 
##' \code{tuned}
##' @param OPT.user.val (\emph{optional, default} \code{NULL}) \cr
##' A \code{list} containing parameters values for some (all) models
##' @param OPT.user.base (\emph{optional, default} \code{bigboss}) \cr A character, 
##' \code{default} or \code{bigboss} used when \code{OPT.strategy = 'user.defined'}. 
##' It sets the bases of parameters to be modified by user defined values.
##' @param OPT.user (\emph{optional, default} \code{TRUE}) \cr  
##' A \code{\link{BIOMOD.models.options}} object returned by the \code{\link{bm_ModelingOptions}} 
##' function
##' 
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
##' \code{Max_error}, \code{Accuracy}, \code{"Recall"}, \code{"Precision"}, \code{"F1"}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param scale.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions should be scaled with a 
##' binomial GLM or not
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' @param do.progress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the progress bar is to be rendered or not
##' 
##' 
##' @return
##' 
##' A \code{\link{BIOMOD.models.out}} object containing models outputs, or links to saved outputs. \cr
##' Models outputs are stored out of \R (for memory storage reasons) in 2 different folders 
##' created in the current working directory :
##' \enumerate{
##'   \item a \emph{models} folder, named after the \code{resp.name} argument of 
##'   \code{\link{BIOMOD_FormatingData}}, and containing all calibrated models for each 
##'   repetition and pseudo-absence run
##'   \item a \emph{hidden} folder, named \code{.BIOMOD_DATA}, and containing outputs related 
##'   files (original dataset, calibration lines, pseudo-absences selected, predictions, 
##'   variables importance, evaluation values...), that can be retrieved with 
##'   \href{https://biomodhub.github.io/biomod2/reference/getters.out.html}{\code{get_[...]}} 
##'   or \code{\link{load}} functions, and used by other \pkg{biomod2} functions, like 
##'   \code{\link{BIOMOD_Projection}} or \code{\link{BIOMOD_EnsembleModeling}}
##' }
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
  .bm_cat("Build Single Models")
  
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
    
    models.chosen.sp <- models.chosen[[sp]]
    # 2. Run MS_EnsembleModeling
    output <- capture.output(em_models <- BIOMOD_EnsembleModeling(bm.mod,
                                                                  models.chosen = models.chosen.sp,
                                                                  em.by = em.by,
                                                                  em.algo =em.algo,
                                                                  metric.select = metric.select,
                                                                  metric.select.thresh = metric.select.thresh,
                                                                  metric.select.table = metric.select.table,
                                                                  metric.select.dataset = metric.select.dataset,
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
    testCI <- grepl(pattern = "EMci", x = em.algo)
    if(any(testCI)){
      em.algo <- em.algo[-which(testCI)]
      em.algo <- c(em.algo, 'EMciInf', 'EMciSup')
    }
    if (ms.mod@data.type != "binary" & 'EMca' %in% em.algo){
      cat ("\n\t EMca is not available with",ms.mod@data.type, "data")
      em.algo <- em.algo[-which(em.algo == "EMca")]
    }
  }

  
  ## 4. Check metric.select ---------------------------------------------------
  metric.select.user = FALSE
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
    } else {
      if ('all' %in% metric.select) {
        metric.select <- unique(get_evaluations(ms.mod, sp = ms.mod@sp.name[1])$metric.eval)
      }
      .fun_testIfIn(TRUE, "metric.select", metric.select, unique(get_evaluations(ms.mod, sp = ms.mod@sp.name[1])$metric.eval))
      ## Remove MPA from metric.select
      if ('MPA' %in% metric.select) {
        metric.select.thresh <- metric.select.thresh[which(metric.select != 'MPA')]
        metric.select <- metric.select[which(metric.select != 'MPA')]
      }
    }
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
  
  ## 6. Check metric.select.thresh --------------------------------------------
  if (!is.null(metric.select)) {
    if (!is.null(metric.select.thresh)) {
      if (!is.numeric(metric.select.thresh)) {
        stop("metric.select.thresh must be NULL or a numeric vector")
      }
      if (length(metric.select) != length(metric.select.thresh)) {
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
      metric.select.thresh <- rep(0, length(metric.select))
    }
  } else {
    metric.select <- 'none'
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

