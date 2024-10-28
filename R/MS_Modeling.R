# MS_Modeling ---------------------------------------------------------
##' @name MS_Modeling
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
##' @param ms.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' @param models.pa (\emph{optional, default} \code{NULL}) \cr 
##' A \code{list} containing for each model a \code{vector} defining which pseudo-absence datasets 
##' are to be used, must be among \code{colnames(ms.format@PA.table)}
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

MS_Modeling <- function(ms.format,
                            modeling.id = as.character(format(Sys.time(), "%s")),
                            models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS'
                                       , 'MAXNET', 'RF', 'RFd', 'SRE', 'XGBOOST'),
                            params.CV = NULL,
                            params.OPT = NULL,
                            weights = NULL,
                            prevalence = NULL,
                            metric.eval = c('KAPPA', 'TSS', 'ROC'),
                            var.import = 0,
                            scale.models = FALSE,
                            nb.cpu = 1,
                            seed.val = NULL)
{
  .bm_cat("Build Single Models")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .MS_Modeling.check.args(
    ms.format = ms.format, 
    modeling.id = modeling.id, 
    models = models, 
    weights = weights, 
    prevalence = prevalence, 
    metric.eval = metric.eval, 
    var.import = var.import, 
    scale.models = scale.models,
    nb.cpu = nb.cpu,
    seed.val = seed.val
  )
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  MSmodels<- new(
    'MS.models.out',
    ms.project = ms.format@ms.project,
    modeling.id = modeling.id,
    dir.name = file.path(ms.format@dir.name, ms.format@ms.project),
    sp.name = ms.format@sp.name,
    data.type = ms.format@data.type,
    expl.var.names = names(ms.format@data.env.var),
    scale.models = scale.models
  )
  
  nameFolder <- file.path(ms.format@dir.name, ms.format@ms.project, ".BIOMOD_DATA", "single.formated.data")
  
  dir.name <- file.path(ms.format@dir.name, ms.format@ms.project)
  
  workflow <- foreach(sp = ms.format@sp.name) %do% {
    cat("\n\t Modeling of", sp)
    # 1. Récupération ms.format 
    sfd <- get(load(file.path(nameFolder, paste0(sp,".sfd"))))
    
    parameter.CV <- params.CV[[sp]]
    parameter.OPT <- params.OPT[[sp]]
    
    # 2. Run MS_Modeling
    output <- capture.output(sfd_models <- BIOMOD_Modeling(sfd,
                                                       modeling.id = modeling.id,
                                                       models = models,
                                                       CV.strategy = parameter.CV$CV.strategy,
                                                       CV.nb.rep = parameter.CV$CV.nb.rep,
                                                       CV.perc = parameter.CV$CV.perc,
                                                       CV.k = parameter.CV$CV.k,
                                                       CV.balance = parameter.CV$CV.balance,
                                                       CV.env.var = parameter.CV$CV.env.var,
                                                       CV.strat = parameter.CV$CV.strat,
                                                       CV.user.table = parameter.CV$CV.user.table,
                                                       CV.do.full.models = TRUE,
                                                       OPT.strategy = parameter.OPT$OPT.strategy,
                                                       OPT.user.val = parameter.OPT$OPT.user.val,
                                                       OPT.user.base = parameter.OPT$OPT.user.base,
                                                       OPT.user = parameter.OPT$OPT.user,
                                                       weights = weights,
                                                       prevalence = prevalence,
                                                       metric.eval = metric.eval,
                                                       var.import = var.import,
                                                       scale.models = scale.models,
                                                       nb.cpu = 1,
                                                       seed.val = NULL))
    
    # 3.Stockage
    summary.models.computed <- data.frame("species" = sp, 
                                          "models.computed" = length(sfd_models@models.computed),
                                          "models.failed" = length(sfd_models@models.failed))
    MSmodels@summary.models.computed <- rbind(MSmodels@summary.models.computed, summary.models.computed)
  }
  
  .bm_cat("Done")
  return(MSmodels)
}


# ---------------------------------------------------------------------------- #

.MS_Modeling.prepare.workdir <- function(dir.name, sp.name, modeling.id)
{
  cat("\nCreating suitable Workdir...\n")
  dir.create(file.path(dir.name, sp.name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, ".BIOMOD_DATA", modeling.id), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, "models", modeling.id), showWarnings = FALSE, recursive = TRUE)
}


# ---------------------------------------------------------------------------- #

.MS_Modeling.check.args <- function(ms.format, modeling.id, models
                                        , weights, prevalence, metric.eval, var.import
                                        , scale.models, nb.cpu, seed.val)
{
  ## 0. Check ms.format and models arguments ----------------------------------
  cat('\n\nChecking Models arguments...\n')
  
  if (!is.character(modeling.id) || length(modeling.id) > 1) { stop("modeling.id must be a 'character' of length 1") }
  
  .fun_testIfInherits(TRUE, "ms.format", ms.format, "MS.formated.data")
  if (!is.character(models)) { stop("models must be a 'character' vector") }
  
  models <- unique(models)
  models.switch.off <- NULL
  
  ## check if model is supported
  if (ms.format@data.type == "binary"){
    avail.models.list <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF','RFd', 'SRE', 'XGBOOST')
  } else if (ms.format@data.type == "ordinal") {
    avail.models.list <- c('CTA', 'FDA', 'GAM', 'GLM', 'MARS', 'RF', 'XGBOOST')
  } else {
    avail.models.list <- c('CTA', 'GAM', 'GBM', 'GLM', 'MARS', 'RF', 'XGBOOST')
  }
  .fun_testIfIn(TRUE, paste0("models with ", ms.format@data.type, " data type"), models, avail.models.list)
  
  
  ## Specific case of one variable with GBM / MAXNET
  if ('GBM' %in% models && ncol(ms.format@data.env.var) == 1) {
    warning('GBM might have issues when only one variable is used. Please be sure to install the following version : devtools::install_github("rpatin/gbm")')
  }
  if ('MAXNET' %in% models && ncol(ms.format@data.env.var) == 1) {
    warning('MAXNET might have issues when only one variable is used. Please be sure to install the following version : devtools::install_github("mrmaxent/maxnet")')
  }
  
  ## 1.1 Remove models not supporting categorical variables --------------------
  categorical_var <- .get_categorical_names(ms.format@data.env.var)
  
  if (length(categorical_var) > 0) {
    models.fact.unsupport <- c("SRE", "XGBOOST")
    models.switch.off <- c(models.switch.off, intersect(models, models.fact.unsupport))
    if (length(models.switch.off) > 0) {
      models <- setdiff(models, models.switch.off)
      cat(paste0("\n\t! ", paste(models.switch.off, collapse = ",")," were switched off because of categorical variables !"))
    }
  }
  
  
  ## 5. Check prevalence arguments --------------------------------------------
  if (!is.null(prevalence)) {
    .fun_testIf01(TRUE, "prevalence", prevalence)
  } else {
    prevalence = 0.5
  }
  
  ## 6. Check weights arguments -----------------------------------------------
  # if (is.null(weights)) {
  #   if (!is.null(prevalence) & ms.format@data.type != "ordinal") {
  #     cat("\n\t> Automatic weights creation to rise a", prevalence, "prevalence")
  #     data.sp <- as.numeric(ms.format@data.species)
  #     if (inherits(ms.format, "BIOMOD.formated.data.PA")) {
  #       weights.pa <- foreach(pa = 1:ncol(ms.format@PA.table), .combine = "cbind") %do%
  #         {
  #           ind.PA <- which(ms.format@PA.table[, pa] == TRUE)
  #           data.sp_pa <- data.sp[ind.PA]
  #           data.sp_pa[which(is.na(data.sp_pa))] <- 0
  #           weights <- .automatic_weights_creation(data.sp_pa, prev = prevalence)
  #           
  #           wei <- rep(NA, length(data.sp))
  #           wei[ind.PA] <- weights
  #           return(matrix(wei, ncol = 1))
  #         }
  #       weights.pa <- cbind(weights.pa, rep(1, nrow(weights.pa)))
  #       colnames(weights.pa) <- c(colnames(ms.format@PA.table), "allData")
  #       weights <- weights.pa
  #     } else {
  #       weights <- .automatic_weights_creation(data.sp, prev = prevalence)
  #       weights <- matrix(weights, nrow = length(weights), ncol = 1)
  #       colnames(weights) <- "allData"
  #     }
  #   } else { ## NEVER OCCURRING NO ?? --> now happen with the abundance
  #     cat("\n\t> No weights : all observations will have the same weight\n")
  #     #weights <- rep(1, length(ms.format@data.species))
  #   }
  # } else {
  #   if (!is.numeric(weights)) { stop("weights must be a numeric vector") }
  #   if (length(weights) != length(ms.format@data.species)) {
  #     stop("The number of 'Weight' does not match with the input calibration data. Simulation cannot proceed.")
  #   }
  #   if (inherits(ms.format, "BIOMOD.formated.data.PA")) {
  #     weights.pa <- foreach(pa = 1:ncol(ms.format@PA.table), .combine = "cbind") %do%
  #       {
  #         wei <- weights
  #         wei[which(ms.format@PA.table[, pa] == FALSE | is.na(ms.format@PA.table[, pa]))] <- NA
  #         return(matrix(wei, ncol = 1))
  #       }
  #     weights.pa <- cbind(weights.pa, rep(1, nrow(weights.pa)))
  #     colnames(weights.pa) <- c(colnames(ms.format@PA.table), "allData")
  #     weights <- weights.pa
  #   } else {
  #     weights <- matrix(weights, nrow = length(weights), ncol = 1)
  #     colnames(weights) <- "allData"
  #   }
  # }
  
  
  ## 7. Check metric.eval arguments -------------------------------------------
  metric.eval <- unique(metric.eval)
  if (ms.format@data.type == "binary") {
    avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                              , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC'
                              , 'BOYCE', 'MPA')
  } else if (ms.format@data.type == "ordinal") {
    avail.eval.meth.list <- c("Accuracy", "Recall", "Precision", "F1")
  } else {
    avail.eval.meth.list <- c('RMSE','MSE',"MAE","Rsquared","Rsquared_aj","Max_error")
  }
  .fun_testIfIn(TRUE, paste0("metric.eval with ", ms.format@data.type, " data type"), metric.eval, avail.eval.meth.list)
  
  
  if (!is.null(seed.val)) { set.seed(seed.val) }
  if (is.null(var.import)) { var.import = 0 }
  
  return(list(models = models,
              weights = weights,
              var.import = var.import,
              metric.eval = metric.eval,
              prevalence = prevalence,
              seed.val = seed.val))
}


# ---------------------------------------------------------------------------- #

.MS_Modeling.summary <- function(ms.format, calib.lines, models, models.pa = NULL)
{
  cat("\n\n")
  .bm_cat(paste(ms.format@sp.name, "Modeling Summary"))
  cat("\n", ncol(ms.format@data.env.var), " environmental variables (", colnames(ms.format@data.env.var), ")")
  nb.eval.rep <- ncol(calib.lines) / ifelse(inherits(ms.format, "BIOMOD.formated.data.PA"), ncol(ms.format@PA.table), 1)
  cat("\nNumber of evaluation repetitions :", nb.eval.rep)
  cat("\nModels selected :", models, "\n")
  if (is.null(models.pa)) {
    nb.runs = ncol(calib.lines) * length(models)
  } else {
    nb.runs = length(which(
      sapply(unlist(models.pa), function(x) grepl(colnames(calib.lines), pattern = x))
    ))
  }
  cat("\nTotal number of model runs:", nb.runs, "\n")
  .bm_cat()
}

