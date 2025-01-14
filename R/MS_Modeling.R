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
##' @param ms.format a \code{\link{MS.formated.data}} object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' 
##' @param params.CV a \code{list} with the species names associated to the parameters of Cross-Validation. See BIOMOD_Modeling
##' @param params.OPT a \code{list} with the species names associated to the options of the algorithms. See BIOMOD_Modeling
##' 
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
##' 
##' 
##' @return
##' 
##' A \code{\link{MS.models.out}} object acting as a proxi for the created \code{BIOMOD.models.out}.
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
    scale.models = scale.models,
    models.computed = list(),
    models.failed = list()
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
    MSmodels@models.computed[[sp]] = sfd_models@models.computed
    MSmodels@models.failed[[sp]] = sfd_models@models.failed
  }
  
  ## SAVE MODEL OBJECT ON HARD DRIVE ----------------------------
  name.OUT = paste0(ms.format@ms.project, '.', modeling.id, '.MS.models.out')
  MSmodels@link <- file.path(dir.name, name.OUT)
  assign(x = name.OUT, value = MSmodels)
  save(list = name.OUT, file = MSmodels@link)
  
  .bm_cat("Done")
  return(MSmodels)
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

