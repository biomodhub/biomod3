###################################################################################################
##' @name bm_SpeciesParameters
##' @author Helene Blancheteau
##' 
##' @title Set the species parameters
##' 
##' @description This function allows to create the differents parameters list for all the species in the same time.
##' 
##' @param resp.name a \code{character vector} of your species name.
##' 
##' @param PA.nb.rep (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection, an \code{integer} corresponding to the number of sets 
##' (repetitions) of pseudo-absence points that will be drawn
##' @param PA.strategy (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection, a \code{character} defining the strategy that will be used to 
##' select the pseudo-absence points. Must be \code{random}, \code{sre}, \code{disk} or 
##' \code{user.defined} (see Details)
##' @param PA.nb.absences (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection, and \code{PA.strategy = 'random'} or \code{PA.strategy = 'sre'} 
##' or \code{PA.strategy = 'disk'}, an \code{integer} corresponding to the number of pseudo-absence 
##' points that will be selected for each pseudo-absence repetition (true absences included). \cr
##' It can also be a \code{vector} of the same length as \code{PA.nb.rep} containing \code{integer} 
##' values corresponding to the different numbers of pseudo-absences to be selected
##' @param PA.sre.quant (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'sre'}, a \code{numeric} between \code{0} 
##' and \code{0.5} defining the half-quantile used to make the \code{sre} pseudo-absence selection 
##' (see Details)
##' @param PA.dist.min (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' minimal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in meters, see Details)
##' @param PA.dist.max (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' maximal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in meters, see Details)
##' @param PA.fact.aggr (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'random'} or \code{strategy = 'disk'}, a \code{integer} defining the factor of aggregation to reduce the resolution
##' @param PA.user.table (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'user.defined'}, a \code{matrix} or 
##' \code{data.frame} with as many rows as \code{resp.var} values, as many columns as 
##' \code{PA.nb.rep}, and containing \code{TRUE} or \code{FALSE} values defining which points 
##' will be used to build the species distribution model(s) for each repetition (see Details)
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
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' @param em.by a \code{character} corresponding to the way kept models will be combined to build 
##' the ensemble models, must be among \code{all}, \code{algo}, \code{PA}, \code{PA+algo}, 
##' \code{PA+run}
##' @param metric.select a \code{vector} containing evaluation metric names to be used together with 
##' \code{metric.select.thresh} to exclude single models based on their evaluation scores 
##' (for ensemble methods like probability weighted mean or committee averaging). Must be among  
##' \code{all} (same evaluation metrics than those of \code{bm.mod}), \code{user.defined} 
##' (and defined through \code{metric.select.table}) or \code{POD}, \code{FAR}, \code{POFD}, 
##' \code{SR}, \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, 
##' \code{ORSS}, \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}, \code{RMSE}, \code{MAE}, 
##' \code{MSE}, \code{Rsquared}, \code{Rsquared_aj}, \code{Max_error}, \code{Accuracy}, \code{Recall},
##' \code{Precision}, \code{F1}
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
##' @param EMci.alpha (\emph{optional, default} \code{0.05}) \cr 
##' A \code{numeric} value corresponding to the significance level to estimate confidence interval
##' @param EMwmean.decay (\emph{optional, default} \code{proportional}) \cr A
##' value defining the relative importance of the weights (if \code{'EMwmean'}
##' was given to argument \code{em.algo}). A high value will strongly
##' discriminate \emph{good} models from the \emph{bad} ones (see Details),
##' while \code{proportional} will attribute weights proportionally to the
##' models evaluation scores
##'   
##' 
##' @return
##' 
##' A \code{list} with the differents set of parameters fo all the species of resp.name :
##' \code{params.PA}, \code{params.CV}, \code{params.OPT} and \code{params.EM}.
##' 
##' @examples
##' library(biomod2)
##' data(DataSpecies)
##' names(DataSpecies)
##' 
##' sp_parameters <- bm_SpeciesParameters(resp.name = names(DataSpecies)[3:8],
##'                                       PA.nb.rep = 3, 
##'                                       PA.nb.absences = c(100,200,500), 
##'                                       PA.strategy = "random",
##'                                       CV.nb.rep = 1, 
##'                                       CV.strategy = 'kfold', 
##'                                       CV.k = 3,
##'                                       OPT.strategy = 'bigboss',
##'                                       models.chosen = "all",
##'                                       metric.select = 'TSS', 
##'                                       metric.select.thresh = 0.7)
##'
##' 
##' @export
##' 
##' 
###################################################################################################

bm_SpeciesParameters <- function(resp.name,
                               PA.nb.rep = NULL, PA.nb.absences = NULL, PA.strategy = NULL, PA.dist.min = NULL, PA.dist.max = NULL,
                               PA.sre.quant = 0.025, PA.fact.aggr = NULL, PA.user.table = NULL,
                               CV.nb.rep = 1, CV.strategy = 'random', CV.perc = 0.7, CV.k = NULL, CV.balance = NULL,
                               CV.env.var = NULL, CV.strat = NULL, CV.user.table = NULL, CV.do.full.models = TRUE,
                               OPT.strategy = 'default', OPT.user.val, OPT.user.base, OPT.user,
                               models.chosen = "all", em.by, metric.select = NULL, metric.select.thresh = NULL, metric.select.table = NULL,
                               metric.select.dataset = NULL, EMci.alpha = 0.05, EMwmean.decay = "proportional") {
  
  .bm_cat("Attribute parameters to species")
  
  ## 1. params.PA
  cat("\n \t Setting of params.PA")
  
  .fun_testIfIn(TRUE, "strategy", PA.strategy, c("random", "sre", "disk", "user.defined"))
  
  ## 4. Check sre.quant argument ----------------------------------------------
  if (PA.strategy == 'SRE' && (PA.sre.quant >= 0.5 || PA.sre.quant < 0)) {
    stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
  }
  
  ## 5. Check dist.min and dist.max arguments ---------------------------------
  if (!is.null(PA.strategy) && PA.strategy == 'disk') {
    if (!is.null(PA.dist.min) && PA.dist.min < 0) {
      dist.min <- 0
    }
    if (!is.null(PA.dist.max) && PA.dist.max < 0) {
      dist.max <- NULL
    }
    if (!is.null(PA.dist.max) && !is.null(PA.dist.min) && PA.dist.min >= PA.dist.max) {
      stop("dist.min >= dist.max")
    }
  }
  
  ## 6. Check nb.absences argument --------------------------------------------
  if (!is.null(PA.strategy) && PA.strategy != "user.defined") {
    if (is.null(PA.nb.absences)) {
      stop("You must give the number of pseudo absences you want")
    } else {
      if (length(PA.nb.absences) > 1) {
        if (length(PA.nb.absences) != PA.nb.rep) {
          stop("You must give one value for pseudo absences, or as many as the number of repetitions")
        } else if (length(unique(PA.nb.absences)) == 1) {
          PA.nb.absences = unique(PA.nb.absences)
        }
      }
    }
  }
  
  ## 7. Check user.table argument --------------------------------------------
  if (!is.null(PA.strategy) && PA.strategy == "user.defined") {
    if (is.null(PA.user.table)) {
      stop("You must give a table defining the pseudo absences you want")
    } else {
      if (!(is.matrix(PA.user.table) | is.data.frame(PA.user.table))) {
        stop("\n PA.user.table must be a matrix or a data.frame")
      }
      colnames(PA.user.table) <- paste0("PA", 1:ncol(PA.user.table))
      PA.nb.absences <- nrow(PA.user.table)
    }
  }
  
  params.PA <- list("PA.nb.rep" = PA.nb.rep, "PA.nb.absences" = PA.nb.absences, 
                     "PA.strategy" = PA.strategy, 
                     "PA.dist.min" = PA.dist.min, "PA.dist.max" = PA.dist.max, "PA.fact.aggr" = PA.fact.aggr,
                     "PA.sre.quant" = PA.sre.quant, 
                     "PA.user.table" = PA.user.table)
  params.PA <- rep(list(params.PA), length(resp.name))
  names(params.PA) <- resp.name
  
  ## 2. params.CV
  cat("\n \t Setting of params.CV")
  
  ## 1. Check strategy argument -------------------------------------
  .fun_testIfIn(TRUE, "strategy", CV.strategy, c("random", "kfold", "block", "strat", "env", "user.defined"))
  
  ## 2.a Check nb.rep / perc argument -------------------------------
  if (CV.strategy %in% c("random", "kfold")) {
    .fun_testIfPosInt(TRUE, "nb.rep", CV.nb.rep)
    if (CV.nb.rep < 1) { stop("nb.rep must be an integer >= 1") }
    
    if (CV.strategy == "random") {
      if (is.null(CV.perc)) {
        stop("perc (or CV.perc) is required when strategy = 'random'")
      }
      .fun_testIf01(TRUE, "perc", CV.perc)
      if (CV.perc < 0.5) {
        warning("You chose to allocate more data to validation than to calibration of your model
                (perc<0.5)\nMake sure you really wanted to do that. \n", immediate. = TRUE)
      } else if (CV.perc == 1) {
        CV.nb.rep <- 0
        warning(paste0("The models will be evaluated on the calibration data only "
                       , "(nb.rep=0 and no independent data) \n\t "
                       , "It could lead to over-optimistic predictive performances.\n")
                , immediate. = TRUE)
      }
    }
  }
  
  ## 2.b Check k argument -------------------------------------------
  if (CV.strategy %in% c("kfold", "strat", "env")) {
    .fun_testIfPosInt(TRUE, "k", CV.k)
    if (CV.k < 2) { stop("k must be an integer >= 2") }
  }
  
  ## 3. Check balance / strat argument ------------------------------
  if (CV.strategy %in% c("strat", "env")) {
    .fun_testIfIn(TRUE, "balance", CV.balance, c("presences","absences"))
    
    if (CV.strategy == "strat") {
      .fun_testIfIn(TRUE, "strat", CV.strat, c("x", "y", "both"))
    }
  }
  
  ## 4. Check user.table argument -----------------------------------
  if (CV.strategy == "user.defined") {
    if (is.null(CV.user.table)) {
      stop("user.table must be a matrix or a data.frame") 
    } else {
      .fun_testIfInherits(TRUE, "user.table", CV.user.table, c("matrix", "data.frame"))
      if (inherits(user.table, 'data.frame')) {
        CV.user.table <- as.matrix(CV.user.table)
      }
    }
  }
  
  params.CV <- list("CV.nb.rep" = CV.nb.rep, "CV.strategy" = CV.strategy,
                    "CV.perc" = CV.perc, "CV.k" = CV.k,
                    "CV.balance" = CV.balance, "CV.env.var" = CV.env.var, "CV.strat" = CV.strat,
                    "CV.user.table" = CV.user.table, "CV.do.full.models" = CV.do.full.models)
  params.CV <- rep(list(params.CV), length(resp.name))
  names(params.CV) <- resp.name
  
  ## 3. params.OPT
  cat("\n \t Setting of params.OPT")
  .fun_testIfIn(TRUE, "OPT.strategy", OPT.strategy, c("bigboss", "default", "tuned", "user.defined"))
  
  if (OPT.strategy == "user.defined"){
    if (is.null(OPT.user) && is.null(OPT.user.val)){
      stop("With user.defined OPT.strategy, you need OPT.user or OPT.user.val")
    }
  } else {
    OPT.user.val <- NULL
    OPT.user.base <- NULL
    OPT.user <- NULL
  }
  
  params.OPT <- list("OPT.strategy" = OPT.strategy,
                                "OPT.user.val" = OPT.user.val,
                                "OPT.user.base" = OPT.user.base,
                                "OPT.user" = OPT.user)
  params.OPT <- rep(list(params.OPT), length(resp.name))
  names(params.OPT) <- resp.name
  
  ## 4. params.EM
  cat("\n \t Setting of params.EM")
  
  em.by.avail <- c('PA', 'algo', 'all', 'PA+run', 'PA+algo')
  if(missing(em.by)){
    em.by <- "all"
    cat("\n \t\t! `em.by` automatically set to 'all'")
  }

  .fun_testIfIn(TRUE, "em.by", em.by, em.by.avail)
  
  metric.select <- unique(metric.select)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC'
                            , 'BOYCE', 'MPA',
                            "Accuracy", "Recall", "Precision", "F1",
                            'RMSE','MSE',"MAE","Rsquared","Rsquared_aj","Max_error")

  .fun_testIfIn(TRUE, "metric.select", metric.select, avail.eval.meth.list)
  
  if (!is.null(metric.select.dataset)) {
    metric.select.dataset.available <- c("calibration", "validation", "evaluation")
    .fun_testIfIn(TRUE, "metric.select.dataset", metric.select.dataset, metric.select.dataset.available)
  }
  
  .fun_testIfPosNum(TRUE, "EMci.alpha", EMci.alpha)
  if (EMci.alpha <= 0 | EMci.alpha >= 0.5) {
    stop("EMci.alpha must be a numeric between 0 and 0.5")
  }
  
  if ((!is.numeric(EMwmean.decay) &&
       !is.character(EMwmean.decay) &&
       !is.function(EMwmean.decay)) ||
      (is.numeric(EMwmean.decay) && EMwmean.decay < 0) ||
      (is.character(EMwmean.decay) && EMwmean.decay != 'proportional')) {
    stop("'EMwmean.decay' should be either 'proportional', a numeric value > 0 or a function")
  }
  
  
  
  params.EM <- list("models.chosen" = models.chosen, "em.by" = em.by, 
                    "metric.select" = metric.select, "metric.select.thresh" = metric.select.thresh,
                    "metric.select.table" = metric.select.table, "metric.select.dataset" = metric.select.dataset, 
                    "EMci.alpha" = EMci.alpha, "EMwmean.decay" = EMwmean.decay)

  params.EM <- rep(list(params.EM), length(resp.name))
  names(params.EM) <- resp.name
  
  
  
  .bm_cat("Done")
  return(list("params.PA" = params.PA,
              "params.CV" = params.CV, 
              "params.OPT" = params.OPT,
              "params.EM" = params.EM))
}
