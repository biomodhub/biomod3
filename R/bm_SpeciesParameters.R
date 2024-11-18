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
##' @return
##' 
##' A \code{list} with the differents set of parameters fo all the species of resp.name :
##' \code{params.PA}, \code{params.CV}, \code{params.OPT} and \code{params.EM}.
##' 
##' 
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
  if (PA.strategy == 'disk') {
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
  if (PA.strategy != "user.defined") {
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
  if (PA.strategy == "user.defined") {
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
