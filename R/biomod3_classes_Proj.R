## --------------------------------------------------------------------------- #
## MS.projection.out --------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name MS.projection.out
##' @aliases MS.projection.out-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_Projection()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_Projection}}, and used by 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' 
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the simulation set
##' @slot proj.name a \code{character} corresponding to the projection name
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory variables
##' @slot coord a 2-columns \code{matrix} or \code{data.frame} containing the corresponding 
##' \code{X} and \code{Y} coordinates used to project the species distribution model(s)
##' @slot scale.models a \code{logical} value defining whether models have been rescaled or 
##' not
##' @slot models.projected a \code{vector} containing names of projected models
##' @slot models.out a \code{\link{BIOMOD.stored.data}} object
##' @slot type a \code{character} corresponding to the class of the \code{val} slot of the 
##' \code{proj.out} slot
##' @slot data.type a \code{character} corresponding to the data type
##' @slot proj.out a \code{\link{BIOMOD.stored.data}} object
##' 
##' @param x a \code{\link{MS.projection.out}} object
##' @param object a \code{\link{MS.projection.out}} object
##' @param coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' @param plot.output (\emph{optional, default} \code{facet}) a character
##'   determining the type of output: with \code{plot.output = 'list'} the
##'   function will return a list of plots (one plot per model) ; with 'facet' ;
##'   with \code{plot.output = 'facet'} the function will return a single plot
##'   with all asked projections as facet.
##' @param do.plot (\emph{optional, default} \code{TRUE}) a boolean determining
##'   whether the plot should be displayed or just returned.
##' @param std (\emph{optional, default} \code{TRUE}) a boolean controlling the
##'   limits of the color scales. With \code{std = TRUE} color scales are
##'   displayed between 0 and 1 (or 1000). With \code{std = FALSE} color scales
##'   are displayed between 0 and the maximum value observed.
##' @param scales (\emph{optional, default} \code{fixed}) a character
##'   determining whether x and y scales are shared among facet. Argument passed
##'   to \code{\link[ggplot2:facet_wrap]{facet_wrap}}. Possible values: 'fixed', 'free_x',
##'   'free_y', 'free'.
##' @param size (\emph{optional, default} \code{0.75}) a numeric determing the
##'   size of points on the plots and passed to
##'   \code{\link[ggplot2:geom_point]{geom_point}}.
##' @param ... additional parameters to be passed to \code{\link{get_predictions}} 
##' to select the models that will be plotted
##'           
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}
##' @family Toolbox objects
##' 
##' 
##' 
##' @importFrom grDevices colorRampPalette colors dev.new gray rainbow
##' @importFrom graphics layout legend par points polygon text
##' @importFrom ggplot2 scale_colour_viridis_c scale_fill_viridis_c
##' 
##' @name MS.projection.out-class
##' @rdname MS.projection.out
##' @export
##' 

# 5.1 Class Definition  -----------------------------------

setClass("MS.projection.out",
         representation(modeling.id = 'character',
                        proj.name = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        coord = 'data.frame',
                        scale.models = 'logical',
                        models.projected = 'list',
                        models.out = 'BIOMOD.stored.data',
                        type = 'character',
                        data.type = 'character',
                        proj.out = 'BIOMOD.stored.data'),
         prototype(modeling.id = '',
                   proj.name = '',
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   coord = data.frame(),
                   scale.models = TRUE,
                   models.projected = list(),
                   type = '',
                   data.type = "binary"),
         validity = function(object){ return(TRUE) })


# Other functions ---------------------------------------------------------
## plot.MS.projection.out -------------------------------------------------
##' 
##' @rdname MS.projection.out
##' @export
##' @importFrom terra global
##' @param maxcell maximum number of cells to plot. Argument transmitted to \code{\link[terra]{plot}}.
##' 

setMethod('plot', signature(x = 'MS.projection.out', y = "missing"),
          function(x,
                   coord = NULL,
                   plot.output, # list or facet
                   do.plot = TRUE, # whether plots are displayed or just returned
                   std = TRUE, # limits between 0 and 1000 or between 0 and max
                   scales, # transmitted to facet_wrap
                   size, # size of points transmitted to geom_point
                   maxcell = 5e5, # max number of cells to plot. Transmitted to terra::plot
                   ...)
          {
            # extraction of projection happens in argument check
            args <- .plot.MS.projection.out.check.args(x,
                                                           coord = coord,
                                                           plot.output = plot.output, # list or facet
                                                           do.plot = do.plot,
                                                           std = std,
                                                           scales = scales,
                                                           size = size,
                                                           maxcell = maxcell,
                                                           ...)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            
            ### Plot SpatRaster ---------------------------------------------------------
            
            if (inherits(proj,"SpatRaster")) {
              maxi <- ifelse(max(global(proj, "max", na.rm = TRUE)$max) > 1, 1000, 1)
              if (x@data.type != "binary") { maxi <- max(global(proj, "max", na.rm = TRUE)$max) }
              if (std) {
                limits <-  c(0, maxi)
              } else {
                limits <- NULL
              }
              
              # if(x@data.type == "ordinal"){
              #   breaks = 1:maxi
              #   labels = c("a","b","c","d")
              # } else {
              #   breaks = waiver()
              #   labels = waiver() 
              # }
              
              if (plot.output == "facet") {
                g <- ggplot() +
                  tidyterra::geom_spatraster(data = proj,
                                             maxcell = maxcell) +
                  scale_fill_viridis_c(NULL, limits = limits) +
                  facet_wrap(~lyr)
              } else if (plot.output == "list") {
                g <- lapply(names(proj), function(thislayer){
                  ggplot() +
                    tidyterra::geom_spatraster(data = subset(proj, thislayer),
                                               maxcell = maxcell) +
                    scale_fill_viridis_c(NULL, limits = limits) +
                    ggtitle(thislayer)
                })
              }
            } else {
              ### Plot data.frame  -----------------------------------------------------
              maxi <- ifelse(max(proj$pred) > 1, 1000, 1)
              if (x@data.type != "binary") { maxi <- max(proj$pred, na.rm = TRUE) }
              if (std) {
                limits <-  c(0,maxi)
              } else {
                limits <- NULL
              }
              plot.df <- merge(proj, coord, by = c("points"))
              if (plot.output == "facet") {
                g <- ggplot(plot.df)+
                  geom_point(aes(x = x, y = y, color = pred), size = size) +
                  scale_colour_viridis_c(NULL, limits = limits) +
                  facet_wrap(~full.name)
              } else if (plot.output == "list"){
                g <- lapply(unique(plot.df$full.name), function(thislayer) {
                  ggplot(subset(plot.df, plot.df$full.name == thislayer)) +
                    geom_point(aes(x = x, y = y, color = pred), size = size) +
                    scale_colour_viridis_c(NULL, limits = limits) +
                    ggtitle(thislayer)
                })
              }
            }
            if (do.plot) {
              show(g)
            } 
            return(g)
          }
)

### .plot.MS.projection.out.check.args ----------------------------------

.plot.MS.projection.out.check.args <- function(x, coord, plot.output # list or facet
                                                   , do.plot, std, scales, size, ...)
{
  proj <- get_predictions(x, ...)
  
  ## 1 - check for tidyterra ----------------------
  if (inherits(proj, "SpatRaster")) {
    if (!requireNamespace("tidyterra")) {
      stop("Package `tidyterra` is missing. Please install it with `install.packages('tidyterra')`.")
    }
  }
  
  ## 2 - plot.output----------------------
  if (missing(plot.output)) {
    plot.output <- "facet"
  } else {
    .fun_testIfIn(TRUE, "plot.output", plot.output, c("facet", "list"))
  }
  
  ## 3 - do.plot ----------------------
  stopifnot(is.logical(do.plot))
  
  ## 4 - std ----------------------
  stopifnot(is.logical(std))
  
  ## 5 - check scales for facet_wrap -------------------------------
  if (missing(scales)) {
    scales <- "fixed"
  } else {
    .fun_testIfIn(TRUE, "scales", scales, c("fixed","free","free_x","free_y"))
  }
  
  ## 6 - check coord if x is a data.frame -------------------------------
  if (inherits(proj, 'data.frame')) {
    npred <- length(unique(proj$points))
    
    if (nrow(x@coord) > 0) {
      if (!is.null(coord)) {
        cat("! ignoring argument `coord` as coordinates were already given to BIOMOD_Projection")
      }
      coord <- x@coord
    }
    
    if (nrow(x@coord) == 0 & is.null(coord)) {
      stop("missing coordinates to plot with a data.frame. Either give argument `coord` to plot or argument `new.env.xy` to BIOMOD_Projection")
    } else if (!inherits(coord, c("data.frame","matrix"))) {
      stop("`coord` must be a data.frame or a matrix.")
    } else if (ncol(coord) != 2) {
      stop("`coord` must have two columns.")
    } else if (nrow(coord) != npred) {
      stop("`coord` must have as many rows as the number of predictions (", npred, ").")
    } else {
      coord <- as.data.frame(coord)
      colnames(coord) <- c("x", "y")
      coord$points <- seq_len(npred)
    }
  }
  
  if (missing(size)) {
    size <- 0.75
  } 
  
  ## 7 - check size -------------------------------
  if (inherits(proj, 'data.frame')) {
    .fun_testIfPosNum(TRUE, "size", size)
  }
  
  return(list(proj = proj,
              coord = coord,
              plot.output = plot.output,
              do.plot = do.plot,
              std = std,
              scales = scales, 
              size = size))
}

## show.MS.projection.out -------------------------------------------------
##' 
##' @rdname MS.projection.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('MS.projection.out'), function(object)
{
  .bm_cat("MS.projection.out")
  cat("\nProjection directory :", paste0(object@dir.name, "/", object@sp.name, "/", object@proj.name), fill = .Options$width)
  cat("\n")
  cat("\nsp.name :", object@sp.name, fill = .Options$width)
  cat("\nexpl.var.names :", object@expl.var.names, fill = .Options$width)
  cat("\n")
  cat("\nmodeling.id :", object@modeling.id , "(", object@models.out@link , ")", fill = .Options$width)
  cat("\nmodels.projected :", toString(object@models.projected), fill = .Options$width)
  df.info <- .extract_projlinkInfo(object)
  if (any(df.info$type == "bin")) {
    available.metric <- unique(subset(df.info, df.info$type == "bin")$metric)
    cat("\navailable binary projection :", toString(available.metric), fill = .Options$width)
  }
  if (any(df.info$type == "filt")) {
    available.metric <- unique(subset(df.info, df.info$type == "filt")$metric)
    cat("\navailable filtered projection :", toString(available.metric), fill = .Options$width)
  }
  .bm_cat()
})

## get_projected_models.MS.projection.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_projected_models", "MS.projection.out",
          function(obj, full.name = NULL, PA = NULL, run = NULL, algo = NULL
                   , merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL)
          {
            out <- obj@models.projected
            if (length(grep("EM|merged", out)) > 0) {
              keep_ind <- .filter_outputs.vec(out, obj.type = "em", subset.list = list(full.name = full.name
                                                                                       , merged.by.PA = merged.by.PA
                                                                                       , merged.by.run = merged.by.run
                                                                                       , merged.by.algo = merged.by.algo
                                                                                       , filtered.by = filtered.by
                                                                                       , algo = algo))
            } else {
              keep_ind <- .filter_outputs.vec(out, obj.type = "mod", subset.list = list(full.name = full.name, PA = PA
                                                                                        , run = run, algo = algo))
            }
            out <- out[keep_ind]
            return(out)
          }
)

## free.MS.projection.out --------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' @importFrom terra rast

setMethod('free', signature('MS.projection.out'), function(obj)
{
  if (inherits(obj@proj.out, "BIOMOD.stored.data.frame")) {
    obj@proj.out@val  <- data.frame()
  } else if (inherits(obj@proj.out, "BIOMOD.stored.SpatRaster")) {
    obj@proj.out@val <- wrap(rast(matrix()))
  } else {
    obj@proj.out@val <- NULL
  }
  obj@proj.out@inMemory <- FALSE
  return(obj)
})

## get_predictions.MS.projection.out ---------------------------------------
# (the method is used for EM as well)
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "MS.projection.out",
          function(obj, metric.binary = NULL, metric.filter = NULL
                   , full.name = NULL, PA = NULL, run = NULL, algo = NULL
                   , merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, 
                   model.as.col = FALSE, ...)
          {
            # extract layers from obj@proj.out@link concerned by metric.filter or metric.binary
            selected.layers <- .extract_selected.layers(obj, 
                                                        metric.binary = metric.binary,
                                                        metric.filter = metric.filter)
            out <- load_stored_object(obj@proj.out, layer = selected.layers)
            
            # subselection of models_selected
            if (obj@type == "SpatRaster") {
              if (length(grep("EM|merged", names(out))) > 0) {
                keep_layers <- .filter_outputs.vec(names(out), obj.type = "em", 
                                                   subset.list = list(full.name =  full.name
                                                                      , merged.by.PA = merged.by.PA
                                                                      , merged.by.run = merged.by.run
                                                                      , merged.by.algo = merged.by.algo
                                                                      , filtered.by = filtered.by
                                                                      , algo = algo))
              } else {
                keep_layers <- .filter_outputs.vec(names(out), obj.type = "mod",
                                                   subset.list = list(full.name =  full.name, PA = PA
                                                                      , run = run, algo = algo))
              }
              out <- subset(out, keep_layers)
            } else {
              if (length(grep("EM|merged", colnames(out))) > 0) {
                keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name
                                                                         , merged.by.PA = merged.by.PA
                                                                         , merged.by.run = merged.by.run
                                                                         , merged.by.algo = merged.by.algo
                                                                         , filtered.by = filtered.by
                                                                         , algo = algo))
              } else {
                keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name, PA = PA
                                                                         , run = run, algo = algo))
              }
              out <- out[keep_lines, ]
              if (model.as.col) {
                out <- .transform_model.as.col(out)
              }
            }
            
            return(out)
          }
)

