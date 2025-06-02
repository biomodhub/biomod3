## --------------------------------------------------------------------------- #
## 1. MS.formated.data ---------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name MS.formated.data
##' @aliases MS.formated.data-class
##' @author Damien Georges
##' 
##' @title \code{MS_FormatingData()} output object class
##' 
##' @description Class returned by \code{\link{MS_FormatingData}}
##' \code{\link{MS_Modeling}}
##' 
##' 
##' @param object a \code{\link{MS.formated.data}} object
##' 
##' @slot data.type a \code{character} corresponding to the data type
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot ms.project a \code{character} correspondinf to the name of the multispecies project
##' @slot sp.name a \code{character vector} corresponding to the species name
##' @slot coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates
##' @slot data.species a \code{data.frame} containing the species observations (\code{0}, \code{1} or 
##' \code{NA})
##' @slot data.env.var a \code{data.frame} containing explanatory variables
##' @slot data.mask a \code{\link[terra:rast]{SpatRaster}} object containing the mask of the 
##' studied area
##' @slot PA.strategy a \code{character vector} corresponding to the pseudo-absence selection strategy
##' @slot PA.table a \code{list} with the \code{data.frames} containing the corresponding table of selected 
##' pseudo-absences (indicated by \code{TRUE} or \code{FALSE}) from the \code{pa.tab} list 
##' element returned by the \code{\link{bm_PseudoAbsences}} function
##' @slot has.data.eval a \code{logical} value defining whether evaluation data is given
##' @slot has.filter.raster a \code{logical} value defining whether filtering have been done or not
##' @slot eval.coord (\emph{optional, default} \code{NULL}) \cr 
##' A 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates for evaluation data
##' @slot eval.data.species (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the species observations (\code{0}, \code{1} or \code{NA}) for 
##' evaluation data
##' @slot eval.data.env.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{data.frame} containing explanatory variables for evaluation data
##' @slot biomod3.version a \code{character} corresponding to the biomod3 version
##' 
##' 
##' @rdname MS.formated.data
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("MS.formated.data",
         representation(data.type = 'character',
                        dir.name = 'character',
                        ms.project = 'character',
                        sp.name = 'character',
                        coord = "data.frame",
                        data.species = "ANY",
                        data.env.var = "data.frame",
                        data.mask = "list",
                        PA.strategy = "character",
                        PA.table = "list",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "ANY",
                        eval.data.env.var = "data.frame",
                        has.filter.raster = "logical",
                        biomod3.version = "character"),
         validity = function(object) { 
           check.data.mask <- suppressWarnings(
             all(sapply(object@data.mask, function(x) inherits(x, "PackedSpatRaster")))
           )
           if (check.data.mask) {
             return(TRUE)
           } else {
             return(FALSE)
           }
         })

### plot.MS.formated.data (doc) --------------------------------------------------
##' 
##' @rdname plot
##' @docType methods
##' @title \code{plot} method for \code{\link{MS.formated.data}} object class
##' 
##' @description Plot the spatial distribution of presences, absences and 
##' pseudo-absences among the different potential dataset (calibration, 
##' validation and evaluation). Available only if coordinates were given to 
##' \code{\link{MS_FormatingData}}.
##' 
##' 
##' @param x a \code{\link{MS.formated.data}} object. Coordinates must be available to be able to use \code{plot}.
##' @param plot.type a \code{character}, either \code{'points'} (\emph{default}) 
##' or \code{'raster'} (\emph{if environmental variables were given as a raster}). 
##' With \code{plot.type = 'points'} occurrences will be represented as points
##' (better when using fine-grained data). With \code{plot.type = 'raster'}
##' occurrences will be represented as a raster (better when using coarse-grained
##' data)
##' @param plot.output a \code{character}, either \code{'facet'} (\emph{default}) 
##' or \code{'list'}. \code{plot.output} determines whether plots are returned
##' as a single facet with all plots or a \code{list} of individual plots
##' (better when there are numerous graphics)
##' @param plot.eval (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether evaluation data should be added to the plot or not
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether the plot is to be rendered or not
##' @param point.size a \code{numeric} to adjust the size of points when
##'  \code{plot.type = 'points'}.
##' 
##' @return a \code{list} with the data used to generate the plot and a
##' \code{ggplot2} object 
##' 
##' @importFrom terra rast minmax crds ext rasterize
##' @importFrom ggplot2 ggplot aes scale_color_manual scale_shape_manual scale_fill_manual guides xlim ylim ggtitle facet_wrap theme guide_legend after_stat xlab ylab element_blank element_rect geom_point
##' 
##' @export




setMethod('plot', signature(x = 'MS.formated.data', y = "missing"),
          function(x,
                   #calib.lines = NULL,
                   plot.type,
                   plot.output, 
                   #PA,
                   #run,
                   plot.eval,
                   point.size = 1.5,
                   do.plot = TRUE){
            
            args <- .plot.BIOMOD.formated.data.check.args(x = x,
                                                          #calib.lines = calib.lines,
                                                          plot.type = plot.type,
                                                          plot.output = plot.output, 
                                                          #PA = PA,
                                                          #run = run,
                                                          plot.eval = plot.eval,
                                                          do.plot = do.plot)
            
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            # 1 - extract data for all required data ----------------------
            
            nameFolder <- file.path(x@dir.name, x@ms.project, ".BIOMOD_DATA", "single.formated.data")
            table <- foreach(sp = x@sp.name, .combine = merge) %do% {
              sfd <- get(load(file.path(nameFolder, paste0(sp,".sfd"))))
              plot_sfd <- biomod2::plot(sfd, do.plot = FALSE)
              df <- plot_sfd$data.plot$data
              names(df) <- c(sp, "dataset", "x", "y")
              return(df)
            }
            
            # 2- define colors and breaks ------------------------------------
            data_breaks <- c(9,10, 11, 12, # presence
                             19,20, 21, 22, # absence
                             29,30, 31,         # pseudo-absences
                             1)              # background
            data_labels <- c("9" = "**Presences**",
                             "10" = "Presences (calibration)",
                             "11" = "Presences (validation)",
                             "12" = "Presences (evaluation)",
                             "19" = "**True Absences**",
                             "20" = "True Absences (calibration)",
                             "21" = "True Absences (validation)",
                             "22" = "True Absences (evaluation)",
                             "29" = "**Pseudo-Absences**",
                             "30" = "Pseudo-Absences (calibration)",
                             "31" = "Pseudo-Absences (validation)",
                             "1" = "Background")
            data_labels_facet <- c("9" = "**Presences**",
                                   "10" = "calibration",
                                   "11" = "validation",
                                   "12" = "evaluation",
                                   "19" = "**True Absences**",
                                   "20" = "calibration",
                                   "21" = "validation",
                                   "22" = "evaluation",
                                   "29" = "**Pseudo-Absences**",
                                   "30" = "calibration",
                                   "31" = "validation",
                                   "1" = NA) # background
            
            data_colors <- c("9" = NA,
                             "10" = "#004488",
                             "11" = "#6699CC",
                             "12" = "#6699CC",
                             "19" = NA,
                             "20" = "#994455",
                             "21" = "#EE99AA",
                             "22" = "#EE99AA",
                             "29" = NA,
                             "30" = "#997700",
                             "31" = "#EECC66",
                             "1" = "grey70")
            
            shape_fit <- 16
            shape_eval <- 17
            data_shape <- c("9" = NA,
                            "10" = shape_fit,
                            "11" = shape_fit,
                            "12" = shape_eval,
                            "19" = NA,
                            "20" = shape_fit,
                            "21" = shape_fit,
                            "22" = shape_eval,
                            "29" = NA,
                            "30" = shape_fit,
                            "31" = shape_fit,
                            "1" = NA)
            data_alpha <- c("9" = 0,
                            "10" = 1,
                            "11" = 1,
                            "12" = 1,
                            "19" = 0,
                            "20" = 1,
                            "21" = 1,
                            "22" = 1,
                            "29" = 0,
                            "30" = 1,
                            "31" = 1,
                            "1"  = 0)
            data_background <- "#FFFFFF00"
            
            
            # 3 - prepare plots -------------------------------------------------------
            this_mask_eval <- rast()
            if(has.mask){
              this_mask <- rast(x@data.mask[["calibration"]])
              this_mask_eval <- this_mask
            } else {
              this_mask <- rast()
            }
            if(has.mask.eval){
              this_mask_eval <- rast(x@data.mask[["evaluation"]])
            }
            if(has.mask | has.mask.eval){
              plot_mask <- foreach(this_dataset = unique(table$dataset), 
                                   .combine = 'c') %do% {
                                     if(this_dataset == "Evaluation dataset"){
                                       return(this_mask_eval)
                                     } else {
                                       return(this_mask)
                                     }
                                   }
              names(plot_mask) <- unique(table$dataset)
            }
            ## 3.1 Raster plot --------------------------------------------------------
            
            if(plot.type == "raster"){
              
              rast.plot <- foreach(this_dataset = unique(table$dataset), .combine = 'c') %do% {
                this_rast  <-
                  rasterize(as.matrix(table[table$dataset == this_dataset, c("x", "y")]), 
                            plot_mask[[this_dataset]],
                            table[, x@sp.name], background = 1)
                names(this_rast) <- x@sp.name
                this_rast*this_mask
              }
              
              g <- ggplot()+
                tidyterra::geom_spatraster(data = rast.plot,
                                           aes(fill = factor(after_stat(value), data_breaks)))+
                facet_wrap(~lyr)+
                scale_fill_manual(
                  NULL,
                  breaks = data_breaks,
                  values = data_colors,
                  labels = data_labels_facet,
                  na.value = data_background, 
                  drop = FALSE)+
                guides(fill = guide_legend(
                  override.aes = list(alpha = data_alpha),
                  ncol = 3))+
                theme(legend.position = "top",
                      legend.key = element_blank(),
                      legend.background = element_rect(fill = "grey90"),
                      legend.text = ggtext::element_markdown())
              
            if(do.plot){
              print(g)
            }
            return(list("data.vect"  = table,
                        "data.rast"  = rast.plot,
                        "data.label" = data_labels,
                        "data.plot"  = g))
          } else {
            ## 3.2 Points plot --------------------------------------------------------
            
            data.df <- as.data.frame(table)
            data.df <- reshape2::melt(data.df, c("dataset", "x", "y"), value.name = "resp", variable.name = "species")
            base_g <-  ggplot(data.df)
            if(has.mask){
              base_g <- base_g +
                tidyterra::geom_spatraster(data = this_mask, aes(fill = factor(after_stat(value))))
            }
            
            
            g <- base_g +      
              geom_point(aes(x = x, y = y, 
                             color = factor(resp, levels = data_breaks[-12]),
                             shape = factor(resp, levels = data_breaks[-12])), 
                         alpha = 1, size = point.size)+
              facet_wrap(~species)+
              scale_color_manual(
                NULL,
                breaks = data_breaks,
                values = data_colors,
                labels = data_labels_facet,
                drop = FALSE)+
              scale_shape_manual(
                NULL,
                breaks = data_breaks,
                values = data_shape,
                labels = data_labels_facet,
                drop = FALSE)+
              scale_fill_manual(
                guide = "none",
                breaks = data_breaks,
                values = data_colors,
                labels = data_labels,
                na.value = data_background)+
              xlab(NULL)+ ylab(NULL)+
              guides(color = guide_legend(override.aes = list(size = 3),
                                          ncol = 3))+
              theme(legend.position = "top",
                    legend.key = element_blank(),
                    legend.background = element_rect(fill = "grey90"),
                    legend.text = ggtext::element_markdown())
            
            
          }
          if(do.plot){
            print(g)
          }
          return(list("data.table"  = table,
                      "data.label" = data_labels,
                      "data.plot"  = g))
          }

)


.plot.BIOMOD.formated.data.check.args <- function(x,
                                                  calib.lines,
                                                  plot.type,
                                                  plot.output, 
                                                  PA,
                                                  run,
                                                  plot.eval,
                                                  do.plot){
  
  
  ## 1 - check x -----------------------------------------
  .fun_testIfInherits(TRUE, "x", x, c("MS.formated.data"))
  
  
  ## 3 - check plot.eval ----------------------
  if (missing(plot.eval)) {
    plot.eval <- x@has.data.eval
  } else {
    stopifnot(is.logical(plot.eval))
    if(plot.eval & !x@has.data.eval){
      plot.eval <- FALSE
      cat('\n  ! Evaluation data are missing and its plot was deactivated')
    }
  }
  
  ## 4 are proper mask available ? -----------------------
  has.mask <- rast.has.values(rast(x@data.mask[["calibration"]]))
  if(plot.eval){
    has.mask.eval <- length(x@data.mask) > 1
  } else {
    has.mask.eval <- FALSE
  }
  if (has.mask | has.mask.eval) {  
    if (!requireNamespace("tidyterra")) {
      stop("Package `tidyterra` is missing. Please install it with `install.packages('tidyterra')`.")
    }
  } 
  ## 5 - check plot.type  ----------------------
  if (missing(plot.type)) {
    plot.type <- "points"
  } else {
    .fun_testIfIn(TRUE, "plot.type", plot.type, c("raster","points"))
    if ( !has.mask & plot.type == "raster") {
      plot.type <- "points"
      cat("\n ! no raster available, `plot.type` automatically set to 'points'\n")
    }
  }
  
  ## 6 - plot.output----------------------
  if (missing(plot.output)) {
    plot.output <- "facet"
  } else {
    .fun_testIfIn(TRUE, "plot.output", plot.output, c("facet","list"))
  }
  
  if(plot.output == "facet"){
    if(!requireNamespace("ggtext")){
      stop("Package `ggtext` is missing. Please install it with `install.packages('ggtext')`.")
    }
  }
  
  ## 7 - do.plot ----------------------
  # do.plot
  stopifnot(is.logical(do.plot))
  
  ##  9 - check that coordinates are available -------------------------------
  if(nrow(x@coord) == 0){
    stop("coordinates are required to plot BIOMOD.formated.data objects")
  }
  
  ## End - return arguments ----------------------------------------------------
  return(list(x = x,
              #calib.lines = calib.lines,
              plot.type = plot.type,
              plot.output = plot.output, 
              #PA = PA,
              #run = run,
              plot.eval = plot.eval,
              do.plot = do.plot,
              has.mask = has.mask,
              has.mask.eval = has.mask.eval))
}


### show.MS.formated.data  --------------------------------------------------
##' 
##' @rdname MS.formated.data
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('MS.formated.data'),
          function(object)
          {
            .bm_cat("MS.formated.data")
            cat("\ndir.name = ", object@dir.name, fill = .Options$width)
            cat("\nname of the project = ", object@ms.project, fill = .Options$width)
            cat("\nnumber of species = ", length(object@sp.name), fill = .Options$width)
            
            ## Summary of data.species ?? Peut être très lourd
            cat("\n\n\t",
                ncol(object@data.env.var),
                'explanatory variables\n',
                fill = .Options$width)
            print(summary(object@data.env.var))
            
            if (object@has.data.eval) {
              cat("\n\nEvaluation data :", fill = .Options$width)
              # cat("\n\t",
              #     sum(object@eval.data.species, na.rm = TRUE),
              #     'presences, ',
              #     sum(object@eval.data.species == 0, na.rm = TRUE),
              #     'true absences and ',
              #     sum(is.na(object@eval.data.species), na.rm = TRUE),
              #     'undefined points in dataset',
              #     fill = .Options$width)
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            .bm_cat()
          }
)


##' 
##' @importMethodsFrom biomod2 summary
##'

setMethod('summary', signature(object = 'MS.formated.data'),
          function(object, calib.lines = NULL) {
            nameFolder <- file.path(object@dir.name, object@ms.project, ".BIOMOD_DATA", "single.formated.data")
            table <- foreach(sp = object@sp.name, .combine = rbind) %do% {
              sfd <- get(load(file.path(nameFolder, paste0(sp,".sfd"))))
              return(data.frame("species" = sp , summary(sfd, calib.lines = calib.lines)))
            }
            table
          }
)

