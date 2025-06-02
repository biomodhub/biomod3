###################################################################################################
##' @name MS_FormatingData
##' @author Helene Blancheteau
##' 
##' @title Format input data for usage in \pkg{biomod3}
##' 
##' @description This function gathers together all input data needed (\emph{xy, 
##' presences/absences, explanatory variables, and the same for evaluation data if available}) to 
##' run \pkg{biomod2} models. It run the function BIOMOD_FormatingData for all the different species.
##' 
##' 
##' @param ms.project.name a \code{character} corresponding to the name of your project for your multispecies modelling
##' @param dir.name (\emph{optional, default} \code{.}) \cr
##' A \code{character} corresponding to the modeling folder
##' @param resp.name a \code{character vector} corresponding to the species name
##' 
##' @param resp.var a \code{vector} or a \code{\link[terra:vect]{SpatVector}} containing your response variable (See Details).
##' @param data.type a \code{character}, corresponding to the response data type to be used, must be either 
##' \code{binary}, \code{count}, \code{multiclass}, \code{ordinal}, \code{relative}, or \code{abundance}. If data.type is not provided,
##' \code{biomod2} will try to guess.
##' 
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' 
##' @param resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be
##' used to build the species distribution model(s)
##' 
##' @param eval.resp.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector}, a \code{\link[terra:vect]{SpatVector}}
##' without associated data (\emph{if presence-only}), 
##' or a \code{\link[terra:vect]{SpatVector}} object containing binary data  
##' (\code{0} : absence, \code{1} : presence, \code{NA} : indeterminate) 
##'  for a single species that will be used to 
##' evaluate the species distribution model(s) with independent data
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##'  \code{SpatialPoints}  (if presence-only) or \code{SpatialPointsDataFrame}
##'  object containing binary data.}
##' @param eval.expl.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables (in 
##' columns or layers) that will be used to evaluate the species distribution model(s) with 
##' independent data.
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' @param eval.resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to evaluate 
##' the species distribution model(s) with independent data
##' 
##' @param params.PA a \code{list} with the species names associated to the parameters of PA
##' 
##' @param single.formated.data a \code{list} of \code{BIOMOD.formated.data} objects to add to the MS.formated.data object
##' 
##' @param filter.raster (\emph{optional, default} \code{FALSE}) \cr 
##' If \code{expl.var} is of raster type, a \code{logical} value defining whether \code{resp.var} 
##' is to be filtered when several points occur in the same raster cell
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' 
##' 
##' 
##' @return 
##' 
##' A \code{\link{MS.formated.data}} object that can be used to build species distribution 
##' model(s) with the \code{\link{MS_Modeling}} function. \cr
##' 
##' 
##' @keywords dataset format evaluation pseudo-absence
##' 
##' 
##' @seealso \code{\link{bm_PseudoAbsences}}, \code{\link{BIOMOD_Modeling}}
##' @family Main functions
##' 
##' @examples
##' library(terra)
##' library(biomod2)
##' 
##' # Load species occurrences (6 species available)
##' data(DataSpecies)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- c("PantheraOnca", "PteropusGiganteus")
##' 
##' # Get corresponding presence/absence data
##' myResp <- DataSpecies[, myRespName]
##' 
##' # Get corresponding XY coordinates
##' myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
##' 
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' data(bioclim_current)
##' myExpl <- terra::rast(bioclim_current)
##' 
##' \dontshow{
##'   myExtent <- terra::ext(0,30,45,70)
##'   myExpl <- terra::crop(myExpl, myExtent)
##' }
##' 
##' 
##' 
##' myMSData <- MS_FormatingData(ms.project.name = "Example_MS",
##'                              resp.name = myRespName,
##'                              resp.var = myResp,
##'                              expl.var = myExpl,
##'                              data.type = "binary",
##'                              resp.xy = myRespXY)
##' 
##' myMSData
##' summary(myMSData)
##' plot(myMSData)
##' 
##' \dontshow{
##'   unlink('Example_MS', recursive = TRUE)
##' }
##' 
##' 
##' 
##' 
##' @importFrom biomod2 BIOMOD_FormatingData
##' @importFrom foreach foreach %do%
##' @importFrom methods .hasSlot
##' @importFrom utils packageVersion
##' 
##' @export
##' 
##' 
###################################################################################################

MS_FormatingData <- function(ms.project.name,
                             dir.name = '.',
                             resp.name = NULL,
                             resp.var = NULL,
                             expl.var = NULL,
                             data.type = NULL,
                             resp.xy = NULL,
                             eval.resp.var = NULL,
                             eval.expl.var = NULL,
                             eval.resp.xy = NULL,
                             params.PA = NULL, 
                             single.formated.data = NULL,
                             filter.raster = FALSE,
                             seed.val = NULL)
{ 
  .bm_cat("MultiSpecies Data Formating")
  
  ## 1. check args ------------------------------------------------------------
  args <- .MS_FormatingData.check.args(resp.name, dir.name, data.type, resp.xy, expl.var, single.formated.data)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 2. Create folder
  dir.create(file.path(dir.name, ms.project.name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, ms.project.name, ".BIOMOD_DATA", "single.formated.data"), showWarnings = FALSE, recursive = TRUE) #ça va être encore super long
  dir.create(file.path(dir.name, ms.project.name, ".BIOMOD_DATA", "output"), showWarnings = FALSE, recursive = TRUE)
  file.txt <- file.path(dir.name, ms.project.name, ".BIOMOD_DATA", "output", "MS_FormatingData.output.txt")
  nameFolder <- file.path(dir.name, ms.project.name, ".BIOMOD_DATA", "single.formated.data")
  
  cat("Creation of ms.formated.data \n\n", file = file.txt, append = FALSE)
  
  ## 3.Create MS.formated.data with the data frame
  new.formated.data <- foreach(sp = resp.name) %do% { ## No parallelisation as it is really quick but can be add if necessary ? 

    parameters <- params.PA[[sp]]
    
    capture.output(formated.data <- BIOMOD_FormatingData(dir.name = file.path(dir.name, ms.project.name),
                                                         resp.var = resp.var[,sp],
                                                         expl.var = expl.var,
                                                         resp.xy = resp.xy,
                                                         resp.name = sp,
                                                         eval.resp.var = eval.resp.var[,sp],
                                                         eval.expl.var = eval.expl.var,
                                                         eval.resp.xy = eval.resp.xy,
                                                         PA.nb.rep = parameters$PA.nb.rep,
                                                         PA.nb.absences = parameters$PA.nb.absences,
                                                         PA.strategy = parameters$PA.strategy,
                                                         PA.dist.min = parameters$PA.dist.min,
                                                         PA.dist.max = parameters$PA.dist.max,
                                                         PA.sre.quant = parameters$PA.sre.quant,
                                                         PA.fact.aggr = parameters$PA.fact.aggr,
                                                         PA.user.table = parameters$PA.user.table,
                                                         na.rm = T,
                                                         filter.raster = filter.raster),
                   file = file.txt, append = TRUE)
    save(formated.data, file = file.path(nameFolder, paste0(sp,".sfd")), compress = TRUE)
    #cat(output, fill = TRUE, file = file.txt, append = TRUE)
    cat("\n\n", file = file.txt, append = TRUE)
    return(formated.data)
  }
  
  names(new.formated.data) <- resp.name
  
  single.formated.data <- c(new.formated.data, single.formated.data)
  
  ## 2.Create MS.formated.data with the first single.formated.data 
  if (!is.null(single.formated.data)){
    sfd <- single.formated.data[[1]]
    single.formated.data <- single.formated.data[-1]
    
    cat("\n\t Data formating of", sfd@sp.name)
    sfd@dir.name = file.path(dir.name, ms.project.name)
    save(sfd, file = file.path(nameFolder, paste0(sfd@sp.name,".sfd")), compress = TRUE)
    
    PA.strategy <- ifelse(.hasSlot(sfd, "PA.strategy"), sfd@PA.strategy, "real")
    PA.table <- list()
    if (PA.strategy != "real"){
      PA.table[[sfd@sp.name]] <- sfd@PA.table
    } 
    
    eval.sp <- list()
    if(sfd@has.data.eval){
      eval.sp[[sfd@sp.name]] <- sfd@eval.data.species
    }
    
    MSFD <- new(
      'MS.formated.data',
      data.type = sfd@data.type,
      dir.name = R.utils::getAbsolutePath(dir.name),
      ms.project = ms.project.name,
      sp.name = c(sfd@sp.name),
      coord = sfd@coord,
      data.species = list(sfd@data.species),
      data.env.var = sfd@data.env.var,
      data.mask = sfd@data.mask,
      PA.strategy = PA.strategy,
      PA.table = PA.table,
      has.data.eval = sfd@has.data.eval,
      eval.coord = sfd@eval.coord,
      eval.data.species = eval.sp,
      eval.data.env.var = sfd@eval.data.env.var,
      has.filter.raster = filter.raster,
      biomod3.version = as.character(packageVersion("biomod2")) #change to biomod3 after
    )
  }
  ## 3.Add single.formated.data to MS.formated.data
  
  if (!is.null(single.formated.data) && length(single.formated.data) != 0){
    for (i in 1:length(single.formated.data)){
      sfd <- single.formated.data[[i]]
      sfd@dir.name = file.path(dir.name, ms.project.name)
      MSFD <- MS_AddData(MSFD, sfd)
      save(sfd, file = file.path(nameFolder, paste0(sfd@sp.name,".sfd")), compress = TRUE)
    }
  }
  
  names(MSFD@data.species) <- MSFD@sp.name
  #names(MSFD@PA.table) <- MSFD@sp.name
  
  .bm_cat("Done")
  return(MSFD)
}


# Argument Check -------------------------------------------------------------

.MS_FormatingData.check.args <- function(resp.name, dir.name, data.type, resp.xy, expl.var, single.formated.data)
{
  ## 0. Checking names (resp.name available ?) --------------------------------
  if (any(grepl('/', resp.name))) {
    stop(paste0("Response variable name must be a character, and not a file path."
                , "\n Please refer to dir.name parameter to set a modeling folder."))
  }
  if (any(grepl('_', resp.name) | grepl(' ', resp.name))) {
    resp.name <- sapply(resp.name, function(x) {paste(unlist(strsplit(x, '_')), collapse = '.')})
    resp.name <- sapply(resp.name, function(x) {paste(unlist(strsplit(x, ' ')), collapse = '.')})
    cat('\n      ! Response variable name was converted into', resp.name)
  }
  if (!dir.exists(dir.name)) {
    stop(paste0("Modeling folder '", dir.name, "' does not exist"))
  }
  
  if (length(single.formated.data) != 0){
    dt <- ifelse(is.null(data.type), single.formated.data[[1]]@data.type, data.type)
    #xy <- ifelse(is.null(resp.xy), single.formated.data[[1]]@coord, resp.xy)
    
    for (i in 1:length(single.formated.data)){
      bm.format <- single.formated.data[[i]]
      if (dt != bm.format@data.type){
        stop("All your species must have the same datatype")
      }

      # if (!identical(MS.format@coord, bm.format@coord) | 
      #     !identical(MS.format@data.env.var,bm.format@data.env.var)){
      #   stop("Your environnement and coord must be the same")
      # }
    }
  }
  
  
  return(list(resp.name = resp.name, dir.name = dir.name))
}


# Add BIOMOD.formated.data to a MS.formated.data ------------------------------

##' 
##' @importFrom methods .hasSlot
##' 

MS_AddData <- function(MS.format, bm.format){
  # Ajouter check si on ouvre cette fonction
  cat("\n\t Data formating of", bm.format@sp.name)
  
  if (MS.format@data.type != bm.format@data.type){
    stop("All your species must have the same datatype")
  }
  #dir.name = R.utils::getAbsolutePath(dir.name) dir.name de MS prend le pas
  MS.format@sp.name <- c(MS.format@sp.name, bm.format@sp.name)
  if (!identical(MS.format@coord, bm.format@coord) | 
      !identical(MS.format@data.env.var,bm.format@data.env.var) | 
      !identical(MS.format@data.mask, bm.format@data.mask)){
    stop("Your environnement and coord must be the same")
  }
  MS.format@data.species[[bm.format@sp.name]] <- bm.format@data.species
  
  if (.hasSlot(bm.format, "PA.strategy")) {
    MS.format@PA.strategy <- c(MS.format@PA.strategy, bm.format@PA.strategy)
    MS.format@PA.table[[bm.format@sp.name]] <- bm.format@PA.table
  }
  
  
  if (MS.format@has.data.eval){
    if (!identical(MS.format@eval.coord, bm.format@eval.coord) | 
        !identical(MS.format@eval.data.env.var,bm.format@eval.data.env.var) ){
      stop("Your eval environnement and eval coord must be the same")
    }
    MS.format@eval.data.species[[bm.format@sp.name]] <- bm.format@eval.data.species
  }
  
  return(MS.format)
}
