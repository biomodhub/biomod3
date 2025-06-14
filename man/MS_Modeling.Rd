% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MS_Modeling.R
\name{MS_Modeling}
\alias{MS_Modeling}
\title{Run a range of species distribution models}
\usage{
MS_Modeling(
  ms.format,
  modeling.id = as.character(format(Sys.time(), "\%s")),
  models = c("ANN", "DNN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF",
    "RFd", "SRE", "XGBOOST"),
  params.CV = NULL,
  params.OPT = NULL,
  weights = NULL,
  prevalence = NULL,
  metric.eval = c("KAPPA", "TSS", "AUCroc"),
  var.import = 0,
  scale.models = FALSE,
  nb.cpu = 1,
  seed.val = NULL
)
}
\arguments{
\item{ms.format}{a \code{\link{MS.formated.data}} object returned by the \code{\link{MS_FormatingData}} function}

\item{modeling.id}{a \code{character} corresponding to the name (ID) of the simulation set
(\emph{a random number by default})}

\item{models}{a \code{vector} containing model names to be computed, must be among
\code{ANN}, \code{CTA}, \code{DNN}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS},
\code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}}

\item{params.CV}{a \code{list} with the species names associated to the parameters of Cross-Validation. See BIOMOD_Modeling}

\item{params.OPT}{a \code{list} with the species names associated to the options of the algorithms. See BIOMOD_Modeling}

\item{weights}{(\emph{optional, default} \code{NULL}) \cr
A \code{vector} of \code{numeric} values corresponding to observation weights (one per
observation, see Details)}

\item{prevalence}{(\emph{optional, default} \code{NULL}) \cr
A \code{numeric} between \code{0} and \code{1} corresponding to the species prevalence to
build '\emph{weighted response weights}' (see Details)}

\item{metric.eval}{a \code{vector} containing evaluation metric names to be used, must
be among \code{AUCroc}, \code{AUCrpg}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD},
\code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{OR},
\code{ORSS}, \code{BOYCE}, \code{MPA}, \code{RMSE}, \code{MAE}, \code{MSE}, \code{Rsquared}, \code{Rsquared_aj},
\code{Max_error}, \code{Accuracy}, \code{"Recall"}, \code{"Precision"}, \code{"F1"}}

\item{var.import}{(\emph{optional, default} \code{NULL}) \cr
An \code{integer} corresponding to the number of permutations to be done for each variable to
estimate variable importance}

\item{scale.models}{(\emph{optional, default} \code{FALSE}) \cr
A \code{logical} value defining whether all models predictions should be scaled with a
binomial GLM or not}

\item{nb.cpu}{(\emph{optional, default} \code{1}) \cr
An \code{integer} value corresponding to the number of computing resources to be used to
parallelize the single models computation}

\item{seed.val}{(\emph{optional, default} \code{NULL}) \cr
An \code{integer} value corresponding to the new seed value to be set}
}
\value{
A \code{\link{MS.models.out}} object acting as a proxi for the created \code{BIOMOD.models.out} objects.
}
\description{
This function allows to calibrate and evaluate a range of modeling techniques
for some given species distribution. It run the function BIOMOD_Modeling for all the different species..
}
\examples{
library(terra)
library(biomod2)

# Load species occurrences (6 species available)
data(DataSpecies)

# Select the name of the studied species
myRespName <- c("PantheraOnca", "PteropusGiganteus")

# Get corresponding presence/absence data
myResp <- DataSpecies[, myRespName]

# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data(bioclim_current)
myExpl <- terra::rast(bioclim_current)




myMSData <- MS_FormatingData(ms.project.name = "Example_MS",
                             resp.name = myRespName,
                             resp.var = myResp,
                             expl.var = myExpl,
                             data.type = "binary",
                             resp.xy = myRespXY)

params.CV <- list("PantheraOnca" = list(CV.strategy = "random", CV.nb.rep = 2, CV.perc = 0.8),
                  "PteropusGiganteus" = list(CV.strategy = "random", CV.nb.rep = 2, CV.perc = 0.8))

params.OPT <- list("PantheraOnca" = list(OPT.strategy = "bigboss"),
                   "PteropusGiganteus" = list(OPT.strategy = "bigboss"))

myMSModelOut <- MS_Modeling(myMSData, 
                            modeling.id = "FirstModels",
                            models = c("GLM", "XGBOOST"),
                            params.CV = params.CV,
                            params.OPT = params.OPT)

myMSModelOut

get_evaluations(myMSModelOut, sp = "PantheraOnca")


\dontshow{
  unlink('Example_MS', recursive = TRUE)
}

}
\author{
Helene Blancheteau
}
