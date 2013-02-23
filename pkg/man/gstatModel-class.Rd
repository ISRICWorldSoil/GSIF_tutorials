\name{gstatModel-class}
\docType{class}
\alias{gstatModel-class}
\alias{validate,gstatModel-method}
\alias{validate}
\title{A class for a geostatistical model}
\description{A class containing fitted parameters of a geostatistical model to be used to run predictions by regression-kriging. It comprises regression model (e.g. GLM), variogram model, and observation locations of sampled values used to fit the model.}
\section{Slots}{
  \describe{
    \item{\code{regModel}:}{object of class \code{"ANY"}; output of fitting a generalized linear model (GLM) or any similar regression model}
    \item{\code{vgmModel}:}{object of class \code{"data.frame"}; the fitted gstat variogram model parameters containing variogram model, nugget, sill, range and the five anisotropy parameters}
    \item{\code{sp}:}{object of class \code{"SpatialPoints"}; observation locations}
  }
}
\section{Methods}{
  \describe{
	\item{summary}{\code{signature(obj = "gstatModel")}: summarize object by showing the mapping accuracy (cross-validation) and the amount of variation explained by the model }
	\item{predict}{\code{signature(obj = "gstatModel")}: makes predictions for a set of given predictionLocations (gridded maps) at block support corresponding to the \code{cellsize} slot in the object of class \code{"SpatialPixelsDataFrame"}; to produce predictions at point support, submit the \code{predictionLocations} as \code{"SpatialPointsDataFrame"} }
	\item{validate}{\code{signature(obj = "gstatModel")}: runs \emph{n}-fold cross-validation of the existing gstatModel (it re-fits the model using existing formula string and then validates the predictions at validation locations)}
  }
}
\details{ Any model passed to the \code{regModel} slot must come with generic functions such as \code{residuals}, \code{fitted.values}, \code{summary} and \code{predict}.
}
\note{ \code{"SpatialPredictions"} saves results of predictions for a single target variable, which can be of type numeric or factor. Multiple variables can be combined into a list. When using \code{nsim} argument with the \code{predict} method, the output result will be of type:\cr

 \code{plotKML::RasterBrickSimulations-class} \cr

i.e. \emph{N} number of equiprobable realizations. To generate an object of type:\cr 

\code{plotKML::SpatialPredictions-class}\cr

set \code{nsim = 0}. 
}
\author{ Tomislav Hengl and Gerard B.M. Heuvelink }
\seealso{ \code{\link{predict.gstatModel}}, \code{plotKML::SpatialPredictions-class}, \code{plotKML::RasterBrickSimulations-class}, \code{gstat::gstat}, \code{stats::glm}
}
\examples{ 
## load observations:
library(plotKML)
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
## load grids:
data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
## fit a model:
omm <- fit.gstatModel(meuse, om~dist+ffreq, family=gaussian(link="log"), meuse.grid)
show(omm@regModel)
## produce SpatialPredictions:
om.rk <- predict(omm, predictionLocations = meuse.grid)
#plotKML(om.rk)
## run a proper cross-validation:
rk.cv <- validate(omm)
## RMSE:
sqrt(mean((rk.cv$validation$var1.pred-rk.cv$validation$observed)^2))
}
\keyword{classes}
