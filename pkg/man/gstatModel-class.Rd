\name{gstatModel-class}
\docType{class}
\alias{gstatModel-class}
\alias{predict,gstatModel-method}
\alias{validate,gstatModel-method}
\alias{validate}
\alias{predict.gstatModel}
\alias{predict,list-method}
\alias{predict.gstatModelList}
\title{A class for a linear geostatistical model}
\description{A class containing fitted parameters of a linear geostatistical model. It comprises regression model (GLM), variogram model, and observation locations of sampled values used to fit the model.}
\section{Slots}{
  \describe{
    \item{\code{regModel}:}{object of class \code{"glm"}; output of fitting a generalized linear model (GLM)}
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
\note{ \code{"SpatialPredictions"} saves results of predictions for a single target variable, which can be of type numeric or factor. Multiple variables can be combined into a list. When using \code{nsim} argument with the \code{predict} method, the output result will be of type:\cr

 \code{plotKML::RasterBrickSimulations-class} \cr

i.e. \emph{N} number of equiprobable realizations. To generate an object of type:\cr 

\code{plotKML::SpatialPredictions-class}\cr

set \code{nsim = 0}. 
}
\author{ Tomislav Hengl and Gerard B.M. Heuvelink }
\seealso{ \code{plotKML::SpatialPredictions-class}, \code{plotKML::RasterBrickSimulations-class}, \code{gstat::gstat}, \code{stats::glm}
}
\examples{ 
# load observations:
library(plotKML)
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
# load grids:
data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
# fit a model:
omm <- fit.gstatModel(meuse, om~dist+ffreq, family=gaussian(link="log"), meuse.grid)
show(omm@regModel)
# produce SpatialPredictions:
om.rk <- predict(omm, predictionLocations = meuse.grid)
# plotKML(om.rk)
# run a proper cross-validation:
rk.cv <- validate(omm)
}
\keyword{classes}
