\name{gstatModel-class}
\docType{class}
\alias{gstatModel-class}
\alias{predict,gstatModel}
\title{A class for a linear geostatistical model}
\description{A class containing fitted parameters of a linear geostatistical model. It contains regression (GLM) and the kriging part (variogram).}
\section{Slots}{
  \describe{
    \item{\code{regModel}:}{object of class \code{"glm"}; output of fitting a generalized linear model (GLM)}
    \item{\code{vgmModel}:}{object of class \code{"data.frame"}; the fitted gstat variogram model parameters containing }
    \item{\code{sp}:}{object of class \code{"SpatialPointsDataFrame"}; observation locations}
  }
}
\section{Methods}{
  \describe{
	\item{summary}{\code{signature(obj = "gstatModel")}: summarize object by showing the mapping accuracy (cross-validation) and the amount of variation explained by the model }
	\item{predict}{\code{signature(obj = "gstatModel")}: makes predictions for a set of given predictionLocations (gridded maps) }
  }
}
\note{ \code{"SpatialPredictions"} saves results of predictions for a single target variable, which can be of type numeric or factor. Multiple variables can be combined into a list. When using \code{nsim} argument with the \code{predict} method, the output result will be of type:\cr

 \code{plotKML::RasterBrickSimulations-class} \cr

i.e. \emph{N} number of equiprobable realizations. To generate an object of type:\cr 

\code{\link{SpatialPredictions-class}}\cr

use \code{nsim = 0}. 
}
\author{ Tomislav Hengl }
\seealso{ \code{plotKML::SpatialPredictions-class}, \code{plotKML::RasterBrickSimulations-class}, \code{gstat::gstat}, \code{stats::glm}
}
\examples{ 
# load observations:
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
# load grids:
data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
\dontrun{# fit a model:
omm <- fit.gstatModel(meuse, formulaString = om~dist, covariates = meuse.grid)
show(omm@regModel)
# produce SpatialPredictions:
om.rk <- predict(omm, predictionLocations = meuse.grid)
# overview of the prediction process:
show(om.rk)
# plot the whole project and open in Google Earth:
library(plotKML)
plotKML(om.rk)
}
}
\keyword{classes}
