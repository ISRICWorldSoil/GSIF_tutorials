\name{glm.sp}
\alias{glm.sp}
\title{Fits a GLM to spatial data}
\description{Fits a GLM and variogram for the residuals using the default settings used in the gstat and glm packages.}
\usage{
glm.sp(formulaString, rmatrix, predictionDomain, family=gaussian, stepwise=TRUE, rvgm=NULL, vgmFun="Exp", type, \dots)
}
\arguments{
  \item{formulaString}{object of class \code{"formula"} --- regression model} 
  \item{rmatrix}{object of class \code{"data.frame"}; regression matrix produced as a result of spatial overlay}
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; spatial domain of interest}
  \item{family}{family to be passed to the \code{glm} (see examples below)}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{rvgm}{residual variogram (assumed to be unknown)}
  \item{vgmFun}{character; variogram function (\code{"Exp"} by default)}
  \item{type}{character; spatial object type (e.g. \code{"SpatialPointsDataFrame"} or \code{"geosamples"})}
  \item{\dots}{other optional arguments that can be passed to \code{glm} and/or \code{fit.variogram}}
}
\details{Produces a list of three objects: (1) fitted GLM, (2) fitted variogram, and (c) object of class \code{"SpatialPoints"} with observation locations. To combine overlay and model fitting operations, use \code{\link{fit.gstatModel}}.}
\author{ Tomislav Hengl }
\seealso{ \code{\link{fit.gstatModel}}, \code{stats::glm}, \code{gstat::fit.variogram} }
\examples{
# universal kriging with Meuse data:
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
# produce a regression matrix:
index <- overlay(meuse.grid, meuse)
sel <- !is.na(index)
ov <- cbind(data.frame(meuse[sel,"om"]), meuse.grid[index[sel],])
# fit a model:
m <- glm.sp(om~dist, rmatrix=ov, meuse.grid, family=gaussian(log), type="SpatialPointsDataFrame")
m$regModel
m$vgmModel
}
