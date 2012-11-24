\name{fit.regModel-methods}
\alias{fit.regModel}
\alias{fit.regModel,formula,data.frame,SpatialPixelsDataFrame,character-method}
\title{Fits a regression model to spatial data}
\description{Fits a regression model (e.g. GLM) and variogram for the residuals using the default settings.}
\usage{
\S4method{fit.regModel}{formula,data.frame,SpatialPixelsDataFrame,character}(formulaString,
      rmatrix, predictionDomain, method = list("GLM", "cubist", "HB")[[1]], 
      dimensions = NULL, family=gaussian, stepwise=TRUE, rvgm=NULL, 
      vgmFun="Exp", subsample = 5000, \dots)
}
\arguments{
  \item{formulaString}{object of class \code{"formula"} --- regression model} 
  \item{rmatrix}{object of class \code{"data.frame"}; regression matrix produced as a result of spatial overlay}
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; spatial domain of interest}
  \item{method}{character; family of methods considered e.g. \code{"GLM"}}
  \item{dimensions}{character; \code{"3D"}, \code{"2D"}, \code{"2D+T"}, \code{"3D+T"} models}
  \item{family}{family to be passed to the \code{glm} (see examples below)}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{rvgm}{residual variogram (assumed to be unknown)}
  \item{vgmFun}{character; variogram function (\code{"Exp"} by default)}
  \item{subsample}{integer; size of the subset}
  \item{\dots}{other optional arguments that can be passed to \code{glm} and/or \code{gstat::fit.variogram}}
}
\details{Produces an object of class \code{"gstatModel"} that contains: (1) fitted GLM, (2) fitted variogram, and (c) object of class \code{"SpatialPoints"} with observation locations. To combine overlay and model fitting operations, use \code{\link{fit.gstatModel}}.}
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
m <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, method="GLM", family=gaussian(log))
m@regModel
m@vgmModel
}
