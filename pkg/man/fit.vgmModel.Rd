\name{fit.vgmModel-methods}
\alias{fit.vgmModel}
\alias{fit.vgmModel,formula,data.frame,SpatialPixelsDataFrame-method}
\title{Fits a 2D or 3D variogram model to spatial data}
\description{Fits a 2D or 3D variogram model based on a regression matrix and spatial domain of interest.}
\usage{
\S4method{fit.vgmModel}{formula,data.frame,SpatialPixelsDataFrame}(formulaString,
      rmatrix, predictionDomain, vgmFun = "Exp", 
      dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]],
      anis = NULL, subsample = nrow(rmatrix), ivgm, \dots)
}
\arguments{
  \item{formulaString}{object of class \code{"formula"} --- regression model} 
  \item{rmatrix}{object of class \code{"data.frame"}; regression matrix produced as a result of spatial overlay}
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; spatial domain of interest}
  \item{vgmFun}{character; variogram function (\code{"Exp"} by default)}
  \item{dimensions}{character; \code{"3D"}, \code{"2D"}, \code{"2D+T"}, \code{"3D+T"} models}
  \item{anis}{vector containing 2, 5 or more anisotropy parameters; see \code{gstat::vgm} for more info}
  \item{subsample}{integer; size of the subset}
  \item{ivgm}{vgm; initial variogram model}
  \item{\dots}{other optional arguments that can be passed to \code{gstat::fit.variogram}}
}
\details{It will try to fit a variogram to multidimensional data. If the data set is large, this process can be time-consuming, hence one way to speed up fitting is to subset the regression matrix using the \code{subsample} argument (random subsetting).}
\author{ Tomislav Hengl }
\seealso{ \code{\link{fit.regModel}}, \code{\link{fit.gstatModel}}, \code{gstat::fit.variogram} }
\examples{
library(sp)

# fit variogram to the Meuse data:
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
# produce a regression matrix:
index <- overlay(meuse.grid, meuse)
sel <- !is.na(index)
ov <- data.frame(meuse.grid[index[sel],])
ov$om <- meuse@data[sel,"om"]
# fit a model:
v <- fit.vgmModel(om~1, rmatrix=ov, meuse.grid, dimensions="2D")
plot(variogram(om ~ 1, meuse[!is.na(meuse$om),]), v$vgm)
}
