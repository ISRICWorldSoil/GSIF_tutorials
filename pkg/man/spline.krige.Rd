\name{spline.krige}
\alias{spline.krige}
\alias{resample.grid}
\title{Kriging combined with splines}
\description{Combines kriging and spline interpolation to speed up kriging. Spline interpolation is implemented via the SAGA GIS function \code{"Cubic Spline Approximation"}.}
\usage{spline.krige(formula, locations, newdata, newlocs = NULL, model, 
    te = as.vector(newdata@bbox), file.name, silent = FALSE, 
    t_cellsize = newdata@grid@cellsize[1], optN = 20, quant.nndist = .5, 
    nmax = 30, \dots)
}
\arguments{
  \item{formula}{formula that defines the dependent variable as a linear model of independent variables; usually in the form \code{z~1}}
  \item{locations}{object of class \code{SpatialPoints}; sampling locations}
  \item{newdata}{object of class \code{SpatialPixels*}; spatial domain of interest}
  \item{newlocs}{object of class \code{SpatialPoints*}; prediction locations produced using the \code{resample.grid} function (if missing it will be generated using the \code{resample.grid} function)}
  \item{model}{variogram model of dependent variable (or its residuals); see \code{gstat::krige}}
  \item{te}{numeric; a vector in the form c(xmin,ymin,xmax,ymax); sets bounding box of the kriging predictions}
  \item{file.name}{character; optional output file name pattern}
  \item{silent}{logical; specifies whether to print out the progress}
  \item{t_cellsize}{numeric; target cell size (output grid)}
  \item{optN}{integer; optimal number of prediction locations per sampling location}    
  \item{quant.nndist}{numeric; threshold probability to determine the shortest distance}
  \item{nmax}{integer; the number of nearest observations that
  should be used for kriging}
  \item{\dots}{other optional arguments that can be passed to function \code{gstat::krige}}
}
\value{
Returns an object of class \code{"SpatialGridDataFrame"}, or output file name.
}
\note{
The function adjusts grid density (prediction locations) to the actual sampling intensity and then fills in prediction gaps using spline interpolation. This allows for speeding up the kriging process in the order of magnitude times without loosing precision. 
Spline interpolation is implemented via the SAGA GIS v2.1 function \code{"Cubic Spline Approximation"} using the default settings. This function is especially suitable for producing predictions for large grids where the sampling locations show high clustering.
}
\author{ Tomislav Hengl }
\examples{
library(plotKML)
library(spatstat)
library(RSAGA)
library(gstat)
library(raster)
data(eberg)
data(eberg_grid25)
library(sp)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
m <- vgm(psill=320, model="Exp", range=1200, nugget=160)
plot(variogram(SNDMHT_A~1, eberg[!is.na(eberg$SNDMHT_A),]), m)
## prediction locations:
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")
grd <- resample.grid(locations=eberg["SNDMHT_A"], newdata=eberg_grid25, 
   t_cellsize=25, optN=5, quant.nndist=.5, nstrata=4)
## plot resampled grid:
plot(grd)
points(eberg, pch=19, col="red", cex=.7)
env <- rsaga.env()
if(env$version=="2.1.0"){
 ## ordinary kriging:
 res.ok <- spline.krige(locations=eberg["SNDMHT_A"], t_cellsize=10, 
      newdata=eberg_grid25, newlocs=grd, model=m)
 plot(raster(res.ok[1]))
}
}