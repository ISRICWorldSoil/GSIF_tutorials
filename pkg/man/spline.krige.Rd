\name{spline.krige}
\alias{spline.krige}
\alias{resample.grid}
\title{Kriging combined with splines}
\description{Combines kriging and spline interpolation to speed up the kriging with minimal loss in precision, whilst reducing generation of artifacts. Spline interpolation is implemented via the SAGA GIS function \code{"Multilevel B-Spline Interpolation"} (SAGA GIS needs to be installed separately).}
\usage{spline.krige(formula, locations, newdata, newlocs = NULL, model, 
    te = as.vector(newdata@bbox), file.name, silent = FALSE, 
    t_cellsize = newdata@grid@cellsize[1], optN = 20, quant.nndist = .5, 
    nmax = 30, predictOnly = FALSE, resample = TRUE, saga.env, 
    saga.lib=c("grid_spline","grid_tools"), saga.module=c(4,0), \dots)
}
\arguments{
  \item{formula}{formula that defines the dependent variable as a linear model of independent variables; usually in the form \code{z~1}}
  \item{locations}{object of class \code{SpatialPoints}; sampling locations}
  \item{newdata}{object of class \code{SpatialPixels*}; spatial domain of interest}
  \item{newlocs}{object of class \code{SpatialPoints*}; prediction locations produced using the \code{resample.grid} function (if missing it will be generated using the \code{resample.grid} function)}
  \item{model}{variogram model of dependent variable (or its residuals); see \code{gstat::krige}}
  \item{te}{numeric; a vector in the form \code{c(xmin,ymin,xmax,ymax)}; sets bounding box of the kriging predictions}
  \item{file.name}{character; optional output file name pattern (without any file extension)}
  \item{silent}{logical; specifies whether to print out the progress}
  \item{t_cellsize}{numeric; target cell size (output grid)}
  \item{optN}{integer; optimal number of prediction locations per sampling location e.g. 1 sampling location is used to predict values for 20 new pixels}    
  \item{quant.nndist}{numeric; threshold probability to determine the search radius (sigma)}
  \item{nmax}{integer; the number of nearest observations that
  should be used for kriging}
  \item{predictOnly}{logical; specifies whether to generate only predictions (\code{var1.pred} column)}
  \item{resample}{logical; specifies whether to down or upscale SAGA GIS grids to match the grid system of \code{newdata}}
  \item{saga.env}{list; path to location of the SAGA binaries (extracted using \code{rsaga.env()})}
  \item{saga.lib}{character; names of the SAGA libraries used}
  \item{saga.module}{integer; corresponding module numbers}
  \item{\dots}{other optional arguments that can be passed to function \code{gstat::krige}}
}
\value{
Returns an object of class \code{"SpatialGridDataFrame"}, or an output file name.
}
\note{
This function adjusts grid density (prediction locations) in reference to the actual local sampling intensity. High resolution grids are created where sampling density is higher and vice versa (\href{http://dx.doi.org/10.1016/j.cageo.2005.11.008}{Hengl, 2006}). Low resolution grids (due to sparse data) are then downscaled to the target resolution using spline interpolation. This allows for speeding up the kriging with minimal loss in precision, whilst reducing generation of artifacts. 
Spline interpolation is implemented via the SAGA GIS v2.1 function \code{"Multilevel B-Spline Interpolation"} using the default settings. This function is especially suitable for producing predictions for large grids where the sampling locations show high spatial clustering. It is NOT intended for predicting using point samples collected using sampling designs with constant spatial sampling intensity e.g. point samples collected using simple random sampling or grid sampling.
}
\references{
\itemize{
\item Hengl T., (2006) \href{http://dx.doi.org/10.1016/j.cageo.2005.11.008}{Finding the right pixel size}. Computers and Geosciences, 32(9): 1283-1298.
\item SAGA GIS (\url{http://saga-gis.org})
\item SpatStat package (\url{http://www.spatstat.org})
}
}
\author{ Tomislav Hengl }
\examples{
\dontrun{
library(plotKML)
library(spatstat)
library(RSAGA)
library(gstat)
library(raster)
data(eberg)
data(eberg_grid)
data(eberg_grid25)
library(sp)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
m <- vgm(psill=320, model="Exp", range=1200, nugget=160)
plot(variogram(SNDMHT_A~1, eberg[!is.na(eberg$SNDMHT_A),]), m)
## prediction locations:
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")
## prepare prediction locations for spline.krige:
grd <- resample.grid(locations=eberg["SNDMHT_A"], t_cellsize=25,
   newdata=eberg_grid25, optN=5, quant.nndist=.9)
## plot resampled grid:
plot(raster(grd$density))
plot(grd$newlocs)
points(eberg, pch=19, col="red", cex=.7)
env <- rsaga.env()
if(exists("env") & env$version=="2.1.0"){
 ## compare processing time:
 system.time( SND.sok <- spline.krige(locations=eberg["SNDMHT_A"], 
      t_cellsize=25, newdata=eberg_grid25, 
      newlocs=grd$newlocs, model=m, nmax=30) )
 system.time( SND.ok <- krige(SNDMHT_A~1, 
      eberg[!is.na(eberg$SNDMHT_A),], 
      newdata=eberg_grid, m, 
      debug.level = -1, nmax=30) )
 system.time( SND.ok25 <- krige(SNDMHT_A~1, 
      eberg[!is.na(eberg$SNDMHT_A),], 
      newdata=eberg_grid25, m, 
      debug.level = -1, nmax=30) )  
 ## compare outputs visually:
 par(mfrow=c(1,3))
 plot(raster(SND.sok[1]), main="spline.krige (25 m)")
 plot(raster(SND.ok25[1]), main="krige (25 m)")
 plot(raster(SND.ok[1]), main="krige (100 m)") 
}
}
## conclusion: spline.krige produces less artifacts, 
## and is at order of magnitude faster than simple 'krige'
}