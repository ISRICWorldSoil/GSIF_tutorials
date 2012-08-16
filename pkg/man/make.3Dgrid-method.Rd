\name{make.3Dgrid}
\docType{methods}
\encoding{latin1}
\alias{make.3Dgrid}
\alias{sp3D}
\alias{sp3D,SpatialPixelsDataFrame-method}
\alias{sp3D,list-method}
\alias{gdalwarp}
\alias{gdalwarp,SpatialPixelsDataFrame-method}
\alias{make.3Dgrid,SpatialPixelsDataFrame-method}
\alias{make.3Dgrid,RasterBrick-method}
\title{Methods to prepare 3D prediction locations}
\description{Generates a list of objects of type \code{"SpatialPixelsDataFrame"} with \code{longitude}, \code{latitude} and \code{altitude} coordinates (these names are used by default for compatibility with the \code{\link{geosamples-class}}).}
\usage{
\S4method{make.3Dgrid}{SpatialPixelsDataFrame}(obj, 
      proj4s = get("ref_CRS", envir = GSIF.opts), 
      pixsize = get("cellsize", envir = GSIF.opts)[2], 
      resampling_method = "bilinear", 
      NAflag = get("NAflag", envir = GSIF.opts), 
      stdepths = get("stdepths", envir = GSIF.opts), 
      tmp.file = TRUE, show.output.on.console = FALSE, \dots)
\S4method{make.3Dgrid}{RasterBrick}(obj, 
      proj4s = get("ref_CRS", envir = GSIF.opts), 
      pixsize = get("cellsize", envir = GSIF.opts)[2], 
      resampling_method = "bilinear", 
      NAflag = get("NAflag", envir = GSIF.opts), 
      stdepths = get("stdepths", envir = GSIF.opts), 
      tmp.file = TRUE, show.output.on.console = FALSE, \dots)
}
\arguments{
  \item{obj}{object of class \code{"SpatialPixelsDataFrame"} or \code{"RasterBrick"}} 
  \item{proj4s}{character; proj4string describing the target coordinate system}
  \item{pixsize}{grid cell size in decimal degrees (set by default at 1/1200 (0.0008333333 or 100 m around equator)} 
  \item{resampling_method}{character; resampling method to be passed the reprojection algorithm}
  \item{NAflag}{character; missing value flag}
  \item{stdepths}{numeric; list of standard depths}
  \item{tmp.file}{logical; specifies whether a temporary file name should be generated}
  \item{show.output.on.console}{logical; specifies whether to print out the progress}
  \item{\dots}{optional arguments that can be passed to the reprojetion algorithm}  
}
\value{The output is list of objects of class \code{"SpatialPixelsDataFrame"} where the number of elements in the list corresponds to the number of standard depths.
}
\note{If the input object is of class \code{"SpatialPixelsDataFrame"}, the method by default uses \code{FWTools} (\code{gdalwarp} command) to resample grids, otherwise the \code{raster::projectRaster} command is passed. \href{http://fwtools.maptools.org/}{FWTools} must be installed separately. \cr 
Note: this operation can be time consuming for large areas (e.g. >> 1e6 pixels).}
\references{
\itemize{
  \item Bivand, R.S., Pebesma, E.J., and \enc{Gómez}{Gomez}-Rubio, V., (2008) \href{http://www.asdar-book.org/}{Applied Spatial Data Analysis with R}. Springer, 378 p.
  \item FWTools (\url{http://fwtools.maptools.org})
  \item Raster package (\url{http://CRAN.R-project.org/package=raster})
}                 
}
\author{ Tomislav Hengl }
\seealso{ \code{\link{spc}}, \code{\link{geosamples-class}}, \code{plotKML::reproject} }
\examples{
## grids Ebergotzen:
library(plotKML)
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
\dontrun{##  convert to spatial components:
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
# create 3D locations in the original coordinate system:
eberg_3Dxy <- sp3D(eberg_spc@predicted)
# wrapper function to create 3D locations in the default WGS84 system:
eberg_3D <- make.3Dgrid(eberg_spc@predicted)
image(eberg_3D[[1]]["PC1"])
# downscale 100 m resolution imagery to 25 m:
data(eberg_grid25)
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")
eberg_grid25@data <- cbind(eberg_grid25@data, 
    gdalwarp(eberg_grid, pixsize=eberg_grid25@grid@cellsize[1], 
    GridTopology=eberg_grid25@grid, resampling_method="cubicspline")@data)
# this function requires FWTools!
}
}
\keyword{methods}
