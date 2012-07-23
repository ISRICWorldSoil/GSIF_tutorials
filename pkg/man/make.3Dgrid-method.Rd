\name{make.3Dgrid}
\docType{methods}
\encoding{latin1}
\alias{make.3Dgrid}
\alias{sp3D}
\alias{make.3Dgrid,SpatialPixelsDataFrame-method}
\alias{make.3Dgrid,RasterBrick-method}
\title{Methods to prepare 3D prediction locations}
\description{Generates a list of objects of type \code{"SpatialPixelsDataFrame"} with \code{longitude}, \code{latitude} and \code{altitude} coordinates (these names are used by default for compatibility with the \code{\link{geosamples-class}}).}
\usage{
make.3Dgrid(obj, ...)
}
\arguments{
  \item{obj}{object of class \code{"SpatialPixelsDataFrame"} or \code{"RasterBrick"}} 
  \item{\dots}{optional arguments that can be passed to the reprojetion algorithm}  
}
\details{Optional arguments include:
\describe{
  \item{proj4s}{character; proj4string describing the target coordinate system}
  \item{pixsize}{grid cell size in decimal degrees (set by default at 1/1200 (0.0008333333 or 100 m around equator)} 
  \item{resample}{character; resampling method to be passed the reprojection algorithm}
  \item{NAflag}{character; missing value flag}
  \item{stdepths}{numeric; list of standard depths}                                             
  \item{tmp.file}{logical; specifies whether a temporary file name should be generated}
  \item{show.output.on.console}{logical; specifies whether to print out the progress}
}
The output is list of objects of class \code{"SpatialPixelsDataFrame"} where the number of elements in the list corresponds to the number of standard depths.}
\note{If the input object is of class \code{"SpatialPixelsDataFrame"}, the method by default uses \code{FWTools} (\code{gdalwarp} command) to resample grids, otherwise the \code{raster::projectRaster} command is passed. \href{http://fwtools.maptools.org/}{FWTools} is a separate application and must be installed separately. This operation can be time consuming for large areas (e.g. continents).}
\references{
\itemize{
  \item Bivand, R.S., Pebesma, E.J., and \enc{Gómez}{Gomez}-Rubio, V., (2008) \href{http://www.asdar-book.org/}{Applied Spatial Data Analysis with R}. Springer, 378 p.
}
}
\author{ Tomislav Hengl }
\seealso{ \code{\link{spc}}, \code{\link{geosamples-class}}, \code{plotKML::reproject} }
\examples{
\dontrun{# prepare prediction locations for Ebergotzen:
library(plotKML)
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
# convert to spatial components:
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
# create 3D locations in the original coordinate system:
eberg_3Dxy <- sp3D(eberg_spc@predicted)
# wrapper function to create 3D locations in the default system:
eberg_3D <- make.3Dgrid(eberg_spc@predicted)
image(eberg_3D[[1]]["PC1"])
}
}
\keyword{methods}
