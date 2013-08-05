\name{gdalwarp}
\docType{methods}
\alias{gdalwarp}
\alias{gdalwarp,SpatialPixelsDataFrame-method}
\alias{gdalwarp,RasterLayer-method}
\title{GDAL warp function from FWTools}
\description{Runs GDAL warp program available via FWTools.}
\usage{
\S4method{gdalwarp}{SpatialPixelsDataFrame}(obj, proj4s = proj4string(obj), 
      GridTopology = NULL, pixsize, 
      resampling_method = "bilinear", 
      NAflag = get("NAflag", envir = GSIF.opts), 
      tmp.file = FALSE, show.output.on.console = FALSE, program)
\S4method{gdalwarp}{RasterLayer}(obj, proj4s = proj4string(obj), 
      GridTopology = NULL, pixsize, 
      resampling_method = "bilinear", 
      NAflag = get("NAflag", envir = GSIF.opts), 
      tmp.file = FALSE, show.output.on.console = FALSE, program)
}
\arguments{
  \item{obj}{object of class \code{"SpatialPixelsDataFrame"} or class \code{"RasterLayer"}} 
  \item{proj4s}{character; proj4string describing the target coordinate system}
  \item{GridTopology}{optional grid topology from sp package}
  \item{pixsize}{grid cell size in decimal degrees} 
  \item{resampling_method}{character; resampling method to be passed the reprojection algorithm}
  \item{NAflag}{character; missing value flag}
  \item{tmp.file}{logical; specifies whether a temporary file name should be generated}
  \item{show.output.on.console}{logical; specifies whether to print out the progress}
  \item{program}{location of the GDAL warp program on the local machine}  
}
\note{\href{http://fwtools.maptools.org/}{FWTools} must be installed separately.
}
\author{ Tomislav Hengl }
\seealso{ \code{\link{make.3Dgrid}}, \code{plotKML::reproject} }
\keyword{methods}