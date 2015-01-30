\name{tile}
\docType{methods}
\alias{tile}
\alias{tile,SpatialPixelsDataFrame-method}
\alias{tile,SpatialPointsDataFrame-method}
\alias{tile,SpatialLinesDataFrame-method}
\alias{tile,SpatialPolygonsDataFrame-method}
\alias{tile,RasterLayer-method}
\title{Tiles (subsets or clips) a spatial object to regular blocks}
\description{Tiles objects of class \code{"Spatial*"} or \code{"RasterLayer*"} into regular blocks.}
\usage{
\S4method{tile}{SpatialPointsDataFrame}(x, y, block.x, \dots)
\S4method{tile}{SpatialPixelsDataFrame}(x, y, block.x, \dots)
\S4method{tile}{SpatialPolygonsDataFrame}(x, y, block.x, tmp.file = TRUE, 
    program, show.output.on.console = FALSE, \dots)
\S4method{tile}{SpatialLinesDataFrame}(x, y, block.x, tmp.file = TRUE, 
    program, show.output.on.console = FALSE, \dots)
\S4method{tile}{RasterLayer}(x, y, block.x, tmp.file = TRUE, 
    program, show.output.on.console = FALSE, \dots)
}
\arguments{
  \item{x}{object of class \code{"Spatial*"} \code{"RasterLayer"}}
  \item{y}{list of \code{"SpatialPolygons"}; if missing will be derived based on \code{block.x}}
  \item{block.x}{numeric; size of block in meters or corresponding mapping units}
  \item{tmp.file}{logical; specifies whether to generate a temporary file}
  \item{program}{character; location of the auxiliary program in the system}
  \item{show.output.on.console}{logical; specifies whether to print the progress of a function}
  \item{\dots}{optional arguments that can be passed to the \code{\link{getSpatialTiles}}}
}
\value{Returns a list of objects of the same class as the input object.}
\details{When working with objects of type \code{"SpatialLinesDataFrame"}, \code{"SpatialPolygonsDataFrame"} and or \code{"RasterLayer"}, the function looks for FWTools binary files \code{ogr2ogr} and \code{warp}. \href{http://fwtools.maptools.org/}{FWTools} is a separate program and must be installed separately.}
\author{ Tomislav Hengl}
\seealso{ \code{\link{getSpatialTiles}} }
\examples{
## spatial pixels:
library(sp)
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
tl <- getSpatialTiles(meuse.grid, block.x=1000)
image(meuse.grid)
lines(as(tl, "SpatialLines"))
## all at once:
pix.lst <- tile(meuse.grid, block.x=1000)
\dontrun{## lines:
library(plotKML)
data(eberg_contours)
line.lst <- tile(eberg_contours, block.x=5000)
spplot(line.lst[[1]][2])
## polygons:
data(eberg_zones)
## this one requires ogr2ogr function:
pol.lst <- tile(eberg_zones, block.x=5000)
spplot(pol.lst[[1]][1])
## raster files via rgdal:
library(rgdal)
fn = system.file("pictures/SP27GTIF.TIF", 
 package = "rgdal")
obj <- GDALinfo(fn)
ras.lst <- getSpatialTiles(obj, block.x=1000)
offset <- c(ras.lst$offset.y[1], ras.lst$offset.x[1])
region.dim <- c(ras.lst$region.dim.y[1], 
  ras.lst$region.dim.x[1])
## read the first tile:
SP27GTIF_T1 <- readGDAL(fn, offset=offset, 
  region.dim=region.dim)
str(SP27GTIF_T1)
}
}
\keyword{methods}
