\name{getID}
\docType{methods}
\alias{getID}
\alias{getID,SpatialPolygons-method}
\title{Derive 1 degree cell IDs}
\description{Derives ID's of the 1 degree cells in the default land mask for a given polygon defining the spatial domain of interest.}
\usage{
\S4method{getID}{SpatialPolygons}(obj, pixsize = 3/3600, empty.tif = FALSE, 
      compress = FALSE, zipname = set.file.extension(tempfile(tmpdir = getwd()), "zip"))
}
\arguments{
  \item{obj}{object of class \code{"SpatialPolygons"}; must be in geographical coordinates (WGS84)}
  \item{pixsize}{grid cell size in decimal degrees (set at 0.0008333333 or 100 m around equator)} 
  \item{empty.tif}{logical; specify whether a GeoTiff mask file should be created}
  \item{compress}{logical; specify whether to compress GeoTiffs}
  \item{zipname}{(optional); zip archive file name}      
}
\value{
The output is a vector of grid cell ID names e.g. \code{W79_N83}. These can be further used to automate digital soil mapping for large areas.
}
\note{This operation can be time consuming for large areas (e.g. continents).}
\author{ Tomislav Hengl }
\seealso{ \code{\link{landmask}}
}
\examples{
library(sp)
## Bounding box for Malawi:
bbox = expand.grid(lon=c(32.67152, 35.915046), lat=c(-17.12721, -9.363796))
bbox[5,] <- bbox[1,]
crs = CRS("+proj=longlat +datum=WGS84")
x <- SpatialPolygons(list(Polygons(list(Polygon(bbox)), ID="1")), proj4string=crs)
ID.lst <- getID(x)
}
\keyword{methods}
