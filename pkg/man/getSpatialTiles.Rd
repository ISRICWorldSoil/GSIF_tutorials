\name{getSpatialTiles}
\docType{methods}
\alias{getSpatialTiles}
\alias{getSpatialTiles,Spatial-method}
\alias{getSpatialTiles,ANY-method}
\alias{makeTiles}
\title{Get a list of tiles (regular blocks)}
\description{Creates a list of tiles (\code{"SpatialPolygons"}) for a given spatial domain i.e. extent. Input can be any object of class \code{"Spatial"} or \code{"GDALobj"}.}
\usage{
\S4method{getSpatialTiles}{Spatial}(obj, block.x, block.y = block.x, 
    overlap.percent = 0, limit.bbox = TRUE, return.SpatialPolygons = TRUE)
\S4method{getSpatialTiles}{ANY}(obj, block.x, block.y = block.x, 
    overlap.percent = 0, limit.bbox = TRUE, return.SpatialPolygons = FALSE)
}
\arguments{
  \item{obj}{object of class \code{"Spatial*"}}
  \item{block.x}{numeric; size of block in x-direction (meters or corresponding mapping units)}
  \item{block.y}{numeric; size of block in y-direction (meters or corresponding mapping units)}
  \item{overlap.percent}{numeric; percentage overlap (must be a positive number)}
  \item{limit.bbox}{logical; specifies whether to limit the extent of tiles to the bounding box only}
  \item{return.SpatialPolygons}{logical; specifies whether to return a list of tiles as \code{"SpatialPolygons"} or a data frame with bounding box coordinates}
}
\value{Returns a list of tiles either as a list of \code{"SpatialPolygons"} or a data frame with with bounding box coordinates.}
\details{The first output tile starts by default at the lower left corner. \code{getSpatialTiles-method} can only be used to generate regular tiles.}
\author{ Tomislav Hengl }
\seealso{ \code{\link{tile}}, \code{sp::spsample} }
\keyword{methods}
