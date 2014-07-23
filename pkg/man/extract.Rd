\name{extract}
\docType{methods}
\alias{extract}
\alias{extract.list}
\alias{extract,SpatialPoints,character-method}
\alias{extract,SpatialPointsDataFrame,character-method}
\title{Extracts values at points from a list of files}
\description{Overlays and extracts values at points from a list of files (GDAL compatible). Extends the \code{extract} function from the \href{http://cran.r-project.org/package=raster}{raster} package. This method is especially suitable for extracting values of a large list of rasters that have not been organized into a mosaick (a virtual stack), for example a list of Landsat scenes.}
\usage{
\S4method{extract}{SpatialPoints,character}(x, y, path=".", ID="SOURCEID", method="simple", is.pattern=FALSE, force.projection = TRUE, NAflag = "", \dots)
\S4method{extract}{SpatialPointsDataFrame,character}(x, y, path=".", ID="SOURCEID", method="simple", is.pattern=FALSE, force.projection = TRUE, NAflag = "", \dots)
}
\arguments{
  \item{x}{object of class \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"}} 
  \item{y}{character; list of files that can be read using the \code{raster} function}
  \item{path}{optional working directory where the files are stored}
  \item{ID}{character; column name for the unique identifier that will be generated} 
  \item{method}{character; resampling method (see \code{raster::extract})}
  \item{is.pattern}{logical; specifies whether the list is a pattern}
  \item{force.projection}{logical; specifies whether the reprojection should be ignored}
  \item{NAflag}{character; missing value flag (all missing values are removed by default)}
  \item{\dots}{additional arguments that can be passed to the \code{raster::extract} function}
}
\note{The method will try to reproject the values to the native coordinate system (if ). Missing values are removed by default. If file list is indicated using a pattern, and if multiple rasters covering the same area are found, the method will aggregate the values and derive a mean value.
}
\author{ Tomislav Hengl }
\seealso{ \code{raster::extract}, \code{\link{warp}} }
\keyword{methods}