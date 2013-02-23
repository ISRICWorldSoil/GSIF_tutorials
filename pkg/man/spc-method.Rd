\name{spc}
\docType{methods}
\alias{spc}
\alias{spc,SpatialPixelsDataFrame,formula-method}
\alias{spc,list,list-method}
\title{Derive Spatial Predictive Components}
\description{Derives Spatial Predictive Components for a given set of covariates. It wraps the \code{stats::prcomp} method and predicts a list principal components for an object of type \code{"SpatialPixelsDataFrame"}.}
\usage{
\S4method{spc}{SpatialPixelsDataFrame,formula}(obj, formulaString, scale. = TRUE, 
      silent = FALSE, \dots)
\S4method{spc}{list,list}(obj, formulaString, scale. = TRUE, 
      silent = FALSE, \dots)
}
\arguments{
  \item{obj}{object of class \code{"SpatialPixelsDataFrame"} (must contain at least two grids) or a list of objects of type \code{"SpatialPixelsDataFrame"}}
  \item{formulaString}{object of class \code{"formula"} or a list of formulas}
  \item{scale.}{object of class \code{"logical"}; specifies whether covariates need to be scaled}
  \item{silent}{object of class \code{"logical"}; specifies whether to print the progress}
  \item{\dots}{additional arguments that can be passed to \code{stats::prcomp}} 
}
\value{\code{spc} returns an object of type \code{"SpatialComponents"}. This is a list of grids with generic names \code{PC1},\dots,\code{PCp}, where \code{p} is the total number of input grids.
}
\note{This method assumes that the input covariates are cross-correlated and hence their overlap can be reduced. The input variables are scaled by default and the missing values will be replaced with 0 values to reduce loss of data due to missing pixels. This operation can be time consuming for large grids.}
\author{ Tomislav Hengl }
\seealso{ \code{stats::prcomp}, \code{\link{SpatialComponents-class}}
}
\examples{
# load data:
library(plotKML)
library(sp)

pal = rev(rainbow(65)[1:48])
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
names(eberg_spc@predicted) # 11 components on the end;
\dontrun{# plot maps:
rd = range(eberg_spc@predicted@data[,1], na.rm=TRUE)
sq = seq(rd[1], rd[2], length.out=48)
spplot(eberg_spc@predicted[1:4], at=sq, col.regions=pal)
}
}
\keyword{methods}
