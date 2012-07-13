\name{spc}
\docType{methods}
\alias{spc}
\alias{spc,SpatialPixelsDataFrame,formula-method}
\title{Derive Soil Predictive Components}
\description{Derives Soil Predictive Components for a given set of covariates. It wraps the \code{stats::prcomp} method and predicts a list principal components for an object of type \code{"SpatialPixelsDataFrame"}.}
\usage{
spc(obj, formulaString, scale.=TRUE, silent = FALSE, ...)
}
\arguments{
  \item{obj}{object of class \code{"SpatialPixelsDataFrame"}; must contain at least two grids}
  \item{formulaString}{object of class \code{"formula"}}
  \item{scale.}{object of class \code{"logical"}; specifies whether covariates need to be scaled}
  \item{silent}{object of class \code{"logical"}; specifies whether to print the progress}
  \item{\dots}{additional arguments that can be passed to \code{stats::prcomp}} 
}
\value{
\code{spc} returns an object of type \code{"SpatialComponents"}.
}
\details{The output of \code{spc} is a list of grids with generic names \code{PC1},\dots,\code{PCp}, where \code{p} is the total number of input grids.}
\note{This method assumes that the input covariates are cross-correlated and hence their overlap can be reduced. The input variables are scaled by default and the missing values will be replaced with 0 values to reduce loss of data due to missing pixels. This operation can be time consuming for large grids.}
\author{ Tomislav Hengl }
\seealso{ \code{stats::prcomp}, \code{\link{SpatialComponents-class}}
}
\examples{
\dontrun{# load data:
library(plotKML)
pal = rev(rainbow(65)[1:48])
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
eberg_spc <- spc(eberg_grid, as.formula("~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6"))
names(eberg_spc@sp) # 11 components on the end;
rd = range(eberg_spc@sp@data[,1], na.rm=TRUE)
spplot(eberg_spc@sp[1:4], at=seq(rd[1], rd[2], length.out=48), col.regions=pal)
}
}
\keyword{methods}
