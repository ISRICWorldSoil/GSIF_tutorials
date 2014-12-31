\name{spsample.prob}
\docType{methods}
\alias{spsample.prob}
\alias{spsample.prob,SpatialPoints,SpatialPixelsDataFrame-method}
\title{Estimate inclusion probabilities of a sampling plan (points)}
\description{Estimates inclusion probabilities as an average between the kernel density estimation (geographical spreading of points) and MaxLike analysis (feature space spreading of points). The output \code{'iprob'} indicates whether the sampling plan  has systematically missed some important locations / features, and can be used as an input for geostatistical modelling.}
\usage{ 
\S4method{spsample.prob}{SpatialPoints,SpatialPixelsDataFrame}(observations, 
   covariates, quant.nndist=.95, n.sigma, \dots)
}
\arguments{
  \item{observations}{object of class \code{SpatialPoints}; sampling locations}
  \item{covariates}{object of class \code{SpatialPixelsDataFrame}; list of covariates of interest}
  \item{quant.nndist}{numeric; threshold probability to determine the search radius (sigma)}
  \item{n.sigma}{numeric; size of sigma used for kernel density estimation (optional)}
  \item{\dots}{other optional arguments that can be passed to function \code{spatstat::density}}
}
\value{
Returns a list of objects where \code{'iprob'} (\code{"SpatialPixelsDataFrame"}) is the map showing the estimated inclusion probabilities.
}
\note{
Inclusion probabilities for geographical space are derived using kernel density estimator. The sampling intensities are converted to probabilities by deviding the sampling intensity by the maximum sampling intensity for the study area (Baddeley, 2008). The inclusion probabilities for feature space are determined using MaxLike algorithm (Royle et al., 2012). The lower the average inclusion probability for the whole study area, the lower the representation efficiency of a sampling plan. \cr 
This function can be time consuming for large grids.
}
\references{
\itemize{
\item Baddeley, A. (2008) \href{http://www.csiro.au/resources/pf16h.html}{Analysing spatial point patterns in R}. Technical report, CSIRO Australia. Version 4.
\item Royle, J.A., Chandler, R.B., Yackulic, C. and J. D. Nichols. (2012) \href{http://dx.doi.org/10.1111/j.2041-210X.2011.00182.x}{Likelihood analysis of species occurrence probability from presence-only data for modelling species distributions}. Methods in Ecology and Evolution.
}
}
\author{ Tomislav Hengl }
\seealso{ \code{maxlike-package}, \code{spatstat-package} }
\examples{
library(plotKML)
library(maxlike)
library(spatstat)
library(maptools)

data(eberg)
data(eberg_grid)
## existing sampling plan:
eberg.xy <- eberg[runif(nrow(eberg)) < .2,c("X","Y")]
coordinates(eberg.xy) <- ~X+Y
proj4string(eberg.xy) <- CRS("+init=epsg:31467")
## covariates:
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
## convert to continuous independent covariates:
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)

## derive inclusion probability:
iprob <- spsample.prob(eberg.xy, eberg_spc@predicted[1:8])
hist(iprob[[1]]@data[,1])
## compare with random sampling:
rnd <- spsample(eberg_grid, type="random",
     n=length(iprob[["observations"]]))
iprob2 <- spsample.prob(rnd, eberg_spc@predicted[1:8])
## compare the two:
par(mfrow=c(1,2))
plot(raster(iprob[[1]]), zlim=c(0,1), col=SAGA_pal[[1]])
points(iprob[["observations"]])
plot(raster(iprob2[[1]]), zlim=c(0,1), col=SAGA_pal[[1]])
points(iprob2[["observations"]])
}