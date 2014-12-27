\name{spsample.prob}
\docType{methods}
\alias{spsample.prob}
\alias{spsample.prob,SpatialPoints,SpatialPixelsDataFrame-method}
\title{Estimate inclusion probabilities of a sampling plan (points)}
\description{Estimates inclusion probabilities as an average between the kernel density estimation (geographical) and MaxEnt analysis (feature space). The output \code{'iprob'} indicates whether the sampling plan  has systematically missed some important locations / features, and can be used to adjust geostatistical modelling.}
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
Inclusion probabilities for geographical space are derived using kernel density estimator. The sampling intensities are converted to probabilities by deviding the sampling intensity by the maximum sampling intensity for the study area. The lower the average inclusion probability for the whole study area, the lower the representation efficiency of a sampling plan. \cr 
This function can be time consuming for large grids. For more information on how to install MaxEnt and use it in R, see \href{http://CRAN.R-project.org/package=dismo}{dismo} package documentation.
}
\references{
\itemize{
\item Baddeley, A. (2008) \href{http://www.csiro.au/resources/pf16h.html}{Analysing spatial point patterns in R}. Technical report, CSIRO Australia. Version 4.
\item Phillips, S.J., Anderson, R.P., Schapire, R.E., (2006) \href{http://dx.doi.org/10.1016/j.ecolmodel.2005.03.026}{Maximum entropy modeling of species geographic distributions}. Ecological Modelling, 190:231-259.
\item MaxEnt software (\url{http://www.cs.princeton.edu/~schapire/maxent/})
}
}
\author{ Tomislav Hengl }
\seealso{ \code{MaxEnt} }
\examples{
library(plotKML)
library(dismo)
library(rJava)
library(spatstat)
library(maptools)

data(eberg)
data(eberg_grid)
# prepare data for spatial analysis:
eberg.xy <- eberg[runif(nrow(eberg)) < .3,c("X","Y")]
coordinates(eberg.xy) <- ~X+Y
proj4string(eberg.xy) <- CRS("+init=epsg:31467")
# format gridded data:
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")

jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if(file.exists(jar)){
  iprob <- spsample.prob(eberg.xy, eberg_grid)
  ## compare with random sampling:
  rnd <- spsample(eberg_grid, type="random",
      n=length(iprob[["observations"]]))
  iprob2 <- spsample.prob(rnd, eberg_grid)
  par(mfrow=c(1,2))
  plot(raster(iprob[[1]]), zlim=c(0,.8))
  points(iprob[["observations"]])
  plot(raster(iprob2[[1]]), zlim=c(0,.8))
  points(iprob2[["observations"]])
}
}