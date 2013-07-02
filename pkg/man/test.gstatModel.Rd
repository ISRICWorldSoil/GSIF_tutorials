\name{test.gstatModel-methods}
\docType{methods}
\alias{test.gstatModel-method}
\alias{test.gstatModel}
\alias{test.gstatModel,SpatialPointsDataFrame,formula,SpatialPixelsDataFrame-method}
\title{Methods to test predictability of a regression-kriging model}
\description{Test predictability of a regression-kriging model. Automates model fitting, cross-validation and prediction and prints out: (1) amount of variance explained under different sampling intensities, (2) number of predictions per second and (3) number of prediction failures (failure = predictions where cross-validation z-scores exceed +/- 1.5 or cross-validation residuals exceed three standard deviations).}
\usage{ 
\S4method{test.gstatModel}{SpatialPointsDataFrame,formula,SpatialPixelsDataFrame}(observations, formulaString, covariates, Ns, 
 predictionLocations = covariates, save.predictions = FALSE, svar, debug.level = 0, \dots)
}
\arguments{
  \item{observations}{object of type \code{"SpatialPointsDataFrame"} or \code{"geosamples-class"}}
  \item{formulaString}{object of type \code{"formula"} or a list of formulas}
  \item{covariates}{object of type \code{"SpatialPixelsDataFrame"}, or list of grids}
  \item{Ns}{vector; list of sampling intensities (maximum should not exceed the total number of samples)}
  \item{predictionLocations}{object of class \code{"SpatialPixelsDataFrame"}; if not specified then the same as the covariates}
  \item{save.predictions}{logical; indicates whether the prediction results should also be saved}
  \item{svar}{numeric; total/global variance (optional)}
  \item{debug.level}{integer; gstat's setting to hide the progress output}
  \item{\dots}{other optional arguments that can be passed to \code{fit.gstatModel}}  
}
\note{Vector of sampling intensities, if not provided, will be estimated as 10 numbers on square root scale (N minimum is determined as 20 + number of covariates times 10). Where no model can be fitted, function returns an empty set. This function can be time consuming for large data sets.
}
\author{ Tomislav Hengl, Gerard B.M. Heuvelink and Bas Kempen }
\seealso{ \code{\link{fit.gstatModel}}, \code{\link{gstatModel-class}} }
\examples{
# 2D model:
library(sp)
## load the Meuse data set:
demo(meuse, echo=FALSE)
## model diagnostics:
t1 <- test.gstatModel(meuse, om~dist+ffreq, meuse.grid, family = gaussian(log), Ns=c(80, 155))
t1
}
\keyword{methods}