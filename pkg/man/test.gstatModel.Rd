\name{test.gstatModel-methods}
\docType{methods}
\alias{test.gstatModel-method}
\alias{test.gstatModel}
\alias{test.gstatModel,SpatialPointsDataFrame,formula,SpatialPixelsDataFrame-method}
\alias{test.gstatModel,geosamples,formula,SpatialPixelsDataFrame-method}
\title{Methods to test predictability of a regression-kriging model}
\description{Tests predictability of a regression-kriging model on a sample data set. Automates model fitting, cross-validation and prediction and prints out: (1) RMSE at validation points under different sampling intensities, (2) number of predictions per second and (3) number of prediction failures (failure = predictions where cross-validation z-scores exceed value of +/- 1.5 or cross-validation residuals exceed three standard deviations of the observed values).}
\usage{ 
\S4method{test.gstatModel}{SpatialPointsDataFrame,formula,SpatialPixelsDataFrame}(observations, formulaString, covariates, Ns, 
    predictionLocations, save.predictions = TRUE, debug.level = 0, nfold = 5, \dots)
\S4method{test.gstatModel}{geosamples,formula,SpatialPixelsDataFrame}(observations, formulaString, covariates, Ns, 
    predictionLocations, save.predictions = TRUE, debug.level = 0, nfold = 5, \dots)
}
\arguments{
  \item{observations}{object of type \code{"SpatialPointsDataFrame"} or \code{"geosamples-class"}}
  \item{formulaString}{object of type \code{"formula"} or a list of formulas}
  \item{covariates}{object of type \code{"SpatialPixelsDataFrame"}, or list of grids}
  \item{Ns}{vector; list of sampling intensities (maximum should not exceed the total number of samples)}
  \item{predictionLocations}{object of class \code{"SpatialPixelsDataFrame"}; if not specified then passes the object \code{covariates}}
  \item{save.predictions}{logical; indicates whether the prediction results should also be saved}
  \item{debug.level}{integer; gstat's setting to hide the progress output}
  \item{nfold}{integer; number of folds for cross-validation}
  \item{\dots}{other optional arguments that can be passed to \code{fit.gstatModel}}  
}
\note{Vector of sampling intensities, if not provided, will be estimated as: sequence of 10 numbers on square root scale (where N minimum is determined as 20 + number of covariates times 10 and N maximum is the total number of observations). Where no model can be fitted, function returns an empty set. This function can be time consuming for large data sets and is hence recommended only for testing a mapping algorithm using sample data.
}
\author{ Tomislav Hengl, Gerard B.M. Heuvelink }
\seealso{ \code{\link{fit.gstatModel}}, \code{\link{gstatModel-class}} }
\examples{
# 2D model:
library(sp)
## load the Meuse data set:
demo(meuse, echo=FALSE)
## model diagnostics:
t1 <- test.gstatModel(meuse, om~dist+ffreq, meuse.grid, 
   fit.family = gaussian(log), Ns=c(80, 155))
t1[[1]]
}
\keyword{methods}