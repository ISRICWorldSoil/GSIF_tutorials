\name{predict.gstatModel-method}
\docType{methods}
\alias{predict,gstatModel-method}
\alias{predict.gstatModel}
\alias{predict,list-method}
\alias{predict.gstatModelList}
\title{Predict from an object of class \code{"gstatModel"}}
\description{Predicts from an object of class \code{\link{gstatModel-class}} using new prediction locations. It combines predictions by regression (e.g. GLM) and interpolation of residuals (kriging).}
\usage{ 
\S4method{predict}{gstatModel}(object, 
     predictionLocations, nmin = 10, nmax = 30, debug.level = -1, 
     predict.method = c("RK", "KED")[1], nfold = 5, verbose = FALSE, 
     nsim = 0, mask.extra = TRUE, block = predictionLocations@grid@cellsize, 
     zmin = -Inf, zmax = Inf, subsample = length(object@sp), 
     coarsening.factor = 1, vgmmodel = object@vgmModel,
     subset.observations = !is.na(object@sp@coords[,1]), betas = c(0,1), \dots)
\S4method{predict}{list}(object, 
     predictionLocations, nmin = 10, nmax = 30, debug.level = -1, 
     predict.method = c("RK", "KED")[1], nfold = 5, verbose = FALSE, 
     nsim = 0, mask.extra = TRUE, block = predictionLocations@grid@cellsize, 
     zmin = -Inf, zmax = Inf, subsample = length(object@sp), \dots)
}
\arguments{
  \item{object}{object of type \code{"gstatModel"}}
  \item{predictionLocations}{object of type \code{"SpatialPixelsDataFrame"} prediction locations (must contain all covariates from the model)}
  \item{nmin}{integer; minimum number of nearest observations sent to \code{gstat::krige}}
  \item{nmax}{integer; maximum number of nearest observations sent to \code{gstat::krige}}
  \item{debug.level}{integer; default debug level mode sent to \code{gstat::krige}}
  \item{predict.method}{character; mathematical implementation of the \code{gstat::krige} interpolation method with covariates: Regression-Kriging (RK) or Kriging with External Drift (KED)}
  \item{nfold}{integer; n-fold cross validation sent to \code{gstat::krige.cv}}  
  \item{verbose}{logical; specifies whether to supress the progress bar of the \code{gstat::krige.cv}}
  \item{nsim}{integer; triggers the geostatistical simulations}
  \item{mask.extra}{logical; specifies whether to mask out the extrapolation pixels (prediction variance exceeding the global variance)}
  \item{block}{numeric; support size (block predictions by default)}
  \item{zmin}{numeric; lower physical limit for the target variable}
  \item{zmax}{numeric; upper physical limit for the target variable}
  \item{subsample}{integer; sub-sample point observations to speed up the processing}
  \item{coarsening.factor}{integer; coarsening factor (1:5) to speed up the processing}
  \item{vgmmodel}{object of class \code{data.frame} corresponding to the \code{gstat::vgm} variogram}
  \item{subset.observations}{logical; vector specifying the subset of observations used for interpolation}
  \item{betas}{numeric; vector of the beta coefficients to be passed to the \code{gstat::krige}}
  \item{\dots}{other optional arguments that can be passed to \code{gstat::krige} and/or \code{predict.glm}}  
}
\details{Selecting \code{predict.method = "KED"} invokes simple kriging with external drift with betas set at 0 (intercept) and 1 (regression predictions i.e. single covariate). This assumes that the regression model already results in an unbiased estimator of the trend model.\cr 
If not speficied otherwise, \code{subset.observations} by default selects only obserations within the spatial domain (bounding box) of the \code{predictionLocations} plus 50\% of the one third of the extent of the area.\cr 
Predictions can be speed up by using a larger \code{coarsening.factor} e.g. 2 to 5, in which case the ordinary kriging on residuals will run at a coarser resolution, and the output would be then downscaled to the original resolution using splines (via the \code{\link{gdalwarp}} method). In the case of \code{predict.method = RK}, the kriging variance is derived as a sum of the GLM variance and the OK variance, which is statistically sub-optimal.
}
\note{Predictions using \code{predict.method = "KED"} (the default gstat setting) can be time consuming for large data set and can result in instabilities (singular matrix problems) if the search radius is small and/or if all covariates contain exactly the same values.}
\author{ Tomislav Hengl, Gerard B.M. Heuvelink and Bas Kempen}
\seealso{ \code{\link{gstatModel-class}}, \code{\link{fit.gstatModel}}
}
\keyword{methods}