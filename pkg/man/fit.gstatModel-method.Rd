\name{fit.gstatModel-methods}
\docType{methods}
\encoding{latin1}
\alias{fit.gstatModel-method}
\alias{fit.gstatModel}
\alias{fit.gstatModel,SpatialPointsDataFrame,formula,SpatialPixelsDataFrame-method}
\alias{fit.gstatModel,geosamples,formula,SpatialPixelsDataFrame-method}
\alias{fit.gstatModel,geosamples,formula,list-method}
\alias{fit.gstatModel,geosamples,list,list-method}
\title{Methods to fit a regression-kriging model}
\description{Tries to automatically fit a 2D or 3D regression-kriging model for a given set of points (object of type \code{"SpatialPointsDataFrame"} or \code{"geosamples"}) and covariates (object of type \code{"SpatialPixelsDataFrame"}). It first fits a Generalized Linear Model (GLM) following the \code{formulaString}, then fits variogram for residuals usign the \code{fit.variogram} method from the \href{http://www.gstat.org}{gstat} package. Creates an output object of class \code{\link{gstatModel-class}}.}
\usage{ 
\S4method{fit.gstatModel}{SpatialPointsDataFrame,formula,SpatialPixelsDataFrame}(observations,
     formulaString, covariates, family = gaussian, stepwise = TRUE, 
     vgmFun = "Exp", rvgm = NULL, \dots)
\S4method{fit.gstatModel}{geosamples,formula,SpatialPixelsDataFrame}(observations,
     formulaString, covariates, methodid, family = gaussian, stepwise = TRUE, 
     vgmFun = "Exp", rvgm = NULL, \dots)
\S4method{fit.gstatModel}{geosamples,formula,list}(observations, formulaString, covariates,
     methodid, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, \dots)
\S4method{fit.gstatModel}{geosamples,list,list}(observations, formulaString, covariates,
     methodid, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, \dots)
}
\arguments{
  \item{observations}{object of type \code{"SpatialPointsDataFrame"} or \code{"geosamples-class"}}
  \item{formulaString}{object of type \code{"formula"} or a list of formulas}
  \item{covariates}{object of type \code{"SpatialPixelsDataFrame"}, or list of grids}
  \item{methodid}{character; method name in the case of using the \code{"geosamples-class"}}
  \item{family}{character string defyning the GLM family (for more info see \code{stats::glm})}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{vgmFun}{variogram function (\code{"Exp"} by default)}
  \item{rvgm}{residual variogram model can be passed from gstat, in which case \code{fit.gstatModel} will not try to fit the variogram}
  \item{\dots}{other optional arguments that can be passed to \code{glm} and/or \code{fit.variogram}}  
}
\details{The method by default assumes that the target variable follows a normal distribution \code{family = gaussian}. Other possible families are:
\describe{
 \item{normal distribution}{\code{family = gaussian} (default setting)}
 \item{log-normal distribution}{\code{family = gaussian(log)}}
 \item{binomial variable}{\code{family = binomial(logit)}}
 \item{variable following a poisson distribution}{\code{family = poisson(log)}}
}
}
\references{
\itemize{
\item Minasny, B., McBratney, A.B. (2007) \href{http://dx.doi.org/10.1016/j.geoderma.2007.04.028}{Spatial prediction of soil properties using EBLUP with the \enc{Matérn}{Matern} covariance function}. Geoderma, 140(4): 324-336.
\item Hengl, T. (2009) \href{http://spatial-analyst.net/book/}{A Practical Guide to Geostatistical Mapping}, 2nd Edt. University of Amsterdam, www.lulu.com, 291 p.
}
}
\author{ Tomislav Hengl and Gerard B.M. Heuvelink}
\seealso{ \code{\link{gstatModel-class}}, \code{\link{glm.sp}}, \code{\link{geosamples-class}}, \code{stats::glm}, \code{gstat::fit.variogram}
}
\examples{
data(meuse)
coordinates(meuse) <- ~x+y
# load grids:
data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
# fit a model with a gaussian log-link:
omm <- fit.gstatModel(meuse, om~dist, meuse.grid, family = gaussian(log))
omm@regModel
om.sp <- SpatialPointsDataFrame(omm@sp, data = omm@regModel$model)
\dontrun{# plot a variogram:
plot(variogram(om~1, om.sp))
}
}
\keyword{methods}