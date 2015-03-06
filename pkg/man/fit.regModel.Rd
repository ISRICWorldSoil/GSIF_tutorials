\name{fit.regModel-methods}
\alias{fit.regModel}
\alias{fit.regModel,formula,data.frame,SpatialPixelsDataFrame,character-method}
\title{Fits a regression model to spatial data}
\description{Fits a regression or a trend model (e.g. a GLM) and, if not available, a variogram for the response residuals using the default settings.}
\usage{
\S4method{fit.regModel}{formula,data.frame,SpatialPixelsDataFrame,character}(formulaString, rmatrix, predictionDomain, 
      method = list("GLM", "rpart", "randomForest", "quantregForest", "lme")[[1]], 
      dimensions = NULL, fit.family = gaussian(), stepwise = TRUE, rvgm, 
      GLS = FALSE, random, steps, subsample, subsample.reg, \dots)
}
\arguments{
  \item{formulaString}{object of class \code{"formula"} --- regression model} 
  \item{rmatrix}{object of class \code{"data.frame"}; regression matrix produced as a result of spatial overlay}
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; spatial domain of interest}
  \item{method}{character; family of methods considered e.g. \code{"GLM"}, \code{"rpart"} (regression trees), \code{"randomForest"} (random forest)}
  \item{dimensions}{character; \code{"3D"}, \code{"2D"}, \code{"2D+T"}, \code{"3D+T"} models}
  \item{fit.family}{family to be passed to the \code{glm} (see examples below)}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{rvgm}{residual variogram (to avoid fitting the variogram set as NULL)}
  \item{GLS}{fit trend model using Generalized Least Squares implemented in the nlme package}
  \item{random}{specify the random function model; see \code{nlme::lme} for more details}
  \item{steps}{integer; the maximum number of steps to be considered for step-wise regression; see \code{stats::step} for more details}
  \item{subsample}{integer; maximum number of observations to be taken for variogram model fitting (to speed up variogram fitting)}
  \item{subsample.reg}{integer; maximum number of observations to be taken for regression model fitting (especially important for randomForest modelling)}
  \item{\dots}{other optional arguments that can be passed to \code{gstat::fit.variogram}}
}
\details{Produces an object of class \code{"gstatModel"} that contains: (1) fitted regression model (e.g. a GLM, cubist model, or randomForest model), (2) fitted variogram, and (c) object of class \code{"SpatialPoints"} with observation locations. To combine overlay and model fitting operations, consider using \code{\link{fit.gstatModel}}.}
\author{ Tomislav Hengl, Gerard B.M. Heuvelink and Bas Kempen }
\seealso{ \code{\link{fit.gstatModel}}, \code{stats::glm}, \code{gstat::fit.variogram}, \code{randomForest::randomForest} }
\examples{
## Meuse data:
library(sp)
library(rpart)
library(nlme)
library(gstat)

## load the Meuse data set:
demo(meuse, echo=FALSE)

## prepare the regression matrix:
ov <- over(meuse, meuse.grid)
ov <- cbind(data.frame(meuse["om"]), ov)
## skip variogram fitting:
m <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, 
    fit.family=gaussian(log), method="GLM", rvgm=NULL)
m@regModel
m@vgmModel
plot(m)
## fit a GLM with variogram:
m1 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, 
    fit.family=gaussian(log), method="GLM")
m1@vgmModel
plot(m1)
## fit a regression tree with variogram:
m2 <- fit.regModel(log1p(om)~dist+ffreq, rmatrix=ov, meuse.grid,
     method="rpart")
plot(m2)
## fit a lme model with variogram:
m3 <- fit.regModel(log1p(om)~dist, rmatrix=ov, meuse.grid,
     method="lme", random=~1|ffreq)
plot(m3)
## fit a randomForest model with variogram
## (no transformation needed!):
m4 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid,
     method="randomForest")
plot(m4)
}
