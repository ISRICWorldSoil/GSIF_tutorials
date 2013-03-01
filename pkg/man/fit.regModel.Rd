\name{fit.regModel-methods}
\alias{fit.regModel}
\alias{fit.regModel,formula,data.frame,SpatialPixelsDataFrame,character-method}
\title{Fits a regression model to spatial data}
\description{Fits a regression or a trend model (e.g. a GLM) and, if not available, a variogram for the response residuals using the default settings.}
\usage{
\S4method{fit.regModel}{formula,data.frame,SpatialPixelsDataFrame,character}(formulaString, rmatrix, predictionDomain, 
      method = list("GLM", "rpart", "randomForest", "HB")[[1]], 
      dimensions = NULL, family=gaussian, stepwise=TRUE, rvgm, \dots)
}
\arguments{
  \item{formulaString}{object of class \code{"formula"} --- regression model} 
  \item{rmatrix}{object of class \code{"data.frame"}; regression matrix produced as a result of spatial overlay}
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; spatial domain of interest}
  \item{method}{character; family of methods considered e.g. \code{"GLM"}, \code{"rpart"} (regression trees), \code{"randomForest"}}
  \item{dimensions}{character; \code{"3D"}, \code{"2D"}, \code{"2D+T"}, \code{"3D+T"} models}
  \item{family}{family to be passed to the \code{glm} (see examples below)}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{rvgm}{residual variogram (to avoid fitting the variogram set as NULL)}
  \item{\dots}{other optional arguments that can be passed to \code{gstat::fit.variogram}}
}
\details{Produces an object of class \code{"gstatModel"} that contains: (1) fitted regression model (e.g. a GLM, cubist model, randomForest model or a hierarchical Bayes model), (2) fitted variogram, and (c) object of class \code{"SpatialPoints"} with observation locations. To combine overlay and model fitting operations, consider using \code{\link{fit.gstatModel}}.}
\author{ Tomislav Hengl, Gerard B.M. Heuvelink and Bas Kempen }
\seealso{ \code{\link{fit.gstatModel}}, \code{stats::glm}, \code{gstat::fit.variogram}, \code{randomForest::randomForest} }
\examples{
## Meuse data:
library(sp)
library(randomForest)
library(rpart)

loadMeuse()
## produce a regression matrix:
index <- overlay(meuse.grid, meuse)
sel <- !is.na(index)
ov <- cbind(data.frame(meuse[sel,"om"]), meuse.grid[index[sel],])
## skip variogram fitting:
m <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, 
    family=gaussian(log), method="GLM", rvgm=NULL)
m@regModel
m@vgmModel
## fit a GLM with variogram:
m1 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, 
    family=gaussian(log), method="GLM")
m1@vgmModel
## fit a regression tree with variogram:
m2 <- fit.regModel(log1p(om)~dist+ffreq, rmatrix=ov, meuse.grid,
     method="rpart")
}
