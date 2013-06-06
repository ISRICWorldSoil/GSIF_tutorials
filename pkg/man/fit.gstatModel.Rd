\name{fit.gstatModel-methods}
\docType{methods}
\encoding{latin1}
\alias{fit.gstatModel-method}
\alias{fit.gstatModel}
\alias{fit.gstatModel,SpatialPointsDataFrame,formula,SpatialPixelsDataFrame-method}
\alias{fit.gstatModel,geosamples,formula,SpatialPixelsDataFrame-method}
\alias{fit.gstatModel,geosamples,formula,list-method}
\alias{fit.gstatModel,geosamples,list,list-method}
\alias{test.gstatModel}
\alias{test.gstatModel,SpatialPointsDataFrame,formula,SpatialPixelsDataFrame-method}
\title{Methods to fit a regression-kriging model}
\description{Tries to automatically fit a 2D or 3D regression-kriging model for a given set of points (object of type \code{"SpatialPointsDataFrame"} or \code{"geosamples"}) and covariates (object of type \code{"SpatialPixelsDataFrame"}). It first fits a regression model (e.g. Generalized Linear Model, regression tree, random forest model or similar) following the \code{formulaString}, then fits variogram for residuals usign the \code{fit.variogram} method from the \href{http://www.gstat.org}{gstat} package. Creates an output object of class \code{\link{gstatModel-class}}.}
\usage{ 
\S4method{fit.gstatModel}{SpatialPointsDataFrame,formula,SpatialPixelsDataFrame}(observations, formulaString, covariates, 
     method = list("GLM", "rpart", "randomForest", "HB")[[1]], 
     dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]],
     family = gaussian, stepwise = TRUE, vgmFun = "Exp", 
     subsample = 5000, \dots)
\S4method{fit.gstatModel}{geosamples,formula,SpatialPixelsDataFrame}(observations, formulaString, covariates, 
     method = list("GLM", "rpart", "randomForest", "HB")[[1]], 
     dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]],
     family = gaussian, stepwise = TRUE, 
     vgmFun = "Exp", subsample = 5000, \dots)
\S4method{fit.gstatModel}{geosamples,formula,list}(observations, formulaString, covariates, 
     method = list("GLM", "rpart", "randomForest", "HB")[[1]], 
     dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]],
     family = gaussian, stepwise = TRUE, 
     vgmFun = "Exp", subsample = 5000, \dots)
\S4method{fit.gstatModel}{geosamples,list,list}(observations, formulaString, covariates, 
     method = list("GLM", "rpart", "randomForest", "HB")[[1]], 
     dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]],
     family = gaussian, stepwise = TRUE, 
     vgmFun = "Exp", subsample = 5000, \dots)
}
\arguments{
  \item{observations}{object of type \code{"SpatialPointsDataFrame"} or \code{"geosamples-class"}}
  \item{formulaString}{object of type \code{"formula"} or a list of formulas}
  \item{covariates}{object of type \code{"SpatialPixelsDataFrame"}, or list of grids}
  \item{method}{character; family of methods considered e.g. \code{"GLM"}}  
  \item{dimensions}{character; \code{"3D"}, \code{"2D"}, \code{"2D+T"}, \code{"3D+T"} models}  
  \item{family}{character string defyning the GLM family (for more info see \code{stats::glm})}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{vgmFun}{variogram function (\code{"Exp"} by default)}
  \item{subsample}{integer; maximum number of observations to be taken for model fitting (to speed up variogram fitting)}
  \item{\dots}{other optional arguments that can be passed to \code{glm} and/or \code{fit.variogram}}  
}
\details{The GLM method by default assumes that the target variable follows a normal distribution \code{family = gaussian}. Other possible families are:
\describe{
 \item{normal distribution}{\code{family = gaussian} (default setting)}
 \item{log-normal distribution}{\code{family = gaussian(log)}}
 \item{binomial variable}{\code{family = binomial(logit)}}
 \item{variable following a poisson distribution}{\code{family = poisson(log)}}
}
}
\note{Residuals (response residuals from the model) will be checked for normality and problems reported by default. The warning messages should be taken with care, as when the sample size is small, even big departures from normality will not be reported; when the sample size is large, even the smallest deviation from normality might lead to a warning. Likewise, if the variogram fitting fails, consider fitting a variogram manually or using the \code{\link{fit.vgmModel}} method.
}
\references{
\itemize{
\item chapter 8 ``Interpolation and Geostatistics'' in Bivand, R., Pebesma, E., Rubio, V., (2008) \href{http://asdar-book.org/}{Applied Spatial Data Analysis with R}. Use R Series, Springer, Heidelberg, pp. 378.
\item Hengl, T. (2009) \href{http://spatial-analyst.net/book/}{A Practical Guide to Geostatistical Mapping}, 2nd Edt. University of Amsterdam, www.lulu.com, 291 p.
}
}
\author{ Tomislav Hengl, Gerard B.M. Heuvelink and Bas Kempen }
\seealso{ \code{\link{gstatModel-class}}, \code{\link{fit.regModel}}, \code{\link{geosamples-class}}, \code{stats::glm}, \code{gstat::fit.variogram}
}
\examples{
# 2D model:
library(sp)
library(boot)
library(aqp)
library(plyr)
library(rpart)
library(splines)

## load the Meuse data set:
demo(meuse, echo=FALSE)

## fit a regression-tree:
omm <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid, method="rpart")
summary(omm@regModel)
## plot a regression-tree:
plot(omm@regModel, uniform=TRUE); text(omm@regModel, use.n=TRUE, all=TRUE, cex=.8)
omm@vgmModel    
## fit a randomForest model:
omm <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid, method="randomForest")
summary(omm@regModel)
## plot the estimated error for number of bootstrapped trees:
plot(omm@regModel)
omm@vgmModel
## fit a GLM with a gaussian log-link:
omm <- fit.gstatModel(meuse, om~dist+ffreq, meuse.grid, family = gaussian(log))
## it was succesful!
summary(omm@regModel)
v = omm@vgmModel; class(v) = c("variogramModel", "data.frame")
om.sp <- SpatialPointsDataFrame(omm@sp, data = omm@regModel$model)
\dontrun{## plot a variogram:
vve = variogram(om~1, om.sp)
vvres = variogram(om~dist+ffreq, om.sp)
plot(x=vvres$dist, y=vvres$gamma, pch=20, xlab='distance',
    cex=1.1, ylab='gamma', ylim = c(0, max(vve$gamma)))
points(x=vve$dist, y=vve$gamma, pch="+", cex=1.1, col = "grey")
vline <- variogramLine(v, maxdist=max(vvres$dist), n=length(vvres$dist))
lines(x=vline$dist, y=vline$gamma)
}
om.rk <- predict(omm, meuse.grid)
\dontrun{## plot the results in Google Earth:
plotKML(om.rk)
## model diagnostics:
t1 <- test.gstatModel(meuse, om~dist+ffreq, meuse.grid, family = gaussian(log), Ns=c(80, 100, 120, 140, 155))
}

## binary variable (0/1):
som <- fit.gstatModel(meuse, I(soil==1)~dist+ffreq, meuse.grid, family = binomial(logit))
summary(som@regModel)
som.sp <- SpatialPointsDataFrame(som@sp, data = som@regModel$model)
som.sp$I.soil = as.numeric(som.sp@data[,1])
\dontrun{# plot a variogram:
v2 = som@vgmModel; class(v2) = c("variogramModel", "data.frame")
plot(variogram(I.soil~dist+ffreq, som.sp), v2)
}
som.rk <- predict(som, meuse.grid)
\dontrun{# plot the results in Google Earth:
plotKML(som.rk)
}

## 3D model:
library(plotKML)
data(eberg)
## list columns of interest:
s.lst <- c("ID", "soiltype", "TAXGRSC", "X", "Y")
h.lst <- c("UHDICM","LHDICM","SNDMHT","SLTMHT","CLYMHT")
sel <- runif(nrow(eberg))<.05
## get sites table:
sites <- eberg[sel,s.lst]
## get horizons table:
horizons <- getHorizons(eberg[sel,], idcol="ID", sel=h.lst)
## create object of type "SoilProfileCollection"
eberg.spc <- join(horizons, sites, type='inner')
depths(eberg.spc) <- ID ~ UHDICM + LHDICM
site(eberg.spc) <- as.formula(paste("~", paste(s.lst[-1], collapse="+"), sep=""))
coordinates(eberg.spc) <- ~X+Y
proj4string(eberg.spc) <- CRS("+init=epsg:31467")
## convert to geosamples:
eberg.geo <- as.geosamples(eberg.spc)
## covariates:
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
glm.formulaString = as.formula(paste("SNDMHT ~ ", 
  paste(names(eberg_grid), collapse="+"), "+ ns(altitude, df=4)"))
SNDMHT.m <- fit.gstatModel(observations=eberg.geo, glm.formulaString, 
  covariates=eberg_grid)
## problems with the variogram?
\dontrun{## remove classes from the PRMGEO6 that are not represented in the model:
fix.c = levels(eberg_grid$PRMGEO6)[!(levels(eberg_grid$PRMGEO6) \%in\% levels(SNDMHT.m@regModel$model$PRMGEO6))]
summary(eberg_grid$PRMGEO6)
for(j in fix.c){
  eberg_grid$PRMGEO6[eberg_grid$PRMGEO6 == j] <- levels(eberg_grid$PRMGEO6)[7]
}
## prepare new locations:
new3D <- sp3D(eberg_grid)
## regression only:
SNDMHT.rk.sd1 <- predict(SNDMHT.m, new3D[[1]], vgmmodel=NULL)
## regression-kriging:
SNDMHT.rk.sd1 <- predict(SNDMHT.m, new3D[[1]])
## plot the results in Google Earth:
plotKML(SNDMHT.rk.sd1, z.lim=c(5,85))
}
}
\keyword{methods}