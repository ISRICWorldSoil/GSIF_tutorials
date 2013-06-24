\name{spfkm}
\docType{methods}
\alias{spfkm}
\alias{spfkm,formula,SpatialPointsDataFrame,SpatialPixelsDataFrame-method}
\title{Supervised fuzzy \emph{k}-means on spatial pixels}
\description{Runs supervised fuzzy \emph{k}-means (\href{http://dx.doi.org/10.1080/13658810310001620924}{Hengl et al., 2004}) using a list of covariates layers provided as \code{"SpatialPixelsDataFrame-class"} object. If class centres and variances are not provided, it first fits a multinomial logistic regression model (\code{\link{spmultinom}}), then predicts the class centres and variances based on the output from the \code{nnet::multinom}.}
\usage{
\S4method{spfkm}{formula,SpatialPointsDataFrame,SpatialPixelsDataFrame}(formulaString,
           observations, covariates, class.c = NULL, class.sd = NULL, fuzzy.e = 1.2)
}
\arguments{
  \item{formulaString}{formula string}
  \item{observations}{object of type \code{"SpatialPointsData"}; occurrences of factors}
  \item{covariates}{object of type \code{"SpatialPixelsData"} or \code{"RasterBrick"}; list of covariate layers}
  \item{class.c}{object of type \code{"matrix"}; class centres (see examples below)}
  \item{class.sd}{object of type \code{"matrix"}; class deviations (see examples below)}
  \item{fuzzy.e}{object of type \code{"numeric"}; fuzzy exponent}
}
\value{Returns an object of type \code{"SpatialMemberships"} with following slots: \code{predicted} (classes predicted either by the multinomial logistic regression or fuzzy \emph{k}-means), \code{model} (the multinomial logistic regression model; if available), \code{mu} (memberships derived using the fuzzy \emph{k}-means), \code{class.c} (submitted or derived class centres), \code{class.sd} (submitted or derived class deviations), \code{confusion} (confusion matrix).}
\note{Although \code{nnet::multinom} is consider to be robust and suited for large data sets, function might not converge in some cases or result in artifacts. If this happens try setting up the class centres and variances manually.}
\author{Tomislav Hengl and Bas Kempen}
\references{
\itemize{
\item Burrough, P. A., Gaans, P.F.M., and Van Hootsmans, R., (1997) \href{http://dx.doi.org/10.1016/S0016-7061(97)00018-9}{Continuous classification in soil survey: spatial correlation, confusion and boundaries}. Geoderma, 77(2-4), 115--135.
\item Hengl T., Walvoort D.J.J., Brown, A., (2004) \href{http://dx.doi.org/10.1080/13658810310001620924}{A double continuous approach to visualisation and analysis of categorical maps}. Int. Jou. of Geographical Information Science, 18(2): 183-202.
}
}
\seealso{ \code{\link{spmultinom}}, \code{\link{SpatialMemberships-class}}, \code{nnet::multinom}}
\examples{
# load data:
library(plotKML)
library(sp)

data(eberg)
# subset to 20%:
eberg <- eberg[runif(nrow(eberg))<.2,]
data(eberg_grid)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
# derive soil predictive components:
eberg_spc <- spc(eberg_grid, ~PRMGEO6+DEMSRT6+TWISRT6+TIRAST6)
# predict memberships:
formulaString = soiltype ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
eberg_sm <- spfkm(formulaString, eberg, eberg_spc@predicted)
\dontrun{# plot memberships:
pal = seq(0, 1, 1/50)
spplot(eberg_sm@mu, col.regions=grey(rev(pal)))
# predict soil properties using memberships:
glm.formulaString = as.formula(paste("SNDMHT_A ~ ", 
    paste(names(eberg_sm@mu), collapse="+"), "-1"))
SNDMHT.m2 <- fit.gstatModel(observations=eberg, glm.formulaString, 
    covariates=eberg_sm@mu)
summary(SNDMHT.m2@regModel)
# Coefficients correspond to the class centres;
}
}
\keyword{methods}