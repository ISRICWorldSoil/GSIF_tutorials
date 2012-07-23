\name{spfkm}
\docType{methods}
\alias{spfkm}
\alias{spmultinom}
\alias{spfkm,formula,SpatialPointsDataFrame,SpatialPixelsDataFrame-method}
\alias{spmultinom,formula,matrix,SpatialPixelsDataFrame-method}
\title{Supervised fuzzy \emph{k}-means on spatial pixels}
\description{Runs supervised fuzzy \emph{k}-means using a list of covariates layers provided as \code{"SpatialPixelsDataFrame"}. If class centres and variances are not provided by the users, it first fits a multinomial logistic regression model (\code{spmultinom}) using the point data and a regression matrix, then predicts the class centres and variances based on the output from the \code{nnet::multinom}.}
\usage{
spfkm(formulaString, observations, covariates, class.c  = NULL, class.sd = NULL, fuzzy.e = 1.2)
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
\author{Tomislav Hengl}
\references{
\itemize{
\item Burrough, P. A., Gaans, P.F.M., and Van Hootsmans, R., (1997) \href{http://dx.doi.org/10.1016/S0016-7061(97)00018-9}{Continuous classification in soil survey: spatial correlation, confusion and boundaries}. Geoderma, 77(2-4), 115--135.
\item Multinomial logistic regression (\url{http://en.wikipedia.org/wiki/Multinomial_logit})
\item Nnet package (\url{http://CRAN.R-project.org/package=nnet})
}
}
\seealso{ \code{nnet::multinom}, \code{\link{SpatialMemberships-class}}}
\examples{
\dontrun{# load data:
library(plotKML)
data(eberg)
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
# plot memberships:
pal = seq(0, 1, 1/50)
spplot(eberg_sm@mu, col.regions=pal)
image(eberg_sm@mu[1], col=pal)
text(eberg@coords, paste(eberg$soiltype), cex=.6, col="black")
# classes predicted:
eberg_sm@predicted$soiltype <- as.factor(eberg_sm@predicted$soiltype)
Ls = length(levels(eberg_sm@predicted$soiltype))
pnts = list("sp.points", eberg, pch="+", cex=.6, col="black")
spplot(eberg_sm@predicted, col.regions=rainbow(Ls)[rank(runif(Ls))], sp.layout=pnts)
spplot(eberg_sm@mu[2], at=pal, col.regions=rev(heat.colors(50)), sp.layout=pnts)
}
}
\keyword{methods}