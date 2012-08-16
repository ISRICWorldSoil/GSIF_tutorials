\name{spmultinom}
\docType{methods}
\alias{spmultinom}
\alias{spmultinom,formula,SpatialPointsDataFrame,SpatialPixelsDataFrame-method}
\title{Multinomial logistic regression on spatial objects}
\description{Runs the multinomial logistic regression via \code{nnet::multinom} to produce spatial predictions of the target factor-type variable. It requires point locations of observed classes and a list of covariate layers provided as \code{"SpatialPixelsDataFrame-class"} object. The resulting predicted classes are then used to estimate class centres and variances per class.}
\usage{
\S4method{spmultinom}{formula,SpatialPointsDataFrame,SpatialPixelsDataFrame}(formulaString,
     observations, covariates, class.stats = TRUE, predict.probs = TRUE, \dots)
}
\arguments{
  \item{formulaString}{formula string}
  \item{observations}{object of type \code{"SpatialPointsData"}; occurrences of factors}
  \item{covariates}{object of type \code{"SpatialPixelsData"}; list of covariate layers}
  \item{class.stats}{logical; species wether to estimate class centres}
  \item{predict.probs}{logical; species wether to predict probabilities per class}
  \item{\dots}{optional arguments}
}
\value{Returns an object of type \code{"SpatialMemberships"} with following slots: \code{predicted} (classes predicted by the multinomial logistic regression, \code{model} (the multinomial logistic regression model), \code{mu} (probabilities derived using the mutinom model), \code{class.c} (derived class centres), \code{class.sd} (derived class deviations), \code{confusion} (confusion matrix).}
\author{Bas Kempen and Tomislav Hengl}
\references{
\itemize{
\item Multinomial logistic regression (\url{http://en.wikipedia.org/wiki/Multinomial_logit})
\item Nnet package (\url{http://CRAN.R-project.org/package=nnet})
}
}
\seealso{ \code{\link{spfkm}}, \code{\link{SpatialMemberships-class}}}
\examples{
\dontrun{# load data:
library(plotKML)
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
eberg_sm <- spmultinom(formulaString, eberg, eberg_spc@predicted)
# plot memberships:
pal = seq(0, 1, 1/50)
spplot(eberg_sm@mu, col.regions=pal)
image(eberg_sm@mu[1], col=pal)
text(eberg@coords, paste(eberg$soiltype), cex=.6, col="black")
# classes predicted:
Ls = length(levels(eberg_sm@predicted$soiltype))
pnts = list("sp.points", eberg, pch="+", cex=.6, col="black")
spplot(eberg_sm@predicted, col.regions=rainbow(Ls)[rank(runif(Ls))], sp.layout=pnts)
}
}
\keyword{methods}