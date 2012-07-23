\name{MaxEnt}
\docType{methods}
\alias{MaxEnt}
\alias{MaxEnt,ppp,SpatialPixelsDataFrame-method}
\title{Prediction and cross-validation Maximum Entropy}
\description{Runs \href{http://www.cs.princeton.edu/~schapire/maxent/}{MaxEnt} algorithm on a set of observations (\code{"ppp"} class from the \href{http://www.spatstat.org}{spatstat} package) and environmental covariates (of class \code{"SpatialPixelsDataFrame"}) and returns predicted probability of occurrence and cross-validation of models with presence/absence data.}
\usage{
MaxEnt(occurrences, covariates, ...)
}
\arguments{
  \item{occurrences}{object of type \code{"ppp"}; occurrences}
  \item{covariates}{object of type \code{"SpatialPixelsData"}; list of covariate layers}
  \item{\dots}{optional arguments}
}
\details{Optional arguments include:
\describe{
  \item{nfold}{object of type \code{"integer"}; number of folds used for cross-validation (see \code{dismo::randomPoints})}
  \item{Npoints}{object of type \code{"integer"}; number of points used for cross-validation (see \code{dismo::randomPoints})}
  \item{sciname}{object of type \code{"character"}; usually species latin name (it can also be a surveyor's team name or a sampling design)}
  \item{period}{object of type \code{"Date"}; sampling period}
  \item{\dots}{additional arguments that can be passed to \code{dismo::predict}}
}
}
\value{Returns an object of type \code{"SpatialMaxEntOutput"} with the following slots: \code{sciname} (usually latin "genus" and "species" name), \code{occurrences} (occurrence-only records), \code{TimeSpan.begin} (begin of sampling), \code{TimeSpan.end} (end of sampling), \code{maxent} (object of class \code{"MaxEnt"} produced as an output of the \code{dismo::maxent} function), \code{sp.domain} (assumed spatial domain), and \code{predicted} (results of prediction produced using the MaxEnt software).}
\note{MaxEnt is one of the standard tools used in ecology for Niche analysis and species distribution modelling. What makes it especially robust is the fact that it can take both continuous and factor data as inputs, and has no requirements considering the distribution of covariates (\href{http://dx.doi.org/10.1016/j.ecolmodel.2005.03.026}{Phillips et al., 2006}). In the example below, I use MaxEnt to analyze sampling bias (i.e. mis-representation of feature space, or sampling preference) for a given soil sampling design. For more information on how to install MaxEnt and use it in R, see \href{http://CRAN.R-project.org/package=dismo}{dismo} package documentation.}
\author{Tomislav Hengl}
\references{
\itemize{
\item Phillips, S.J., Anderson, R.P., Schapire, R.E., (2006) \href{http://dx.doi.org/10.1016/j.ecolmodel.2005.03.026}{Maximum entropy modeling of species geographic distributions}. Ecological Modelling, 190:231-259.
\item MaxEnt software (\url{http://www.cs.princeton.edu/~schapire/maxent/})
\item Dismo package (\url{http://CRAN.R-project.org/package=dismo})
}
}
\seealso{ \code{dismo::maxent}, \code{plotKML::SpatialMaxEntOutput-class}}
\examples{
\dontrun{# load data:
library(plotKML)
data(eberg)
data(eberg_grid)
# prepare data for spatial analysis:
eberg.xy <- eberg[runif(nrow(eberg)) < .3,]
coordinates(eberg.xy) <- ~X+Y
proj4string(eberg.xy) <- CRS("+init=epsg:31467")
# format gridded data:
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
# convert to a "ppp" object:
library(spatstat)
mg_owin <- as.owin(eberg_grid[1])
eberg.ppp <- ppp(x=coordinates(eberg.xy)[,1], y=coordinates(eberg.xy)[,2], window=mg_owin)
# run MaxEnt analysis (we are interested in possible sampling bias or mis-representation):
me.eberg <- MaxEnt(occurrences=eberg.ppp, covariates=eberg_grid)
## NOTE: MaxEnt can be time-consuming!
# plot the results:
par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,0.5), oma=c(0,0,0,0))
image(as(me.eberg@predicted, "SpatialPixelsDataFrame"), col=rev(heat.colors(25)), xlab="", ylab="")
points(me.eberg@occurrences, pch="+", cex=.7)
image(me.eberg@sp.domain, col="grey", xlab="", ylab="")
}
}
\keyword{methods}