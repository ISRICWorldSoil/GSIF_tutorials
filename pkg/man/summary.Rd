\name{summary-methods}
\alias{summary}
\alias{summary,SpatialPredictions-method}
\alias{show,SpatialPredictions-method}
\title{Summarize an object of class \code{"SpatialPredictions"}}
\description{Derives a statistical summary for an object of class \code{"SpatialPredictions"}.}
\usage{
\S4method{summary}{SpatialPredictions}(object)
}
\arguments{
  \item{object}{object of class \code{"SpatialPredictions"}} 
}
\value{
The summary returns a data.frame with the following columns:
\describe{
  \item{\code{"variable"}}{variable name}
  \item{\code{"minium"}}{lowest value observed}
  \item{\code{"maximum"}}{largest value observed}
  \item{\code{"npoints"}}{number of observations}
  \item{\code{"area"}}{lowest value observed}
  \item{\code{"area.units"}}{area units either square-m or square-arcdegrees}
  \item{\code{"covariates"}}{list of covariates used}
  \item{\code{"family"}}{GLM family (if applicable)}
  \item{\code{"RMSE"}}{RMSE derived using cross-validation}
  \item{\code{"tvar"}}{variance percent explained by the model using the cross-validation}
  \item{\code{"npixels"}}{total number of produced pixels}
  \item{\code{"breaks"}}{breaks based on the half RMSE}
  \item{\code{"bonds"}}{lower and upper boundaries for effective classes}
  \item{\code{"Bytes"}}{effective bytes produced (see \href{http://dx.doi.org/10.1016/j.jag.2012.02.005}{Hengl et al (2012)} for more details)}
  \item{\code{"compress"}}{compression algorithm used}
}
}
\details{The function creates a summary table with standard column names. These tell us what is the summary accuracy of the spatial predictions and what are the effective bytes of information produced.}
\author{ Tomislav Hengl }
\references{
\itemize{
\item Hengl, T., Nikolic, M., MacMillan, R.A., (2013) \href{http://dx.doi.org/10.1016/j.jag.2012.02.005}{Mapping efficiency and information content}. International Journal of Applied Earth Observation and Geoinformation, special issue Spatial Statistics Conference, 22: 127--138.
}
}
\seealso{ \code{plotKML::SpatialPredictions-class} }
\examples{
## load observations:
library(sp)
library(rgdal)
library(gstat)
demo(meuse, echo=FALSE)
## fit a model:
omm <- fit.gstatModel(meuse, om~dist, 
  fit.family=gaussian(link="log"), meuse.grid)
show(omm@regModel)
## produce SpatialPredictions:
om.rk <- predict(omm, predictionLocations = meuse.grid)
x = summary(om.rk)
str(x)
}
