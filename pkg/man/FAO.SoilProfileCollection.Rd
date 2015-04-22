\name{FAO.SoilProfileCollection-class}
\docType{class}
\encoding{latin1}
\alias{FAO.SoilProfileCollection-class}
\alias{FAO.SoilProfileCollection}
\alias{munsell}
\title{A class for FAO SoilProfileCollection}
\description{A class for harmonized (FAO) soil profile records. Extends the \code{"SoilProfileCollection"} class from the \href{http://cran.r-project.org/package=aqp?}{aqp} package.}
\section{Slots}{
  \describe{    
    \item{\code{idcol}:}{object of class \code{"character"}; column name containing IDs}
    \item{\code{depthcols}:}{object of class \code{"character"}; two element vector with column names for horizon top, bottom depths}
    \item{\code{metadata}:}{object of class \code{"data.frame"}; metadata table}
    \item{\code{horizons}:}{object of class \code{"data.frame"}; table containing observations at different depths}
    \item{\code{site}:}{object of class \code{"data.frame"}; table containing observations at site locations}
    \item{\code{sp}:}{object of class \code{"SpatialPoints"}; locations of profiles}
    \item{\code{diagnostic}:}{object of class \code{"data.frame"}; table containing diagnostic properties}
}
Data of class \code{"FAO.SoilProfileCollection"} must satisfy all of the following requirements (class validity):
  \itemize{
   \item All variable names must be registered in the \href{https://www.google.com/fusiontables/DataSource?docid=1VVfBtsRijz_cRM8p5nj5hnijx6AQB9hbYzfRIZQ}{Global Soil Data Registry};
   \item All variable domains must correspond to the FAO Guidelines (2006 or later) for soil description or similar;
   \item All values must pass the validity checks i.e. numeric values must be within physical limits defined in the SoilGrids Global Soil Data Registry;
  }
}
\references{
\itemize{
  \item Beaudette, D. E., Roudier, P., & O'Geen, A. T. (2013) \href{http://dx.doi.org/10.1016/j.cageo.2012.10.020}{Algorithms for quantitative pedology: A toolkit for soil scientists}. Computers & Geosciences, 52, 258-268.
  \item FAO (2006) Guidelines for Soil Description. Food and Agriculture Organization of the United Nations, 4th Ed.
}
}
\author{ Tomislav Hengl }
\seealso{ \code{\link{SoilGrids-class}}, \code{\link{SpatialComponents-class}}, \code{\link{geosamples-class}} 
}
\examples{
library(aqp)
library(sp)

LONWGS84 = 3.90
LATWGS84 = 7.50 
UHDICM = 0
LHDICM = 30
SOURCEID = "ISRIC:NG0017"
SOURCEDB = "AfSP DB"
SPDFAO = "3"
TEXMHT = "SCL"
DCOMNS = "7.5YR_3_2"

sp1 <- new("FAO.SoilProfileCollection", 
  depthcols=c('UHDICM','LHDICM'),
  metadata=soil.vars, 
  horizons=data.frame(SOURCEID, UHDICM, LHDICM, TEXMHT, DCOMNS),
  site=data.frame(SOURCEID, SPDFAO, SOURCEDB),
  sp=SpatialPoints(data.frame(LONWGS84, LATWGS84), 
     proj4string=CRS("+proj=longlat +datum=WGS84"))
)
str(sp1)
}
\keyword{classes}