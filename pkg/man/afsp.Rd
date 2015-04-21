\name{afsp}
\docType{data}
\encoding{latin1}
\alias{afsp}
\title{Africa Soil Profiles Database}
\description{A merge of the \href{http://africasoils.net/data/legacyprofile}{Africa Soil Profiles Database} (AFSP) with 17,000+ geo-referenced legacy soil profile records, and AfSIS Sentinel Site database with 9000+ sampling locations.}
\usage{data(afsp)}
\format{
The \code{afsp} data set contains two data frames --- \code{sites} and \code{horizons}. Sites table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique label to help a user identify a particular site (\code{ProfileID} in the AFSP)}
  \item{\code{SOURCEDB}}{factor; source data base}
  \item{\code{LONWGS84}}{numeric; longitude in decimal degrees on the WGS84 datum (\code{X_LonDD} in the AFSP)}
  \item{\code{LATWGS84}}{numeric; latitude in decimal degrees on the WGS84 datum (\code{Y_LatDD} in the AFSP)}
  \item{\code{TIMESTRR}}{character; the date on which this particular soil was described or sampled (\code{T_Year} in the AFSP)}
  \item{\code{TAXGWRB}}{factor; abbreviated soil group based on the WRB classification system (\code{WRB06rg} in the AFSP)}
  \item{\code{TAXNUSDA}}{factor; Keys to Soil Taxonomy taxon name e.g. \code{"Plinthic Udoxic Dystropept"} (\code{USDA} in the AFSP)}
  \item{\code{BDRICM}}{numeric; depth to bedrock in cm}
  \item{\code{DRAINFAO}}{factor; drainage class based on the FAO guidelines for soil description: \code{E} (excessively drained), \code{S} (somewhat excessively drained), \code{W} (well drained), \code{M} (moderately well drained), \code{I} (somewhat poorly drained) and \code{V} (very poorly drained)}
}
Horizons table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; a short label to help a user identify a particular site}
  \item{\code{UHDICM}}{numeric; upper horizon depth from the surface in cm}
  \item{\code{LHDICM}}{numeric; lower horizon depth from the surface in cm}
  \item{\code{MCOMNS}}{factor; Munsell color moist}
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permilles}
  \item{\code{PHIHOX}}{numeric; pH index measured in water solution}
  \item{\code{SNDPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm)}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm)}
  \item{\code{CRFVOL}}{numeric; volume percentage of coarse fragments (> 2 mm)}
  \item{\code{BLD}}{numeric;bulk density in tonnes per cubic-meter}
  \item{\code{CEC}}{numeric; Cation exchange capacity (fine earth fraction) in cmolc/kg}
  \item{\code{NTO}}{numeric; total N content in permille or g/kg}
  \item{\code{EMGX}}{numeric; exchangable Mg in cmolc/kg}
}
}
\author{The Africa Soil Profiles Database have been prepared by Johan Leenaars <johan.leenaars@wur.nl>. This is a subset of the original database that can be downloaded via \href{http://www.isric.org/content/africa-soil-profiles-database}{www.isric.org}. The AfSIS Sentinel Site database is one of the main deliverables of the Africa Soil Information Service project.}
\references{
\itemize{
\item Leenaars, J.G.B. (2014) \href{http://www.isric.org/content/publications}{Africa Soil Profiles Database, Version 1.2. A compilation of geo-referenced and standardized legacy soil profile data for Sub Saharan Africa (with dataset)}. ISRIC report 2012/03. Africa Soil Information Service (AfSIS) project and ISRIC --- World Soil Information, Wageningen, the Netherlands.
\item Africa Soil Information Service (\url{http://africasoils.net}) 
}
}
\examples{
\dontrun{
library(rgdal)
library(aqp)
library(sp)

data(afsp)
sites <- afsp$sites
coordinates(sites) <- ~ LONWGS84 + LATWGS84
proj4string(sites) <- "+proj=longlat +datum=WGS84"
## obtain country borders:
library(maps)
country.m = map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
proj4string(country) = "+proj=longlat +datum=WGS84" 
## overlay and plot points and maps:
plot(country, col="darkgrey", xlim=c(-25.3,57.8), ylim=c(-34.8, 37.4))
points(sites, pch=21, bg="white", cex=.6, col="black")
}
}
\keyword{datasets}
