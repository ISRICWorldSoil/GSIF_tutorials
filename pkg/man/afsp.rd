\name{afsp}
\docType{data}
\alias{afsp}
\title{Africa Soil Profiles Database}
\description{A subset of the \href{http://africasoils.net/data/legacyprofile}{Africa Soil Profiles Database} (AFSP) that contains over 12,000 geo-referenced legacy soil profile records for 37 countries.}
\usage{data(afsp)}
\format{
Sites table. Contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique label to help a user identify a particular site (\code{ProfileID} in the AFSP)}
  \item{\code{RLBISRIC}}{factor; reliability class according to the ISRIC quality standards (\code{Reliab} in the AFSP)}
  \item{\code{LONWGS84}}{numeric; longitude in decimal degrees on the WGS84 datum (\code{X_LonDD} in the AFSP)}
  \item{\code{LATWGS84}}{numeric; latitude in decimal degrees on the WGS84 datum (\code{Y_LatDD} in the AFSP)}
  \item{\code{STDXYZ}}{numeric; profile location accuracy in decimal degrees (\code{XYAccur} in the AFSP)}
  \item{\code{TIMESTRR}}{character; the date on which this particular soil was described or sampled (\code{T_Year} in the AFSP)}
  \item{\code{TAXNWRB}}{factor; full soil classification name based on the WRB classification system (\code{ProfileID} in the AFSP)}
  \item{\code{WRB06}}{factor; abbreviated soil group based on the WRB classification system (\code{WRB06rg} in the AFSP)}
  \item{\code{TAXNUSDA}}{factor; Keys to Soil Taxonomy taxon name e.g. \code{"Plinthic Udoxic Dystropept"} (\code{USDA} in the AFSP)}
}
Horizons table. Contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; a short label to help a user identify a particular site (\code{ProfileID} in the AFSP)}
  \item{\code{LSQINT}}{integer; a layer sequence number 1 to N (\code{LayerNr} in the AFSP)}
  \item{\code{UHDICM}}{numeric; upper horizon depth from the surface in cm (\code{UpDpth} in the AFSP)}
  \item{\code{LHDICM}}{numeric; lower horizon depth from the surface in cm (\code{LowDpth} in the AFSP)}
  \item{\code{MCOMNS}}{factor; Munsell color moist (\code{ColorM} in the NSCD)}
  \item{\code{CRFVOL}}{numeric; volume percentage of coarse fragments (> 2 mm; \code{CfPc} in the AFSP)}
  \item{\code{SNDPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm; \code{Sand} in the AFSP)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm; \code{Silt} in the AFSP)}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm; \code{Clay} in the AFSP)}
  \item{\code{BLD}}{bulk density in kg per cubic meter (\code{BlkDens} in the AFSP)}
  \item{\code{PHIHO5}}{numeric; pH index measured in water solution(\code{PHH2O} in the AFSP)}
  \item{\code{PHIKCL}}{numeric; pH index measured in KCl solution (\code{PHKCl} in the AFSP)}
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permille (\code{OrgC} in the AFSP)}
}
}
\author{ The Africa Soil Profiles Database have been prepared by Johan Leenaars <johan.leenaars@wur.nl>. This is a subset of the original database that can be downloaded via \href{http://www.isric.org/sites/default/files/private/datasets/AfSP01Qry_ISRIC.zip}{www.isric.org}. The processing steps used to prepare this data are described in \href{http://gsif.r-forge.r-project.org/rw_AFSP.R}{this script}.}
\references{
\itemize{
\item Leenaars, J.G.B. (2012) \href{http://www.isric.org/content/publications}{Africa Soil Profiles Database, Version 1.0. A compilation of geo-referenced and standardized legacy soil profile data for Sub Saharan Africa (with dataset)}. ISRIC report 2012/03. Africa Soil Information Service (AfSIS) project and ISRIC – World Soil Information, Wageningen, the Netherlands.
\item Africa Soil Information Service (\url{http://africasoils.net}) 
}
}
\seealso{ \code{\link{mw.grids}}}
\examples{
\dontrun{# load data and convert to SPC:
data(afsp)
sites <- afsp$sites
coordinates(sites) <- ~ LONWGS84 + LATWGS84
proj4string(sites) <- "+proj=latlong +datum=WGS84"
# obtain country borders:
library(maps)
country.m = map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
# overlay and plot points and maps:
plot(country, col="darkgrey", xlim=c(-25.3,57.8), ylim=c(-34.8, 37.4))
points(sites, pch=21, bg="white", cex=.6, col="black")
}
}
\keyword{datasets}
