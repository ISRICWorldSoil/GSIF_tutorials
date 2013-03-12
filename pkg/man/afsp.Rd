\name{afsp}
\docType{data}
\encoding{latin1}
\alias{afsp}
\title{Africa Soil Profiles Database}
\description{A subset of the \href{http://africasoils.net/data/legacyprofile}{Africa Soil Profiles Database} (AFSP) that contains over 12,000 geo-referenced legacy soil profile records for 37 countries.}
\usage{data(afsp)}
\format{
The \code{afsp} data set contains two data frames --- sites and horizons. Sites table contains the following columns:
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
Horizons table contains the following columns:
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
  \item{\code{CEC}}{numeric; Cation Exchange Capacity in percent (\code{CecSoil} in the AFSP)}
}
The \code{afgrid20} \code{"SpatialPixelsDataFrame"} object contains a list of soil covariates prepared for Africa at 20 km resolution (45,997 grid nodes) obtained from WorldGrids.org:
  \describe{
  \item{\code{LMTGSH0}}{integer; land mask based on \href{http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html}{GSHHS}}
  \item{\code{SMKMOD0}}{integer; mask map indicating soil productive areas}
  \item{\code{DEMSRE0}}{integer; global relief model based on SRTM 30+ and ETOPO DEM at 1/120 arc degrees}
  \item{\code{SLPSRT0}}{integer; slope map in percent derived using the DEMSRE}
  \item{\code{TWISRE0}}{numeric; SAGA Topographic Wetness Index based on the SRTM DEM}
  \item{\code{G01ESA0}}{integer; percentage of land cover class 1 based on the \href{http://ionia1.esrin.esa.int}{GlobCov v2.2}}
  \item{\code{\dots}}{other land cover classes}
  \item{\code{IFLGRE0}}{integer; intact forest landscapes based on Greenpeace}
  \item{\code{L3POBI0}}{factor; physiographic landform units (\href{http://dx.doi.org/10.1016/j.geomorph.2011.12.001}{SCALA project})}
  \item{\code{PREGSM0}}{numeric; mean monthly precipitation (annual)}
  \item{\code{TDHMOD0}}{integer; maximum value of the 8-day MODIS day-time LST time series data}
  \item{\code{TDMMOD0}}{integer; mean value the 8-day MODIS day-time LST time series data}
  \item{\code{TDSMOD0}}{integer; standard deviation of the 8-day MODIS day-time LST time series data}
  \item{\code{TNMMOD0}}{integer; mean value the 8-day MODIS night-time LST time series data}
  \item{\code{TNSMOD0}}{integer; standard deviation of the 8-day MODIS night-time LST time series data}
  \item{\code{EVMMOD0}}{numeric; mean value of the monthly MODIS EVI time series data}
  \item{\code{EVSMOD0}}{integer; standard deviation of the monthly MODIS EVI time series data}
  \item{\code{SGEUSG0}}{factor; parent material class based on the \href{http://pubs.usgs.gov/of/1997/ofr-97-470/OF97-470A/}{USGS Surface Geology map}}      
}
}
\note{The soil covariate layers were selected to represent the standard soil forming factors: climate, organism, relief, and parent materials (CLORPT). Listed gridded layers follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{DEM} (digital elevation model); the next three letters represent the data source or collection method e.g. \code{SRE} (SRTM/ETOPO); the 6th character is the effective scale e.g. \code{0} indicates grid size of i.e. 1/5 decimal degrees (about 20 km).}
\author{The Africa Soil Profiles Database have been prepared by Johan Leenaars <johan.leenaars@wur.nl>. This is a subset of the original database that can be downloaded via \href{http://www.isric.org/sites/default/files/private/datasets/AfSP01Qry_ISRIC.zip}{www.isric.org}. The processing steps used to prepare this data are described in \href{http://gsif.r-forge.r-project.org/rw_AFSP.R}{this script}.}
\references{
\itemize{
\item Leenaars, J.G.B. (2012) \href{http://www.isric.org/content/publications}{Africa Soil Profiles Database, Version 1.0. A compilation of geo-referenced and standardized legacy soil profile data for Sub Saharan Africa (with dataset)}. ISRIC report 2012/03. Africa Soil Information Service (AfSIS) project and ISRIC --- World Soil Information, Wageningen, the Netherlands.
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
proj4string(sites) <- "+proj=latlong +datum=WGS84"
## obtain country borders:
library(maps)
country.m = map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
proj4string(country) = "+proj=latlong +datum=WGS84" 
## overlay and plot points and maps:
plot(country, col="darkgrey", xlim=c(-25.3,57.8), ylim=c(-34.8, 37.4))
points(sites, pch=21, bg="white", cex=.6, col="black")
## load gridded data for whole of Africa:
download.file("http://worldgrids.org/rda/afgrid20.rda", "afgrid20.rda")
load("afgrid20.rda")
## plot maps:
col = grey(rev(seq(0.1,0.975,0.025)))
names(afgrid20)
lns = list("sp.lines", spTransform(country, afgrid20@proj4string))
spplot(afgrid20["TDMMOD0a"], col.regions=col, sp.layout=lns)
## write to a GIS format:
writeGDAL(afgrid20["TDMMOD0a"], "TDMMOD0a.tif", "GTiff", mvFlag=-99999)
}
}
\keyword{datasets}
