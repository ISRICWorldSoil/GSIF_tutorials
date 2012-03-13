\name{NASIS}
\docType{data}
\alias{NASIS}
\title{NASIS pedon data for Indiana State}
\description{A subset of the \href{http://soils.usda.gov/technical/nasis/}{National Soil Information System} (NASIS) pedon data covering the Indiana State. Contains a total of 489 representative soil profiles.}
\usage{data(NASIS)}
\format{
Sites table. Contains the following columns:
  \describe{
  \item{\code{siteiid}}{factor; unique name that identifies a particular site}
  \item{\code{obsdate}}{POSIXct; sampling date}
  \item{\code{usiteid}}{factor; unique name set by the surveyor e.g. "92IN025001"}
  \item{\code{LON}}{numeric; longitude in decimal degrees (NAD83 datum)}
  \item{\code{LAT}}{numeric; latitude in decimal degrees (NAD83 datum)}
  \item{\code{pedonpurpose}}{factor; intended purpose of the profile description e.g. "full pedon description"}
  \item{\code{pedontype}}{factor; the description represents in relation to a series, component, etc e.g. "modal pedon for series"}
  \item{\code{geomposhill}}{factor; geomorphic component - hills e.g. "interfluve"}
  \item{\code{hillslopeprof}}{factor; two dimensional slope segments of a hillslope that have similar geometric, erosional, or depositional characteristics e.g. "footslope"}
  \item{\code{shapeacross}}{factor; two dimensional shape profile e.g. "concave"}
  \item{\code{shapedown}}{factor; the longitudinal shape of the slope e.g. "concave"}
  \item{\code{drainagecl}}{factor; the natural drainage conditions of the soil e.g. "very poorly"}
  \item{\code{siteperm}}{factor; site permeability e.g. "moderately slow"}
  \item{\code{runoff}}{factor; potential runoff class for the soil e.g. "very low"}
  \item{\code{soinmassamp}}{factor; soil name as sampled e.g. "Zanesville"}
  \item{\code{taxorder}}{factor; the Keys to Soil Taxonomy order class e.g. "alfisols"}
  \item{\code{taxsuborder}}{factor; the Keys to Soil Taxonomy suborder class e.g. "udalfs"}
  \item{\code{taxgrtgroup}}{factor; the Keys to Soil Taxonomy group class e.g. "fraglossudalfs"}
  \item{\code{taxsubgrp}}{factor; the Keys to Soil Taxonomy subgroup class e.g. "oxyaquic fraglossudalfs"}
  \item{\code{taxpartsize}}{factor; the Keys to Soil Taxonomy family / particle class e.g. "fine-silty"}
  \item{\code{taxtempregime}}{factor; the soil temperature regime e.g. "mesic"}
}
Horizons table. Contains the following columns:
  \describe{
  \item{\code{siteiid}}{factor; unique name that identifies a particular site}
  \item{\code{phiid}}{factor; unique name that identifies a particular horizon}
  \item{\code{hzname}}{factor; pedological horizon designation e.g. "Btx1"}
  \item{\code{seqnum}}{integer; horizon sequence number}
  \item{\code{hzdept}}{integer; horizon top depth in cm}
  \item{\code{hzdepb}}{integer; horizon bottom depth in cm}
  \item{\code{texture}}{factor; texture by hand class determined on the field e.g. "SIL"}
  \item{\code{rupresblkmst}}{integer; rupture moist class}
  \item{\code{phfield}}{numeric; pH as measured on the field}
  \item{\code{phdetermeth}}{factor; pH determination method}
  \item{\code{bounddistinct}}{integer; boundary distinctness grade}
  \item{\code{boundtopo}}{integer; horizontal shape of the interface between adjacent soil horizons}
}
}
\author{ The NASIS server is located at the USDA Kansas City data center. This subset has been prepared for the purpose of testing various DSM methods by Tomislav Hengl (tom.hengl@wur.nl). The processing steps are described in \href{http://gsif.r-forge.r-project.org/rw_NASIS.R}{this script}.}
\references{
\itemize{
\item National Soil Information System (\url{http://soils.usda.gov/technical/nasis/})
\item soilDB package (\url{http://cran.r-project.org/package=soilDB})
}
}
\note{The original NASIS database contains about 25 thousand profiles. Locations of samples were selected by experienced surveyors and are often clustered in some parts of the state (see plot below).\cr 
The variables listed above were selected to test various Digital Soil Mapping algorithms and represent only a subset of the NASIS data. Please refer to the NASIS technical documentation for more information about the specific columns.
}
\seealso{ \code{\link{NCSS}}, \code{\link{hoosier}}, \code{\link{geochm}}}
\examples{
library(aqp)
data(NASIS)
plot(NASIS$sites$LON, NASIS$sites$LAT)
aea.prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
NASIS.spc <- join(NASIS$horizons, NASIS$sites, type='inner')
depths(NASIS.spc) <- siteiid ~ hzdept + hzdepb
site(NASIS.spc) <- ~ LON + LAT + taxorder
# generate SpatialPoints:
coordinates(NASIS.spc) <- ~ LON + LAT
proj4string(NASIS.spc@sp) <- CRS("+proj=longlat +datum=NAD83")
# soil covariate layers:
data(indiana)
border <- SpatialLines(list(Lines(list(Line(indiana)), "border")))
proj4string(border) <- aea.prj
par(mar=c(0,0,0,0))
plot(spTransform(border, NASIS.spc@sp@proj4string))
points(NASIS.spc@sp, pch="+", cex=.7)
}
\keyword{datasets}
