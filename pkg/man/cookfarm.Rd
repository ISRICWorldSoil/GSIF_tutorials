\name{cookfarm}
\docType{data}
\encoding{latin1}
\alias{cookfarm}
\title{The Cook Agronomy Farm data set}
\description{The R.J. Cook Agronomy Farm (\code{cookfarm}) is a Long-Term Agroecosystem Research Site operated by Washington State University, located near Pullman, Washington, USA. Contains spatio-temporal (3D+T) measurements of three soil properties and a number of spatial and temporal regression covariates.}
\usage{data(cookfarm)}
\format{
The \code{cookfarm} data set contains four data frames. The \code{readings} data frame contains measurements of volumetric water content (cubic-m/cubic-m), temperature (degree C) and bulk electrical conductivity (dS/m), measured at 42 locations using 5TE sensors at five standard depths (0.3, 0.6, 0.9, 1.2, 1.5 m) for the period "2011-01-01" to "2012-12-31":
\describe{
  \item{\code{SOURCEID}}{factor; unique station ID}
  \item{\code{Date}}{date; observation day}
  \item{\code{Port*VW}}{numeric; volumetric water content measurements at five depths}
  \item{\code{Port*C}}{numeric; soil temperature measurements at five depths}
  \item{\code{Port*EC}}{numeric; bulk electrical conductivity measurements at five depths}
}
The \code{profiles} data frame contains soil profile descriptions from 142 sites:
\describe{
  \item{\code{SOURCEID}}{factor; unique station ID}
  \item{\code{Easting}}{numeric; x coordinate in the local projection system}
  \item{\code{Northing}}{numeric; y coordinate in the local projection system}
  \item{\code{TAXNUSDA}}{factor; Keys to Soil Taxonomy taxon name e.g. \code{"Caldwell"}}
  \item{\code{HZDUSD}}{factor; horizon designation}
  \item{\code{UHDICM}}{numeric; upper horizon depth from the surface in cm}
  \item{\code{LHDICM}}{numeric; lower horizon depth from the surface in cm}
  \item{\code{BLD}}{bulk density in tonnes per cubic-meter}
  \item{\code{PHIHOX}}{numeric; pH index measured in water solution}
}
The \code{grids} data frame contains values of regression covariates at 10 m resolution:
\describe{
  \item{\code{DEM}}{numeric; Digital Elevation Model}
  \item{\code{TWI}}{numeric; SAGA GIS Topographic Wetness Index}
  \item{\code{MUSYM}}{factor; soil mapping units e.g. \code{"Thatuna silt loam"}}
  \item{\code{NDRE.M}}{numeric; mean value of the Normalized Difference Red Edge Index (time series of 11 RapidEye images)}
  \item{\code{NDRE.sd}}{numeric; standard deviation of the Normalized Difference Red Edge Index (time series of 11 RapidEye images)}
  \item{\code{Cook_fall_ECa}}{numeric; apparent electrical conductivity image from fall}
  \item{\code{Cook_spr_ECa}}{numeric; apparent electrical conductivity image from spring}
  \item{\code{X2011}}{factor; cropping system in 2011}  
  \item{\code{X2012}}{factor; cropping system in 2012}
}
The \code{weather} data frame contains daily temperatures and rainfall from the nearest meteorological station:
\describe{
  \item{\code{Date}}{date; observation day}
  \item{\code{Precip_wrcc}}{numeric; observed precipitation in mm}
  \item{\code{MaxT_wrcc}}{numeric; observed maximum daily temperature in degree C}
  \item{\code{MinT_wrccc}}{numeric; observed minimum daily temperature in degree C}
}
}
\note{The farm is 37 ha, stationed in the hilly Palouse region, which receives an annual average of 550 mm of precipitation, primarily as rain and snow in November through May. Soils are deep silt loams formed on loess hills; clay silt loam horizons commonly occur at variable depths. Farming practices at Cook Farm are representative of regional dryland annual cropping systems (direct-seeded cereal grains and legume crops).}
\author{Caley Gasch, Tomislav Hengl and David J. Brown}
\references{
\itemize{
\item Gasch, C., Hengl, T., \enc{Gräler}{Graeler}, B., Meyer, H., Magney, T., Brown, D.J., 2015. Spatio-temporal interpolation of soil water, temperature, and electrical conductivity in 3D+T: the Cook Agronomy Farm data set. Spatial Statistics Journal, accepted. 
}
}
\examples{
## An example for 3D+T modelling applied to the cookfarm data set can be assesed via 
## demo(cookfarm_3DT_kriging)
## demo(cookfarm_3DT_RF)
## Please note that the demo's might take 10-15 minutes to complete.
library(rgdal)
library(sp)
library(spacetime)
library(aqp)
library(splines)
library(randomForest)
library(plyr)
library(plotKML)
data(cookfarm)

## gridded data:
grid10m <- cookfarm$grids
gridded(grid10m) <- ~x+y
proj4string(grid10m) <- CRS(cookfarm$proj4string)
spplot(grid10m["DEM"], col.regions=SAGA_pal[[1]])

## soil profiles:
profs <- cookfarm$profiles
levels(cookfarm$profiles$HZDUSD)
## Bt horizon:
sel.Bt <- grep("Bt", profs$HZDUSD, ignore.case=FALSE, fixed=FALSE)
profs$Bt <- 0
profs$Bt[sel.Bt] <- 1
depths(profs) <- SOURCEID ~ UHDICM + LHDICM
site(profs) <- ~ TAXSUSDA + Easting + Northing
coordinates(profs) <- ~Easting + Northing
proj4string(profs) <- CRS(cookfarm$proj4string)
profs.geo <- as.geosamples(profs)

## fit model for Bt horizon:
m.Bt <- GSIF::fit.gstatModel(profs.geo, Bt~DEM+TWI+MUSYM+Cook_fall_ECa
   +Cook_spr_ECa+ns(altitude, df = 4), grid10m, fit.family = binomial(logit))
plot(m.Bt)

## fit model for soil pH:
m.PHI <- fit.gstatModel(profs.geo, PHIHOX~DEM+TWI+MUSYM+Cook_fall_ECa
    +Cook_spr_ECa+ns(altitude, df = 4), grid10m)
plot(m.PHI)
}
\keyword{datasets}
