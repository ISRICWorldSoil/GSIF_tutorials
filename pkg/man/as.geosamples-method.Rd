\name{as.geosamples}
\docType{methods}
\alias{as.geosamples}
\alias{as.geosamples,SoilProfileCollection-method}
\alias{as.geosamples,SpatialPointsDataFrame-method}
\title{Converts an object to geosamples class}
\description{Converts an object of class \code{"SoilProfileCollection"} or \code{"SpatialPointsDataFrame"} to an object of class \code{"geosamples"} with all measurements broken into individual records. Geosamples are standardized spatially and temporally referenced samples from the Earth's surface.}
\usage{
\S4method{as.geosamples}{SoilProfileCollection}(obj, 
    registry = as.character(NA), sample.area = 1, mxd = 2, TimeSpan.begin, TimeSpan.end)
\S4method{as.geosamples}{SpatialPointsDataFrame}(obj, 
    registry = as.character(NA), sample.area = 1, mxd = 2, TimeSpan.begin, TimeSpan.end)
}
\arguments{
  \item{obj}{object of class \code{"SoilProfileCollection"}}
  \item{\dots}{optional arguments}
  \item{registry}{URI specifying the metadata registry (web-service that carries all metadata connected to the certain method ID and/or sample ID)}
  \item{sample.area}{standard sample area in square meters (assumed to be 1 by 1 m)}
  \item{mxd}{maximum depth of interest in meters}
  \item{TimeSpan.begin}{vector of class \code{"POSIXct"}; begin of the measurement period}
  \item{TimeSpan.end}{vector of class \code{"POSIXct"}; end of the measurement period}  
}
\value{Returns an object of type \code{"geosamples"}. Many columns required by the \code{"geosamples"} class might be not available and will result in \code{NA} values. To ensure compatibility, when building an object of type \code{"SoilProfilesCollection"}, use some standard naming convention to attach attributes to each measurement (horizons and sites slots in the \code{"SoilProfileCollection-class"}):
\describe{
  \item{\code{"locationError"}}{can be used to attach location errors in meters to each spatial location}
  \item{\code{"sampleArea"}}{can be used to attach spatial support to each measurement (usually 1 by 1 meter)}
  \item{\code{"measurementError"}}{can be used to attach specific measurement errors to each measurement in both site and horizons table}
  \item{\code{"IGSN"}}{can be used to attach the unique identifier (\href{http://en.wikipedia.org/wiki/International_Geo_Sample_Number}{International Geo Sample Number}) to each specific observation (corresponds to the \code{"observationid"} column)}  
}
}
\author{ Tomislav Hengl and Hannes I. Reuter }
\seealso{ \code{\link{geosamples-class}}, \code{\link{as.data.frame}}, \code{aqp::SoilProfileCollection} }
\examples{
library(aqp)
library(plyr)
# sample profile from Nigeria:
lon = 3.90; lat = 7.50; time = as.POSIXct("1978", format="\%Y") 
id = "ISRIC:NG0017"; TAXNFAO8 = "LXp" 
top = c(0, 18, 36, 65, 87, 127) 
bottom = c(18, 36, 65, 87, 127, 181)
ORCDRC = c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
methodid = c("TAXNFAO8", "ORCDRC")
description = c("FAO 1988 classification system group", 
    "Method of Walkley-Black (Org. matter = Org. C x 1.72)")
units = c("FAO 1988 classes", "permille")
detectionLimit = c(as.character(NA), "0.1")
# prepare a SoilProfileCollection:
prof1 <- join(data.frame(id, top, bottom, ORCDRC), 
    data.frame(id, lon, lat, time, TAXNFAO8), type='inner')
depths(prof1) <- id ~ top + bottom
site(prof1) <- ~ lon + lat + time + TAXNFAO8 
coordinates(prof1) <- ~ lon + lat + time
proj4string(prof1) <- CRS("+proj=longlat +datum=WGS84")
# add measurement errors:
attr(prof1@horizons$ORCDRC, "measurementError") <- c(1.5, 0.5, 0.5, 0.5, 0.5, 0.5)
attr(prof1@sp@coords, "locationError") <- 1500
# add the metadata:
prof1@metadata <- data.frame(methodid, description, units, detectionLimit)
# convert to geosamples:
x <- as.geosamples(prof1)
x
# print only the sampled values of ORCDRC:
ORCDRC <- subset(x, "ORCDRC")
ORCDRC[,c("sampleid", "altitude", "observedValue")]

# convert object of type SpatialPointsDataFrame:
data(meuse)
# prepare columns:
names(meuse)[which(names(meuse)=="x")] = "longitude"
names(meuse)[which(names(meuse)=="y")] = "latitude"
meuse$altitude = -.15
meuse$time = unclass(as.POSIXct("1992-01-01"))
coordinates(meuse) <- ~ longitude + latitude + altitude + time
proj4string(meuse) <- CRS("+init=epsg:28992")
library(plotKML)
hm <- reproject(meuse[,c("zinc", "copper")])
hm.geo <- as.geosamples(hm)
hm.geo
}
\keyword{methods}
