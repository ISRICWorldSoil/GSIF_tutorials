\name{as.data.frame}
\docType{methods}
\alias{as.data.frame}
\alias{as.data.frame,SoilProfileCollection-method}
\alias{getHorizons}
\title{Converts an object of class \code{"SoilProfileCollection"} to a data frame}
\description{Converts an object of class \code{"SoilProfileCollection"} to an object of class \code{"data.frame"} with both site and horizon data sorted in one row. Each original column name in the horizons table receives a sufix \code{*_A, B, \dots, Z} where alphabetic letters represent horizon sequence. }
\usage{
\S4method{as.data.frame}{SoilProfileCollection}(x, row.names = NULL, optional = FALSE, \dots)
}
\arguments{
  \item{x}{object of class \code{"SoilProfileCollection"}}
  \item{row.names}{character; giving the row names for the data frame (missing values are not allowed)}
  \item{optional}{logical; if `TRUE', setting row names and converting column names (to syntactic names: see `make.names') is optional}
  \item{\dots}{optional arguments}                                   
}
\details{The advantage of converting the \code{"SoilProfileCollection"} data to a single table is that, once both tables have been merged to a single data frame, it can be more easily exported and visualized in a GIS and/or imported into a data base.}
\note{Few profiles with a large number of horizons can make the whole data frame become large. Consider removing such locations or aggregating measured values per horizon to a lower number of horizons.}
\author{ Tomislav Hengl and Brendan Malone }
\seealso{ \code{\link{as.geosamples}}, \code{\link{mpspline}} }
\examples{
library(aqp)
library(plyr)
library(rgdal)
library(sp)
## sample profile from Nigeria:
lon = 3.90; lat = 7.50; id = "ISRIC:NG0017"; FAO1988 = "LXp" 
top = c(0, 18, 36, 65, 87, 127) 
bottom = c(18, 36, 65, 87, 127, 181)
ORCDRC = c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
## prepare a SoilProfileCollection:
prof1 <- join(data.frame(id, top, bottom, ORCDRC), 
    data.frame(id, lon, lat, FAO1988), 
    type='inner')
depths(prof1) <- id ~ top + bottom
site(prof1) <- ~ lon + lat + FAO1988 
coordinates(prof1) <- ~ lon + lat
proj4string(prof1) <- CRS("+proj=longlat +datum=WGS84")
## convert to a simple table:
x <- as.data.frame(prof1)
str(x)
## horizons only
horizons <- getHorizons(x, idcol="id", sel=c("top", "bottom", "ORCDRC"))
horizons
}
\keyword{methods}
