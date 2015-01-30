\name{landmask}
\docType{data}
\alias{landmask}
\alias{landmask20km}
\title{Global coarse resolution land / soil mask maps}
\description{Land mask showing the 1-degree cells (about 19 thousand in total) in the geographical coordinates, and the productive soils mask (areas with a positive Leaf Area Index at least once in the period 2002--2011). The land mask is based on the \href{http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html}{Global Self-consistent, Hierarchical, High-resolution Shoreline Database} data (GSHHS 2.1), the productive soils mask on the MODIS Leaf Area Index monthtly product (\href{https://lpdaac.usgs.gov/products/modis_products_table}{MOD15A2}), and the water mask is based on the \href{http://glcf.umd.edu/data/watermask/}{MOD44W} product. The map of the Keys to Soil Taxonomy soil suborders of the world at 20 km is based on the \href{http://www.nrcs.usda.gov/}{USDA-NRCS map of the global soil regions}.}
\usage{data(landmask)}
\format{
\code{landmask} data set is a data frame with the following columns:
  \describe{
  \item{\code{mask}}{percent; land mask value}
  \item{\code{soilmask}}{boolean; soil mask value}
  \item{\code{watermask}}{percent; water mask value}
	\item{\code{Lon_it}}{indication of the longitude quadrant (W or E)}
  \item{\code{Lat_it}}{indication of the latitude quadrant (S or N)}
	\item{\code{cell_id}}{cell id code e.g. \code{W79_N83}}
	\item{\code{x}}{longitudes of the center of the grid nodes}
	\item{\code{y}}{latitudes of the center of the grid nodes}  	
}
\code{landmask20km} data set is an object of class \code{SpatialGridDataFrame} with the following columns:
  \describe{
  \item{\code{mask}}{percent; land mask value}
  \item{\code{suborder}}{factor; Keys to Soil Taxonomy suborder class e.g. Histels, Udolls, Calcids, \dots}
  \item{\code{soilmask}}{factor; global soil mask map based on the land cover classes (see: \href{http://worldgrids.org/doku.php?id=wiki:smkisr}{SMKISR3})}
}
}
\references{
\itemize{
\item Carroll, M., Townshend, J., DiMiceli, C., Noojipady, P., Sohlberg, R. (2009) \href{http://dx.doi.org/10.1080/17538940902951401}{A New Global Raster Water Mask at 250 Meter Resolution}. International Journal of Digital Earth, 2(4).
\item Global Self-consistent, Hierarchical, High-resolution Shoreline Database (\url{http://en.wikipedia.org/wiki/GSHHS})
\item USDA-NRCS Global Soil Regions Map (\url{http://www.nrcs.usda.gov/})
\item Savtchenko, A., D. Ouzounov, S. Ahmad, J. Acker, G. Leptoukh, J. Koziana, and D. Nickless, (2004) \href{http://dx.doi.org/10.1016/j.asr.2004.03.012}{Terra and Aqua MODIS products available from NASA GES DAAC}. Advances in Space Research 34(4), 710-714.
\item Wessel, P., Smith, W.H.F., (1996) \href{http://dx.doi.org/10.1029/96JB00104}{A Global Self-consistent, Hierarchical, High-resolution Shoreline Database}. Journal of Geophysical Research, 101, 8741-8743.
}
}
\note{The land mask has been generated from the layer \code{GSHHS_shp/h/GSHHS_h_L1.shp} (level-1 boundaries).
}
\seealso{ \code{rworldmap::rworldmapExamples}, \code{maps::map}}
\examples{
library(rgdal)
library(sp)

data(landmask)
gridded(landmask) <- ~x+y
proj4string(landmask) <- "+proj=longlat +datum=WGS84"
\dontrun{ ## plot maps:
library(maps)
country.m = map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
library(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
spplot(landmask["mask"], col.regions="grey", sp.layout=list("sp.lines", country))
spplot(landmask["soilmask"], col.regions="grey", sp.layout=list("sp.lines", country))
}
## also available in the Robinson projection at 20 km grid:
data(landmask20km)
image(landmask20km[1])
summary(landmask20km$suborder)
summary(landmask20km$soilmask)
}
\keyword{datasets}
