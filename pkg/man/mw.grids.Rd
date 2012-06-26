\name{mw.grids}
\docType{data}
\alias{mw.grids}
\title{Malawi study area grids at 1/120 arcdegrees (1 km) resolution}
\description{A list of soil covariates prepared for Malawi (94,215 grid nodes).}
\usage{str(mw.grids)}
\format{
Data frame contains the following columns:
  \describe{
  \item{\code{DEMSRT3}}{numeric; SRTM DEM}  
  \item{\code{EV1MOD3}}{numeric; first principal component of the monthly MODIS EVI time series data}
  \item{\code{EV2MOD3}}{numeric; second principal component of the monthly MODIS EVI time series data}
  \item{\code{GLCESA3}}{factor; general land cover class based on the \href{http://ionia1.esrin.esa.int}{GlobCov v2.2}}
  \item{\code{INSSRT3}}{numeric; potential incoming solar radiation derived in SAGA GIS}
  \item{\code{SLPSRT3}}{numeric; slope in degrees based on the SRTM DEM}
  \item{\code{TAXWRB3}}{factor; dominant soil class (group) following the World Reference Base soil classification system derived from the general soil polygon map of Malawi}
  \item{\code{TD1MOD3}}{numeric; first principal component of the 8-day MODIS day-time LST time series data}         
  \item{\code{TN1MOD3}}{numeric; first principal component of the 8-day MODIS night-time LST time series data} 
  \item{\code{TWISRT2}}{numeric; SAGA Topographic Wetness Index based on the SRTM DEM}
  \item{\code{x}}{numeric; x-coordinate in the UTM 36S coordiante system with WGS84 ellipsoid}
  \item{\code{y}}{numeric; y-coordinate in the UTM 36S coordiante system with WGS84 ellipsoid}
}
}
\details{The landscape of Malawi is a series of dissected plateaus and mountainous areas along the Rift Valley i.e. lake Malawi. The size of the country is 118,000 square kilometers (the size of the rectangular block is 349 by 860 km). Dominat soil type map (\code{TAXWRB3}) was derived from the 1:600k scale soils map of Malawi and rasterized to 1 km resolution.  }
\note{Listed layer follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{DEM} (digital elevation model); the next three letters represent the data source or collection method e.g. \code{SRT} (SRTM); the 6th character is the effective scale e.g. \code{3} indicates the 2nd standard scale i.e. 1/120 decimal degrees (about 1 km).}
\author{ Tomislav Hengl & Markus Walsh}
\references{
\itemize{
\item Ministry of Agriculture Malawi / FAO, (1992) Soils map of Malawi. Land Resources Evaluation Project /MLW/85/011. Land Resources Conservation Board, Lilongwe.
\item Africa Soil Information Service (\url{http://africasoils.net})
}
}
\examples{
\dontrun{# Load the data from GSIF page:
download.file("http://worldgrids.org/rda/mw.grids.rda", "mw.grids.rda")
load("mw.grids.rda")
gridded(mw.grids) <- ~x+y
proj4string(mw.grids) <- CRS(utm.csy)
# plot maps:
col = grey(rev(seq(0.1,0.975,0.025)))
names(mw.grids)
spplot(mw.grids[1], col.regions=col)
# write to a GIS format:
writeGDAL(mw.grids["EV1MOD3"], "EV1MOD3.tif", "GTiff", mvFlag=-99999)
}
}
\keyword{datasets}