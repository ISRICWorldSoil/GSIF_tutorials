\name{ilin.grids}
\docType{data}
\alias{ilin.grids}
\title{Illinois and Indiana State soil covariates at 3/120 arcdegrees (2.5 km) resolution}
\description{A list of soil covariates prepared for the Illinois and Indiana States (40383 grid nodes).}
\usage{str(ilin.grids)}
\format{
Data frame contains the following columns:
  \describe{
  \item{\code{CRPUSG2}}{numeric; percentage (0--1) of crop fields based on the National Land cover map from 2006}
  \item{\code{FORUSG2}}{numeric; percentage (0--1) of deciduous forest based on the National Land cover map from 2006}
  \item{\code{DINSRT2}}{numeric; diffuse insolation derived in SAGA GIS}  
  \item{\code{INSSRT2}}{numeric; potential incoming solar radiation derived in SAGA GIS}
  \item{\code{LNDUSG2}}{numeric; land mask (0--1)}
  \item{\code{TRITOX2}}{numeric; kernel density of the Toxics Release Inventory (EPA) locations}
  \item{\code{TWISRT2}}{numeric; SAGA Topographic Wetness Index based on the SRTM DEM}
  \item{\code{SLPSRT2}}{numeric; slope in degrees based on the SRTM DEM}
  \item{\code{SORGSM2}}{factor; dominant suborder (Keys to Soil Taxonomy) class from the General Soil Map STATSGO2 (at scale 1:250,000) obtained from the \href{http://datagateway.nrcs.usda.gov}{NRCS data gateway}}         
  \item{\code{LON}}{numeric; x-coordinate (decimal degrees) in the WGS84 system}
  \item{\code{LAT}}{numeric; y-coordinate (decimal degrees) in the WGS84 system}
}
}
\details{ Majority of the gridded layers listed were obtained from the \href{http://seamless.usgs.gov/}{USGS Seamless data warehouse}. The percentage of crop fields and the land mask (0--100 percent) were derived using the area weighted averaging from the 30 m resolution National Land cover map from 2006. The land mask (0--100 percent) shows percentage of water bodies, medium and high intensity build-up (urban) areas. To create a boolean mask map, consider using some treshold value (e.g. 0.6). Density of the industrial facilities / pollutants (Toxics Release Inventory) were obtained from the \href{http://toxmap.nlm.nih.gov/toxmap/}{U.S. National Library of Medicine}. Land surface parameters diffuse insolation, potential incoming solar radiation, SAGA Topographic Wetness Index and slope map were derived in SAGA GIS using the 250 m resolution SRTM DEM obtained from the Consultative Group on International Agricultural Research (CGIAR), then resampled to 3/120 arcdegrees (2.5 km resolution). Soil suborder layer was derived as dominant taxonomic components (highest percentage) in the General Soil Map STATSGO2 (at scale 1:250,000). The processing steps used to prepare the map of dominatant taxonomic compenents are described in \href{http://gsif.r-forge.r-project.org/gsmsoilmu_a_us.R}{this script}. There are cases where the dominant suborder is less than 30\% of the map unit, which means that selecting a single suborer is not likely all that informative. To obtain additional covariate layers, e.g. SRTM DEM, MODIS, GlobCov and similar publicly available global data sets refer to the WorldGrids.org repository. }
\note{Listed layer follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{CRP} (percentage (0--1) of crop fields); the next three letters represent the data source or collection method e.g. \code{USG} (USGS); the 6th character is the effective scale e.g. \code{2} indicates the 3rd standard scale i.e. 3/120 decimal degrees (about 2.5 km).}
\author{ Tomislav Hengl}
\references{
\itemize{
\item Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. U.S. General Soil Map (STATSGO2). Available online at \url{http://soildatamart.nrcs.usda.gov}.
\item Xian, G., and Homer C., (2010) \href{http://dx.doi.org/10.1016/j.rse.2010.02.018}{Updating the 2001 National Land Cover Database impervious surface products to 2006 using Landsat imagery change detection methods}. Remote Sensing of Environment, 114, 1676-1686.
\item soilDB package (\url{http://cran.r-project.org/package=soilDB})
\item WorldGrids repository of soil covariate layers (\url{http://worldgrids.org})
}
}
\seealso{ \code{\link{ilin}}, \code{\link{geochm}}}
\examples{
\dontrun{# Load the data from GSIF page:
download.file("http://gsif.r-forge.r-project.org/ilin.grids.rda", "ilin.grids.rda")
load("ilin.grids.rda")
gridded(ilin.grids) <- ~LON+LAT
proj4string(ilin.grids) <- CRS("+proj=longlat +datum=WGS84")
# obtain state borders:
require(maps)
states.m = map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(states.m$names, ":"), function(x) x[1])
require(maptools)
states <- as(map2SpatialPolygons(states.m, IDs=IDs), "SpatialLines")
# plot some layers:
lst = c("CRPUSG2","FORUSG2")
col = grey(rev(seq(0.1,0.975,0.025)))
spplot(ilin.grids[lst], col.regions=col, sp.layout=list("sp.lines", states))
spplot(ilin.grids["SORGSM2"], sp.layout=list("sp.lines", states))
# write to a GIS format:
writeGDAL(ilin.grids["FORUSG2"], "FORUSG2.tif", "GTiff", mvFlag=-99999)
}
}
\keyword{datasets}