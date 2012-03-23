\name{landmask}
\docType{data}
\alias{landmask}
\title{1-degree cell land mask for Global soil mapping purposes}
\description{Land mask showing the 1-degree cells (about 19 thousand in total) in the geographical coordinates, and the productive soils mask (areas with a positive Leaf Area Index at least once in the period 2002 till today). The land mask is based on the \href{http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html}{Global Self-consistent, Hierarchical, High-resolution Shoreline Database} data (GSHHS 2.1), and the productive soils mask on the MODIS Leaf Area Index monthtly product (\href{https://lpdaac.usgs.gov/lpdaac/products/modis_products_table/leaf_area_index_fraction_of_photosynthetically_active_radiation/8_day_l4_global_1km/mcd15a2}{MOD15A2}).}
\usage{data(landmask)}
\format{
Data frame with the following columns:
  \describe{
  \item{\code{mask}}{boolean; land mask value}
  \item{\code{soilmask}}{boolean; soil mask value}
	\item{\code{Lon_it}}{indication of the longitude quadrant (W or E)}
  \item{\code{Lat_it}}{indication of the latitude quadrant (S or N)}
	\item{\code{cell_id}}{cell id code e.g. \code{W79_N83}}
	\item{\code{x}}{longitudes of the center of the grid nodes}
	\item{\code{y}}{latitudes of the center of the grid nodes}  	
}
}
\author{ Tomislav Hengl }
\references{
\itemize{
\item Wessel, P., Smith, W.H.F., (1996) \href{http://dx.doi.org/10.1029/96JB00104}{A Global Self-consistent, Hierarchical, High-resolution Shoreline Database}. Journal of Geophysical Research, 101, 8741-8743.
\item Savtchenko, A., D. Ouzounov, S. Ahmad, J. Acker, G. Leptoukh, J. Koziana, and D. Nickless, (2004) \href{http://dx.doi.org/10.1016/j.asr.2004.03.012}{Terra and Aqua MODIS products available from NASA GES DAAC}. Advances in Space Research 34(4), 710–714.
\item Global Self-consistent, Hierarchical, High-resolution Shoreline Database (\url{http://en.wikipedia.org/wiki/GSHHS})
}
}
\note{The land mask has been generated from the layer \code{GSHHS_shp/h/GSHHS_h_L1.shp} (level-1 boundaries). For practical purposes, all latitudes beyond 65S and 84N (polar land masses) have been masked out.
}
\examples{
data(landmask)
gridded(landmask) <- ~x+y
proj4string(landmask) <- "+proj=longlat +datum=WGS84"
\dontrun{
spplot(landmask["mask"], col.regions="grey")
spplot(landmask["soilmask"], col.regions="grey")
}
}
\keyword{datasets}
