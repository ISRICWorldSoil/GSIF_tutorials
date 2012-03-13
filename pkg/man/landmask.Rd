\name{landmask}
\docType{data}
\alias{landmask}
\title{1-degree cell land mask for Global soil mapping purposes}
\description{Land mask showing the 1-degree cells (about 19 thousand in total) in the geographical coordinates. Derived using the \href{http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html}{Global Self-consistent, Hierarchical, High-resolution Shoreline Database} data (GSHHS 2.1).}
\usage{data(landmask)}
\format{
Data frame with the following columns:
  \describe{
  \item{\code{mask}}{boolean; land mask value}
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
\item Global Self-consistent, Hierarchical, High-resolution Shoreline Database (\url{http://en.wikipedia.org/wiki/GSHHS})
}
}
\note{The land mask has been generated using the layer \code{GSHHS_shp/h/GSHHS_h_L1.shp} (level-1 boundaries). For practical purposes, all latitudes beyond 65S and 84N (polar land masses) have been masked out.
}
\examples{
data(landmask)
gridded(landmask) <- ~x+y
proj4string(landmask) <- "+proj=longlat +datum=WGS84"
\dontrun{
spplot(landmask["mask"], col.regions="grey")
}
}
\keyword{datasets}
