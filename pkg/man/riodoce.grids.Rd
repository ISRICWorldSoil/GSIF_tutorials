\name{riodoce.grids}
\docType{data}
\encoding{latin1}
\alias{riodoce.grids}
\title{Rio Doce study area grids at 1/120 arcdegrees (1 km) resolution}
\description{A list of soil covariates prepared for the Rio Doce catchment in Brazil (71,492 grid nodes).}
\usage{data(riodoce.grids)}
\format{
Data frame contains the following columns:
  \describe{
  \item{\code{TAXBRC3}}{factor; dominant soil class from a 1:600,000 scale polgyon map, and based on the \href{http://www.cnps.embrapa.br/sibcs/}{Brazilian Soil Classification System}}
  \item{\code{GEOGLS3}}{factor; parent material class from the geological map}
  \item{\code{DEMSRT3}}{numeric; SRTM DEM}  
  \item{\code{SLPSRT3}}{numeric; slope in degrees based on the SRTM DEM}
  \item{\code{TWISRT2}}{numeric; SAGA Topographic Wetness Index based on the SRTM DEM}
  \item{\code{INSSRT3}}{numeric; potential incoming solar radiation derived in SAGA GIS}
  \item{\code{EV1MOD3}}{numeric; first principal component of the monthly MODIS EVI time series data}
  \item{\code{EV2MOD3}}{numeric; second principal component of the monthly MODIS EVI time series data}
  \item{\code{TD1MOD3}}{numeric; first principal component of the 8-day MODIS day-time LST time series data}         
  \item{\code{TN1MOD3}}{numeric; first principal component of the 8-day MODIS night-time LST time series data} 
  \item{\code{x}}{numeric; x-coordinate in the UTM 23S coordiante system with GRS67 ellipsoid}
  \item{\code{y}}{numeric; y-coordinate in the UTM 23S coordiante system with GRS67 ellipsoid}
}
}
\details{Soil classes ("AR", "CH", "CX", "CY", "LV", "LVA", "PV", "PVA", "RL"; Brazilian Soil Classification System) were obtained from the General soil map of Minas Gerais at scale 1:600,000. Land surface parameters: potential incoming solar radiation, SAGA Topographic Wetness Index and slope map were derived in SAGA GIS using the 1 km resolution SRTM DEM obtained from the Consultative Group on International Agricultural Research (CGIAR). Principal components for the Enhanced Vegetation Index (EVI, MOD13A3 product) and Land Surface Temperature (LST, MOD11A2 product) images were derived using the 8-day and monthly time-series MODIS images. }
\note{Listed layer follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{DEM} (digital elevation model); the next three letters represent the data source or collection method e.g. \code{SRT} (SRTM); the 6th character is the effective scale e.g. \code{3} indicates the 2nd standard scale i.e. 1/120 decimal degrees (about 1 km).}
\author{ Tomislav Hengl & Eliana de Souza}
\references{
\itemize{
\item EMBRAPA (2006) \href{http://www.cnps.embrapa.br/sibcs/}{Sistema Brasileiro de \enc{Classificação}{Classificacao} de Solos}. Centro Nacional de Pesquisa de Solos. 2nd Ed, Rio de Janeiro, Brasil.
\item WorldGrids repository of soil covariate layers (\url{http://worldgrids.org})
}
}
\seealso{ \code{\link{riodoce}}}
\examples{
\dontrun{# Load the data from GSIF page:
download.file("http://gsif.r-forge.r-project.org/riodoce.grids.rda", "riodoce.grids.rda")
load("riodoce.grids.rda")
gridded(riodoce.grids) <- ~x+y
utm.csy <- "+proj=utm +zone=23 +south +ellps=GRS67 +units=m +no_defs"
proj4string(riodoce.grids) <- CRS(utm.csy)
# plot maps:
col = grey(rev(seq(0.1,0.975,0.025)))
names(riodoce.grids)
spplot(riodoce.grids[9:10], col.regions=col)
# write to a GIS format:
writeGDAL(riodoce.grids["EV1MOD3"], "EV1MOD3.tif", "GTiff", mvFlag=-99999)
}
}
\keyword{datasets}