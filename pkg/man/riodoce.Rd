\name{riodoce}
\docType{data}
\encoding{latin1}
\alias{riodoce}
\alias{riodoce.grids}
\title{Soil samples for the Rio Doce basin}
\description{Soil samples (456) collected for the Rio Doce basin in Brazil.}
\usage{data(riodoce)}
\format{
The \code{riodoce} data set contains two data frames --- sites and horizons. Sites table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique name that identifies a particular site}
  \item{\code{LATWGS84}}{numeric; latitude in decimal degrees on the NAD83 datum}
  \item{\code{LONWGS84}}{numeric; longitude in decimal degrees on the NAD83 datum}
  \item{\code{TAXBRC}}{factor; original soil class based on the \href{http://www.cnps.embrapa.br/sibcs/}{Brazilian Soil Classification System}}
}
Horizons table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique name that identifies a particular site}
  \item{\code{HRDLOC}}{character; local pedological designation of the horizon e.g. "Bt1"}
  \item{\code{UHDICM}}{numeric; upper horizon depth from the surface in cm}
  \item{\code{LHDICM}}{numeric; lower horizon depth from the surface in cm}
  \item{\code{PHIHO5}}{numeric; pH index measured in water solution}
  \item{\code{PHIKCL}}{numeric; pH index measured in KCl solution}
  \item{\code{ORCDRC}}{numeric; weight percentage in permilles of the total organic carbon}
  \item{\code{SNDPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm)}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm)}
  \item{\code{BSABCL}}{numeric; base saturation in percent}
  \item{\code{CEXBCL}}{numeric; Cation Exchange Capacity in cmol(+)/kg (summary)}
  \item{\code{BLDVOL}}{bulk density oven dry (105 C) in kg per cubic meter}
  \item{\code{AWAIMM}}{total available water capacity in mm}
  \item{\code{HORUID}}{factor; unique name that identifies a particular horizon}
}
The \code{riodoce.grids} data frame contains a list of soil covariates prepared for the Rio Doce catchment in Brazil (71,492 grid nodes):
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
\note{Soil classes in \code{TAXBRC3} ("AR", "CH", "CX", "CY", "LV", "LVA", "PV", "PVA", "RL"; Brazilian Soil Classification System) were obtained from the General soil map of Minas Gerais at scale 1:600,000. Land surface parameters: potential incoming solar radiation, SAGA Topographic Wetness Index and slope map were derived in SAGA GIS using the 1 km resolution SRTM DEM obtained from the Consultative Group on International Agricultural Research (CGIAR). Principal components for the Enhanced Vegetation Index (EVI, MOD13A3 product) and Land Surface Temperature (LST, MOD11A2 product) images were derived using the 8-day and monthly time-series MODIS images. Listed gridded layers follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{DEM} (digital elevation model); the next three letters represent the data source or collection method e.g. \code{SRT} (SRTM); the 6th character is the effective scale e.g. \code{3} indicates the 2nd standard scale i.e. 1/120 decimal degrees (about 1 km).}
\author{This data set has been prepared for the purpose of testing various soil mapping methods by Eliana de Souza and Tomislav Hengl (see \href{http://gsif.r-forge.r-project.org/rw_riodoce.R}{processing steps} for more info).}
\references{
\itemize{
\item Benedetti, M.M. and Curi, N. and Sparovek, G. and de Carvalho Filho, A. and Silva, S.H.G., (2011) Updated Brazilian's Georeferenced Soil Database –-- an Improvement for International Scientific Information Exchanging. In: Principles, Application and Assessment in Soil Science (Edited by E. Burcu \enc{Özkaraova}{Ozkaraova} \enc{Güngör}{Gungor}), InTech, chapter 16, pp. 307-330.
\item EMBRAPA (2006) \href{http://www.cnps.embrapa.br/sibcs/}{Sistema Brasileiro de \enc{Classificação}{Classificacao} de Solos}. Centro Nacional de Pesquisa de Solos. 2nd Ed, Rio de Janeiro, Brasil.
\item Filho, E.F. et al. (2011) Levantamento de solos e \enc{aptidão}{aptidao} \enc{agrícola}{agricola} das terras da bacia do Rio Doce, Estado de Minas Gerais. Belo Horizonte/MG, 453 p. \url{www.meioambiente.mg.gov.br}
\item Cooper, M., Mendes. L.M.S., Silva, W.L.C, and Sparovek, G., (2005) \href{https://www.soils.org/publications/sssaj/articles/69/3/0649}{A National Soil Profile Database for Brazil Available to International Scientists}. Soil Science Society of America Journal, 69:649-652.
}
}
\examples{
\dontrun{
library(aqp)
data(riodoce)
# Promote to SoilProfileCollection
riodoce <- join(riodoce$horizons, riodoce$sites, type='inner')
# specify point IDs and depths:
depths(riodoce) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(riodoce) <- ~ LONWGS84 + LATWGS84 + TAXBRC  
# generate SpatialPoints
coordinates(riodoce) <- ~ LONWGS84 + LATWGS84
proj4string(riodoce) <- "+proj=latlong +datum=WGS84"
str(riodoce)
# load the gridded data:
download.file("http://worldgrids.org/rda/riodoce.grids.rda", "riodoce.grids.rda")
load("riodoce.grids.rda")
gridded(riodoce.grids) <- ~x+y
utm.csy <- "+proj=utm +zone=23 +south +ellps=GRS67 +units=m +no_defs"
proj4string(riodoce.grids) <- CRS(utm.csy)
# plot maps:
names(riodoce.grids)
image(riodoce.grids["DEMSRT3"])
points(spTransform(riodoce@sp, CRS(utm.csy)), pch="+")
# write to a GIS format:
writeGDAL(riodoce.grids["EV1MOD3"], "EV1MOD3.tif", "GTiff", mvFlag=-99999)
}
}
\keyword{datasets}
