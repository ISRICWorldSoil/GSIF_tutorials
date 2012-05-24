\name{riodoce}
\docType{data}
\encoding{latin1}
\alias{riodoce}
\title{Soil samples for the Rio Doce basin}
\description{Soil samples (456) collected for the Rio Doce basin in Brazil.}
\usage{data(riodoce)}
\format{
Sites table. Contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique name that identifies a particular site}
  \item{\code{LAT}}{numeric; latitude in decimal degrees on the NAD83 datum}
  \item{\code{LON}}{numeric; longitude in decimal degrees on the NAD83 datum}
  \item{\code{TAXBRC}}{factor; original soil class based on the \href{http://www.cnps.embrapa.br/sibcs/}{Brazilian Soil Classification System}}
}
Horizons table. Contains the following columns:
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
}
\author{This soil profile data set has been prepared for the purpose of testing various soil mapping methods by Eliana de Souza and Tomislav Hengl. The processing steps are described in \href{http://gsif.r-forge.r-project.org/rw_riodoce.R}{this script}.}
\references{
\itemize{
\item Benedetti, M.M. and Curi, N. and Sparovek, G. and de Carvalho Filho, A. and Silva, S.H.G., (2011) Updated Brazilian's Georeferenced Soil Database – an Improvement for International Scientific Information Exchanging. In: Principles, Application and Assessment in Soil Science (Edited by E. Burcu \enc{Ozkaraova}{Özkaraova} \enc{Gungor}{Güngör}), InTech, chapter 16, pp. 307-330.
\item Filho, E.F. et al. (2011) Levantamento de solos e \enc{aptidao}{aptidão} \enc{agricola}{agrícola} das terras da bacia do Rio Doce, Estado de Minas Gerais. Belo Horizonte/MG, 453 p. \url{www.meioambiente.mg.gov.br}
\item Cooper, M., Mendes. L.M.S., Silva, W.L.C, and Sparovek, G., (2005) \href{https://www.soils.org/publications/sssaj/articles/69/3/0649}{A National Soil Profile Database for Brazil Available to International Scientists}. Soil Science Society of America Journal, 69:649-652.
}
}
\seealso{ \code{\link{riodoce.grids}}}
\examples{
\dontrun{
require(aqp)
data(riodoce)
depths(riodoce) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(riodoce) <- ~ LONWGS84 + LATWGS84 + TAXBRC  
# generate SpatialPoints
coordinates(riodoce) <- ~ LONWGS84 + LATWGS84
proj4string(riodoce) <- "+proj=latlong +datum=WGS84"
str(riodoce)
}
}
\keyword{datasets}
