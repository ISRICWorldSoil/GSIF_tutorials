\name{geochm}
\docType{data}
\alias{geochm}
\title{NGS database samples for Indiana State}
\description{A subset of the National Geochemical Survey (NGS) samples covering the Indiana State. Contains a total of 879 point samples.}
\usage{data(geochm)}
\format{
Data frame; contains the following columns:
  \describe{
  \item{\code{REC_NO}}{factor; unique record identifier}
  \item{\code{DATASET}}{factor; abbreviated dataset group e.g. "AK+MI"}
  \item{\code{TYPEDESC}}{factor; abbreviated description of sample type: stream, pond, spring, soil etc}
  \item{\code{COLL_DATE}}{integer; sampling date}
  \item{\code{LONGITUDE}}{numeric; longitude in decimal degrees (NAD27 datum)}
  \item{\code{LATITIUDE}}{numeric; latitude in decimal degrees (NAD27 datum)}
  \item{\code{DATUM}}{factor; geodetic datum if different from NAD83}
  \item{\code{RELIEF}}{factor; relief in drainage basin from which sample was collected}
  \item{\code{FORMATION}}{factor; code or name of geologic formation in which sample area was located}
  \item{\code{ROCK_TYPE}}{factor; rock type in area of sample collection e.g. "carbonate"}
  \item{\code{SOIL_HORIZ}}{factor; soil horizon from which the sample was collected}
  \item{\code{COLOR}}{factor; observed color of powdered sample during splitting}
  \item{\code{MEDIUM}}{factor; sample medium --- rock, sediment, standard, or unknown}
  \item{\code{SOURCE}}{factor; geological source of the sample medium that was collected e.g. "Beach"}
  \item{\code{AS_ICP40}}{numeric; As (ppm) by Inductively Coupled Plasma Spectrometry (ICP) after acid dissolution}
  \item{\code{CD_ICP40}}{numeric; Cd (ppm)}
  \item{\code{CR_ICP40}}{numeric; Cr (ppm)}
  \item{\code{CU_ICP40}}{numeric; Cu (ppm)}
  \item{\code{NI_ICP40}}{numeric; Ni (ppm)}
  \item{\code{ZN_ICP40}}{numeric; Zn (ppm)}
  \item{\code{AS_AA}}{numeric; As (ppm) by Hydride Atomic Absorption}
  \item{\code{HG_AA}}{numeric; Hg (ppm) by Hydride Atomic Absorption}
  \item{\code{PB_ICP40}}{numeric; Pb (ppm)}
  \item{\code{C_TOT}}{numeric; total carbon (weight percentage) by combustion}
  \item{\code{C_ORG}}{numeric; organic carbon (weight percentage) as a difference between \code{C_TOT} and \code{C_CO3}}
  \item{\code{C_CO3}}{numeric; carbonate carbon (weight percentage) by Coulometric Titration}
  \item{\code{S_TOT}}{numeric; total sulfur (weight percentage) by combustion}
}
}
\author{ National Geochemical Survey database is maintaned by the USGS (contact: Peter Schweitzer). This subset has been prepared for the purpose of testing various geostatistical mapping algoriths by Tomislav Hengl (tom.hengl@wur.nl).}
\references{
\itemize{
\item The National Geochemical Survey Team, (2008) \href{http://tin.er.usgs.gov/geochem/doc/home.htm}{The National Geochemical Survey: database and documentation}. U.S. Geological Survey Open-File Report 2004-1001, U.S. Geological Survey, Reston VA.
\item National Geochemical Survey database (\url{http://tin.er.usgs.gov/geochem/})
}
}
\note{Negative values of the heavy metal concetrations indicate a determination that is below the limit of detection for the analytic method used. The magnitude of the negative number indicates the detection limit. For example, -10 ppm means the result should be regarded as < 10 ppm.
}
\examples{
# Load the NGS data:
data(geochm)
coordinates(geochm) <- ~LONGITUDE+LATITUDE
proj4string(geochm) <- CRS("+proj=longlat +ellps=clrk66 +datum=NAD27 +no_defs")
data(indiana)
border <- SpatialLines(list(Lines(list(Line(indiana)), "border")))
proj4string(border) <- aea.prj
par(mar=c(0,0,0,0))
plot(spTransform(border, geochm@proj4string))
points(geochm, pch="+", cex=.7)
}
\keyword{datasets}
