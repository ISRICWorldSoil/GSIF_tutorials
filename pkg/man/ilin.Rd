\name{ilin}
\docType{data}
\alias{ilin}
\title{NSCD soil samples for Illinois and Indiana States}
\description{A subset of the \href{http://soils.usda.gov/survey/nscd/}{National Cooperative Soil Characterization Database} (NSCD) covering the Illinois and Indiana States (5193 soil profiles). Analytical procedures and methods of soil preparation are identified by methods codes in the column headings of the standard NSSC data sheets and are described in \href{http://soils.usda.gov/survey/nscd/lim/}{Soil Survey Laboratory Information Manual No. 42}.}
\usage{data(ilin)}
\format{
Sites table. Contains the following columns:
  \describe{
  \item{\code{site_key}}{factor; unique name that identifies a particular site}
  \item{\code{user_site_id}}{character; a short label to help a user identify a particular site}
  \item{\code{LAT}}{numeric; latitude in decimal degrees on the NAD83 datum}
  \item{\code{LON}}{numeric; longitude in decimal degrees on the NAD83 datum}
  \item{\code{horizontal_datum_name}}{character; horizontal datum name (if different from NAD83)}
  \item{\code{observation_date}}{character; the date on which this particular soil was described or sampled}
  \item{\code{sampled_taxon_name}}{character; KST soil series or taxadjunct name e.g. "Paxton"}
  \item{\code{correlated_taxon_name}}{character; correlated name for soil series}
  \item{\code{correlated_taxon_kind}}{character; taxon type --- soil series or taxadjunct}
  \item{\code{correlated_class_name}}{character; correlated Keys to Soil Taxonomy taxon complete name e.g. "Sandy, mixed, thermic Ustic Torriorthent"}
  \item{\code{SSL_class_name}}{characted; original Keys to Soil Taxonomy taxon designation}
}
Horizons table. Contains the following columns:
  \describe{
  \item{\code{site_key}}{factor; unique name that identifies a particular site}
  \item{\code{natural_key}}{character; an integer number that uniquely identifies a particular horizon}
  \item{\code{layer_sequence}}{integer; a layer sequence number 1 to N}
  \item{\code{hzn_desgn}}{character; pedological designation of the horizon e.g. "Bt1"}
  \item{\code{hzn_bot}}{numeric; lower horizon depth from the surface in cm}
  \item{\code{hzn_top}}{numeric; upper horizon depth from the surface in cm}
  \item{\code{hzn_vert_subdvn}}{integer; an integer number that identifies horizon subdivision e.g. "Bt1" and "Bt2" both receive number "2"}
  \item{\code{clay_tot_psa}}{numeric; weight percentage of the clay particles (<0.0002 mm)}
  \item{\code{sand_tot_psa}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm)}
  \item{\code{silt_tot_psa}}{numeric; weight percentage of the sand particles (0.05--2 mm)}
  \item{\code{wpg2}}{numeric; weight percentage of coarse fragments (> 2 mm)}
  \item{\code{caco3}}{numeric; weight percentage of carbonates in the soil}
  \item{\code{oc}}{numeric; weight percentage of organic carbon}
  \item{\code{c_tot}}{numeric; weight percentage of the total organic and inorganic carbon}
  \item{\code{base_sum}}{numeric; base saturation in percent}
  \item{\code{cec_sum}}{numeric; Cation Exchange Capacity in cmol(+)/kg (summary)}
  \item{\code{cec_nh4}}{numeric; Cation Exchange Capacity in cmol(+)/kg (NH4OAc)}
  \item{\code{ph_h2o}}{numeric; pH index measured in water solution}
  \item{\code{ph_kcl}}{numeric; pH index measured in KCl solution}
  \item{\code{db_od}}{bulk density oven dry (105 C) in kg per cubic meter}
  \item{\code{layer_ID}}{factor; unique name that identifies a particular horizon}
}
}
\author{ The NSCD data is maintained by the NRCS NSSC Soil Survey Lab, Lincoln (NE). This subset has been prepared for the purpose of testing various DSM methods by Tomislav Hengl (tom.hengl@wur.nl). The processing steps are described in \href{http://gsif.r-forge.r-project.org/rw_NCSS.R}{this script}.}
\references{
\itemize{
\item Soil Data Subcommittee, (1997) \href{http://www.fgdc.gov/standards/projects/FGDC-standards-projects/soils/soil997.PDF}{Soil Geographic Data Standard}. Federal Geographic Data Committee, p. 286.
\item United States Department of Agriculture, Natural Resources Conservation Service, (2010) \href{http://soils.usda.gov/technical/classification/tax_keys/}{Keys to Soil Taxonomy}, Eleventh Edition. Soil Survey Staff.
\item National Cooperative Soil Characterization Database (\url{http://soils.usda.gov/survey/nscd/})
}
}
\note{ This data set is only a sample of the complete laboratory data. The original database contains over 30 thousand profiles. The variables listed above were selected for testing algorithms to produce soil property maps according to the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}. The total organic carbon for alkaline soils (pH in water > 7) needs to be corrected e.g. by removing 12 percent of weight percentage of CaCO3. \cr
Disclaimer: The assessment of the accuracy and applicability of the NSCD is strictly a USER RESPONSIBILITY. NRCS and NCSS take no responsibility for problems that arise from use of this database.
}
\seealso{ \code{\link{ilin.grids}}, \code{\link{geochm}}}
\examples{
\dontrun{
require(aqp)
data(ilin)
coordinates(ilin$sites) <- ~ LON + LAT
proj4string(ilin$sites) <- CRS("+proj=longlat +datum=NAD83")
# obtain state borders:
require(maps)
states.m = map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(states.m$names, ":"), function(x) x[1])
require(maptools)
states <- as(map2SpatialPolygons(states.m, IDs=IDs), "SpatialLines")
# overlay and plot points and maps:
plot(spTransform(ilin$sites, CRS("+proj=longlat +datum=WGS84")), pch="+", cex=.6, col="black")
lines(states)
}
}
\keyword{datasets}
