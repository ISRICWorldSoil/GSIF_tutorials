\name{ilin}
\docType{data}
\alias{ilin}
\title{NSCD soil samples for Illinois and Indiana States}
\description{A subset of the \href{http://soils.usda.gov/survey/nscd/}{National Cooperative Soil Characterization Database} (NSCD) covering the Illinois and Indiana States (5189 soil profiles). Analytical procedures and methods of soil preparation are identified by methods codes in the column headings of the standard NSSC data sheets and are described in \href{http://soils.usda.gov/survey/nscd/lim/}{Soil Survey Laboratory Information Manual No. 42}.}
\usage{data(ilin)}
\format{
Sites table. Contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique label to help a user identify a particular site (\code{user_site_id} in the NSCD)}
  \item{\code{LATNAD83}}{numeric; latitude in decimal degrees on the NAD83 datum}
  \item{\code{LONNAD83}}{numeric; longitude in decimal degrees on the NAD83 datum}
  \item{\code{TIMESTRR}}{Date; the date on which this particular soil was described or sampled (\code{observation_date} in the NSCD)}
  \item{\code{TAXSUSDA}}{factor; correlated name for soil series (\code{correlated_taxon_name} in the NSCD)}
  \item{\code{TAXKUSDA}}{factor; taxon type --- soil series or taxadjunct (\code{correlated_taxon_kind} in the NSCD)}
  \item{\code{TAXNUSDA}}{factor; correlated Keys to Soil Taxonomy taxon complete name e.g. \code{"Sandy, mixed, thermic Ustic Torriorthent"} (\code{correlated_class_name} in the NSCD)}
}
Horizons table. Contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; a short label to help a user identify a particular site (\code{user_site_id} in the NSCD)}
  \item{\code{LSQINT}}{integer; a layer sequence number 1 to N (\code{layer_sequence} in the NSCD)}
  \item{\code{HZDUSD}}{character; pedological designation of the horizon e.g. \code{"Bt1"} (\code{hzn_desgn} in the NSCD)}
  \item{\code{UHDICM}}{numeric; lower horizon depth from the surface in cm (\code{hzn_bot} in the NSCD)}
  \item{\code{LHDICM}}{numeric; upper horizon depth from the surface in cm (\code{hzn_top} in the NSCD)}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm; \code{clay_tot_psa} in the NSCD)}
  \item{\code{SNDPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm; \code{sand_tot_psa} in the NSCD)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm; \code{silt_tot_psa} in the NSCD)}
  \item{\code{CRFMAS}}{numeric; weight percentage of coarse fragments (> 2 mm; \code{wpg2} in the NSCD)}
  \item{\code{PHIHO5}}{numeric; pH index measured in water solution(\code{ph_h2o} in the NSCD)}
  \item{\code{PHIKCL}}{numeric; pH index measured in KCl solution (\code{ph_kcl} in the NSCD)}
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permille (see below)}
  \item{\code{BLDODR}}{bulk density oven dry (105 C) in kg per cubic meter (\code{db_od} in the NSCD)}
}
}
\author{ The NSCD data is maintained by the NRCS NSSC Soil Survey Lab, Lincoln (NE). This subset has been prepared for the purpose of testing various DSM methods by Tomislav Hengl <tom.hengl@wur.nl>. The processing steps are described in \href{http://gsif.r-forge.r-project.org/rw_NCSS.R}{this script}.}
\references{
\itemize{
\item Soil Data Subcommittee, (1997) \href{http://www.fgdc.gov/standards/projects/FGDC-standards-projects/soils/soil997.PDF}{Soil Geographic Data Standard}. Federal Geographic Data Committee, p. 286.
\item United States Department of Agriculture, Natural Resources Conservation Service, (2010) \href{http://soils.usda.gov/technical/classification/tax_keys/}{Keys to Soil Taxonomy}, Eleventh Edition. Soil Survey Staff.
\item National Cooperative Soil Characterization Database (\url{http://soils.usda.gov/survey/nscd/})
}
}
\note{ This data set is only a sample of the complete laboratory data. The original database contains over 30 thousand profiles. The variables listed above were selected for testing algorithms to produce soil property maps according to the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}. \code{ORCDRC} values are based on the \code{c_tot} and \code{oc} values in the NSCD. The total organic carbon in permille for alkaline soils (pH in water > 7) has been corrected by removing 12 percent of weight percentage of CaCO3. \cr
Disclaimer: The assessment of the accuracy and applicability of the NSCD is strictly a USER RESPONSIBILITY. NRCS and NCSS take no responsibility for problems that arise from use of this database.
}
\seealso{ \code{\link{ilin.grids}}, \code{\link{geochm}}}
\examples{
data(ilin)
coordinates(ilin$sites) <- ~ LONNAD83 + LATNAD83
proj4string(ilin$sites) <- CRS("+proj=longlat +datum=NAD83")
\dontrun{# obtain state borders:
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
