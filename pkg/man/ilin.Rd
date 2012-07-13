\name{ilin}
\docType{data}
\alias{ilin}
\alias{ilin.grids}
\title{NSCD soil samples for Illinois and Indiana States}
\description{A subset of the \href{http://soils.usda.gov/survey/nscd/}{National Cooperative Soil Characterization Database} (NSCD) covering the Illinois and Indiana States (5189 soil profiles). Analytical procedures and methods of soil preparation are identified by methods codes in the column headings of the standard NSSC data sheets and are described in \href{http://soils.usda.gov/survey/nscd/lim/}{Soil Survey Laboratory Information Manual No. 42}.}
\usage{data(ilin)}
\format{
The \code{ilin} data set contains two data frames --- sites and horizons. Sites table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique label to help a user identify a particular site (\code{user_site_id} in the NSCD)}
  \item{\code{LATNAD83}}{numeric; latitude in decimal degrees on the NAD83 datum}
  \item{\code{LONNAD83}}{numeric; longitude in decimal degrees on the NAD83 datum}
  \item{\code{TIMESTRR}}{Date; the date on which this particular soil was described or sampled (\code{observation_date} in the NSCD)}
  \item{\code{TAXSUSDA}}{factor; correlated name for soil series (\code{correlated_taxon_name} in the NSCD)}
  \item{\code{TAXKUSDA}}{factor; taxon type --- soil series or taxadjunct (\code{correlated_taxon_kind} in the NSCD)}
  \item{\code{TAXNUSDA}}{factor; correlated Keys to Soil Taxonomy taxon complete name e.g. \code{"Sandy, mixed, thermic Ustic Torriorthent"} (\code{correlated_class_name} in the NSCD)}
}
Horizons table contains the following columns:
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
The \code{ilin.grids} data frame contains a list of covariates prepared for the Illinois and Indiana States (40,383 grid nodes):
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
\details{The NSCD data set is only a sample of the complete laboratory data. The original database contains over 30 thousand profiles. The variables listed above were selected for testing algorithms to produce soil property maps according to the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}. \code{ORCDRC} values are based on the \code{c_tot} and \code{oc} values in the NSCD. The total organic carbon in permille for alkaline soils (pH in water > 7) has been corrected by removing 12 percent of weight percentage of CaCO3. Disclaimer: The assessment of the accuracy and applicability of the NSCD is strictly a USER RESPONSIBILITY. NRCS and NCSS take no responsibility for problems that arise from use of this database. 
}
\note{Majority of the gridded layers listed were obtained from the \href{http://seamless.usgs.gov/}{USGS Seamless data warehouse}. The percentage of crop fields and the land mask (0--100 percent) were derived using the area weighted averaging from the 30 m resolution National Land cover map from 2006. The land mask (0--100 percent) shows percentage of water bodies, medium and high intensity build-up (urban) areas. To create a boolean mask map, consider using some treshold value (e.g. 0.6). Density of the industrial facilities / pollutants (Toxics Release Inventory) were obtained from the \href{http://toxmap.nlm.nih.gov/toxmap/}{U.S. National Library of Medicine}. Land surface parameters diffuse insolation, potential incoming solar radiation, SAGA Topographic Wetness Index and slope map were derived in SAGA GIS using the 250 m resolution SRTM DEM obtained from the Consultative Group on International Agricultural Research (CGIAR), then resampled to 3/120 arcdegrees (2.5 km resolution). Soil suborder layer was derived as dominant taxonomic components (highest percentage) in the General Soil Map STATSGO2 (at scale 1:250,000). The processing steps used to prepare the map of dominatant taxonomic compenents are described in \href{http://gsif.r-forge.r-project.org/gsmsoilmu_a_us.R}{this script}. There are cases where the dominant suborder is less than 30\% of the map unit, which means that selecting a single suborer is not likely all that informative. To obtain additional covariate layers, e.g. SRTM DEM, MODIS, GlobCov and similar publicly available global data sets refer to the WorldGrids.org repository. Listed gridded layers follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{CRP} (percentage (0--1) of crop fields); the next three letters represent the data source or collection method e.g. \code{USG} (USGS); the 6th character is the effective scale e.g. \code{2} indicates the 3rd standard scale i.e. 3/120 decimal degrees (about 2.5 km).}
\author{ The NSCD data is maintained by the NRCS NSSC Soil Survey Lab, Lincoln (NE). This subset has been prepared for the purpose of testing various DSM methods by Tomislav Hengl <tom.hengl@wur.nl>. The processing steps are described in \href{http://gsif.r-forge.r-project.org/rw_NCSS.R}{this script}.}
\references{
\itemize{
\item Soil Data Subcommittee, (1997) \href{http://www.fgdc.gov/standards/projects/FGDC-standards-projects/soils/soil997.PDF}{Soil Geographic Data Standard}. Federal Geographic Data Committee, p. 286.
\item Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. U.S. General Soil Map (STATSGO2). Available online at \url{http://soildatamart.nrcs.usda.gov}.
\item Xian, G., and Homer C., (2010) \href{http://dx.doi.org/10.1016/j.rse.2010.02.018}{Updating the 2001 National Land Cover Database impervious surface products to 2006 using Landsat imagery change detection methods}. Remote Sensing of Environment, 114, 1676-1686.
\item United States Department of Agriculture, Natural Resources Conservation Service, (2010) \href{http://soils.usda.gov/technical/classification/tax_keys/}{Keys to Soil Taxonomy}, Eleventh Edition. Soil Survey Staff.
\item National Cooperative Soil Characterization Database (\url{http://soils.usda.gov/survey/nscd/})
\item soilDB package (\url{http://cran.r-project.org/package=soilDB})
}
}
\seealso{ \code{\link{geochm}}}
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
# download the gridded data:
download.file("http://worldgrids.org/rda/ilin.grids.rda", "ilin.grids.rda")
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
