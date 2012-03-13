\name{NCSS}
\docType{data}
\alias{NCSS}
\title{NCSS soil samples for Indiana State}
\description{A subset of the \href{http://ssldata.nrcs.usda.gov/}{National Cooperative Soil Characterization Database} (NCSS) covering the Indiana State (2236 soil profiles). Analytical procedures and methods of soil preparation are identified by methods codes in the column headings of the standard NSSC data sheets and are described in \href{http://soils.usda.gov/survey/nscd/lim/}{Soil Survey Laboratory Information Manual No. 42}.}
\usage{data(NCSS)}
\format{
Sites table. Contains the following columns:
  \describe{
  \item{\code{site_key}}{factor; unique name that identifies a particular site}
  \item{\code{user_site_id}}{character; a short label to help a user identify a particular site}
  \item{\code{LAT}}{numeric; latitude in decimal degrees on the NAD83 datum}
  \item{\code{LON}}{numeric; longitude in decimal degrees on the NAD83 datum}
  \item{\code{horizontal_datum_name}}{character; horizontal datum name (if different from NAD83)}
  \item{\code{observation_date}}{character; the date on which this particular soil was described or sampled, expressed as month/day/year}
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
  \item{\code{ORCDRC}}{numeric; kg per cubic meter of organic carbon (corrected for CaCO3)}
}
}
\author{ The NCSS data is maintained by the NRCS NSSC Soil Survey Lab, Lincoln (NE). This subset has been prepared for the purpose of testing various DSM methods by Tomislav Hengl (tom.hengl@wur.nl). The processing steps are described in \href{http://gsif.r-forge.r-project.org/rw_NCSS.R}{this script}.}
\references{
\itemize{
\item Soil Data Subcommittee, (1997) \href{http://www.fgdc.gov/standards/projects/FGDC-standards-projects/soils/soil997.PDF}{Soil Geographic Data Standard}. Federal Geographic Data Committee, p. 286.
\item United States Department of Agriculture, Natural Resources Conservation Service, (2010) \href{http://soils.usda.gov/technical/classification/tax_keys/}{Keys to Soil Taxonomy}, Eleventh Edition. Soil Survey Staff.
\item National Cooperative Soil Characterization Database (\url{http://ssldata.nrcs.usda.gov/})
}
}
\note{ This data set is only a sample of the complete laboratory data. The original database contains over 30 thousand profiles. The variables listed above were selected for testing algorithms to produce soil property maps according to the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}.\cr
\code{ORCDRC} (in kg per cubic meter) has been estimated using bulk dens and gravel content measurements for each horizon. The total organic carbon is corrected (for pH in water > 7) by removing 12 percent of weight percentage of CaCO3. \cr
Disclaimer: The assessment of the accuracy and applicability of the NCSS is strictly a USER RESPONSIBILITY and the NRCS and NCSS take no responsibility for problems that arise from use of this database.
}
\seealso{ \code{\link{NASIS}}, \code{\link{hoosier}}, \code{\link{geochm}}}
\examples{
\dontrun{
require(aqp)
aea.prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
data(NCSS)
NCSS.spc <- join(NCSS$horizons, NCSS$sites[,c(1,3,4,7,8)], type='inner')
depths(NCSS.spc) <- site_key ~ hzn_top + hzn_bot
site(NCSS.spc) <- ~ LON + LAT + sampled_taxon_name + correlated_taxon_name 
# generate SpatialPoints:
coordinates(NCSS.spc) <- ~ LON + LAT
proj4string(NCSS.spc@sp) <- CRS("+proj=longlat +datum=NAD83")
# soil covariate layers:
data(hoosier)
gridded(hoosier) <- ~s1+s2
data(indiana)
border <- SpatialLines(list(Lines(list(Line(indiana)), "border")))
# overlay and plot points and maps:
pnts <- list("sp.points", spTransform(NCSS.spc@sp, CRS(aea.prj)), pch="+", col="black")
lns <- list("sp.lines", border, col="black")
plt1 <- spplot(hoosier["DEMSRT3"], col.regions=grey(rev(seq(0.25,0.975,0.025))), sp.layout = list(pnts, lns))
plt2 <- spplot(hoosier["C04SUR3"], col.regions=grey(rev(seq(0,0.975,0.025))), sp.layout = list(lns))
BRK <- hoosier["BRKUSG3"]
BRK$BRKUSG3 <- as.integer(BRK$BRKUSG3)
plt3 <- spplot(BRK, col.regions=grey(runif(length(levels(hoosier$BRKUSG3)))), sp.layout = list(lns))
print(plt1, split=c(1,1,3,1), more=TRUE)
print(plt2, split=c(2,1,3,1), more=TRUE)
print(plt3, split=c(3,1,3,1))
}
}
\keyword{datasets}
