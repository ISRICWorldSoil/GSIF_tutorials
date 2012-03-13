\name{hoosier}
\docType{data}
\alias{hoosier}
\alias{indiana}
\title{Hoosier National Forest soil covariates at 1 km resolution}
\description{A list of standard soil mapping covariates prepared for the souther Indiana state (Hoosier National Forest area; 33,581 grid nodes). Projected in the Albers equal-area projection system that is commonly used to visualize the contiguous 48-state area.}
\usage{data(hoosier)}
\format{
Data frame contains the following columns:
  \describe{
  \item{\code{BRKUSG3}}{factor; USGS-derived bedrock type map}
  \item{\code{CRPUSG3}}{numeric; percentage (0--1) of crop fields based on the \code{LC6USGS}}
  \item{\code{DEMSRT3}}{numeric; SRTM DEM elevations in meters}
  \item{\code{DINSRT3}}{numeric; diffuse insolation derived in SAGA GIS}  
  \item{\code{EV1MOD3}}{numeric; MODIS time series (year 2011) Enhanced Vegetation Index PC1}
  \item{\code{EV2MOD3}}{numeric; MODIS time series (year 2011) Enhanced Vegetation Index PC2}
  \item{\code{GLCMER3}}{factor; GlobCov land cover values}
  \item{\code{INSSRT3}}{numeric; potential incoming solar radiation derived in SAGA GIS}
  \item{\code{LC6USG3}}{factor; USGS National Land Cover Database values (year 2006)}
  \item{\code{LNDUSG3}}{numeric; land mask for the Indiana state (0--1)}
  \item{\code{NL1DMS3}}{numeric; lights at night time series (1998--2010) images PC1}
  \item{\code{NL2DMS3}}{numeric; lights at night time series (1998--2010) images PC2}
  \item{\code{TD1MOD3}}{numeric; MODIS time series (year 2011) daytime Land Surface Temperature PC1}
  \item{\code{TD2MOD3}}{numeric; MODIS time series (year 2011) daytime Land Surface Temperature PC2}
  \item{\code{TN1MOD3}}{numeric; MODIS time series (year 2011) nighttime Land Surface Temperature PC1}
  \item{\code{TN2MOD3}}{numeric; MODIS time series (year 2011) nighttime Land Surface Temperature PC2}
  \item{\code{TRITOX3}}{numeric; kernel density of the Toxics Release Inventory (EPA) locations}
  \item{\code{TWISRT3}}{numeric; SAGA Topographic Wetness Index based on the SRTM DEM}
  \item{\code{SLPSRT3}}{numeric; slope in degrees based on the SRTM DEM}
  \item{\code{SMUSTS3}}{factor; mapping unit key from the General Soil Map STATSGO (at scale 1:250,000) obtained from the \href{http://datagateway.nrcs.usda.gov}{NRCS data gateway}}         
  \item{\code{SRGUSG3}}{factor; USGS-derived surface geology map}
  \item{\code{C01SUR3}}{numeric; percentage (0--1) of "Epiaqualfs"}
  \item{\code{C02SUR3}}{numeric; percentage (0--1) of "Endoaquolls"}
  \item{\code{C03SUR3}}{numeric; percentage (0--1) of "Hapludalfs"}
  \item{\code{C04SUR3}}{numeric; percentage (0--1) of "Eutrudepts"}
  \item{\code{C05SUR3}}{numeric; percentage (0--1) of "Haplosaprists"}
  \item{\code{C06SUR3}}{numeric; percentage (0--1) of "Hapludolls"}
  \item{\code{C07SUR3}}{numeric; percentage (0--1) of "Argiaquolls"}
  \item{\code{C09SUR3}}{numeric; percentage (0--1) of "Udorthents"}
  \item{\code{C10SUR3}}{numeric; percentage (0--1) of "Humaquepts"}
  \item{\code{C11SUR3}}{numeric; percentage (0--1) of "Fragiudalfs"}
  \item{\code{C12SUR3}}{numeric; percentage (0--1) of "Endoaqualfs"}
  \item{\code{C13SUR3}}{numeric; percentage (0--1) of "Fragiaqualfs"}
  \item{\code{C14SUR3}}{numeric; percentage (0--1) of "Dystrudepts"}
  \item{\code{C15SUR3}}{numeric; percentage (0--1) of "Fluvaquents"}
  \item{\code{C16SUR3}}{numeric; percentage (0--1) of "Hapludults"}
  \item{\code{C17SUR3}}{numeric; percentage (0--1) of "Argiudolls"}
  \item{\code{C18SUR3}}{numeric; percentage (0--1) of "Endoaquepts"}
  \item{\code{C19SUR3}}{numeric; percentage (0--1) of "Fragiudults"}
  \item{\code{C20SUR3}}{numeric; percentage (0--1) of "Glossaqualfs"}
  \item{\code{C21SUR3}}{numeric; percentage (0--1) of "Udipsamments"}
  \item{\code{C22SUR3}}{numeric; percentage (0--1) of "Paleudults"}
  \item{\code{C23SUR3}}{numeric; percentage (0--1) of "Paleudalfs"}
  \item{\code{C24SUR3}}{numeric; percentage (0--1) of "Udifluvents"}
  \item{\code{C25SUR3}}{numeric; percentage (0--1) of "Udarents"}
  \item{\code{C26SUR3}}{numeric; percentage (0--1) of "Epiaquepts"}
  \item{\code{C27SUR3}}{numeric; percentage (0--1) of "Psammaquents"}
  \item{\code{C31SUR3}}{numeric; percentage (0--1) of "Fragiaquults"}
  \item{\code{C33SUR3}}{numeric; percentage (0--1) of "Epiaquolls"}          
  \item{\code{s1}}{numeric; x-coordinate (m) in the Albers equal-area projection system}
  \item{\code{s2}}{numeric; y-coordinate (m) in the Albers equal-area projection system}
}
}
\author{ Tomislav Hengl and Zamir Libohova }
\details{ These predictors were selected following the CLORPT concept of factors of soil formation, first introduced by Hans Jenny (1941), where C stands for climate, L for landform, O for organisms, R for relief, P for parent material and T for time. The study area spans four 1--degree cells i.e. from 38--40 degrees latitude and from -87 to -85 degrees longitude.\cr
The percentage of crop fields and the land mask (0--100 percent) were derived using the area weighted averaging from the 30 m resolution land cover map (\code{LC6USG3}). The land mask (0--100 percent) shows percentage of water bodies, medium and high intensity build-up (urban) areas. To create a boolean mask map, consider using some treshold value (e.g. 0.6). SRTM DEM, MODIS and GlobCov are publicly available global data sets commonly used as \emph{standard global thematic layers} for Digital Soil Mapping. These layers are avialable for download from the WorldGrids.org repository. Soil taxonomic classes (\code{C**SUR3} (at scale 1:24,000), mapping units \code{SMUSTS3}), high resolution land cover data (\code{LC6USG3}), and geological maps (\code{SRGUSG3} and \code{BRKUSG3}) are available only for this study area (USA National layers; obtained from the \href{http://seamless.usgs.gov/}{USGS Seamless data warehouse}). Lights at night images have been obtained from the \href{http://ngdc.noaa.gov/dmsp/downloadV4composites.html}{National Geophysical Data Centre}. Density of the industrial facilities / pollutants (Toxics Release Inventory) were obtained from the \href{http://toxmap.nlm.nih.gov/toxmap/}{U.S. National Library of Medicine}. The first component from any time-series data usually shows long term mean value, while the higher order component show seasonal variations.}
\note{Each listed layer follows a standard naming convention used by Worlgrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. TD1 (daily surface temperature PC1); the next three letters represent the data source or collection method e.g. MOD (MODIS day time imagery); the 6th character is the effective scale e.g. 3 indicates the third standard scale (1/120 decimal degrees or about 1 km).\cr
The great groups from the Soil Survey Geographic database (SSURGO) might show in some areas, i.e. for some soil taxa, artificial jumps (border areas between two soil surveys). SSURGO mapsheets are typically harmonized only after the reports are completed and such artifacts are common for regional soil surveys that are run independently one from the other.}
\references{
\itemize{
\item Savtchenko, A. and Ouzounov, D. and Ahmad, S. and Acker, J. and Leptoukh, G. and Koziana, J. and Nickless, D., (2004) \href{http://dx.doi.org/10.1016/j.asr.2004.03.012}{Terra and Aqua MODIS products available from NASA GES DAAC}. Advances in Space Research, 34(4): 710-714. 
\item Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. Soil Survey Geographic (SSURGO) Database for Indiata State. Available online at \url{http://soildatamart.nrcs.usda.gov}.
\item Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. U.S. General Soil Map (STATSGO2). Available online at \url{http://soildatamart.nrcs.usda.gov}.
\item Xian, G., and Homer C., (2010) \href{http://dx.doi.org/10.1016/j.rse.2010.02.018}{Updating the 2001 National Land Cover Database impervious surface products to 2006 using Landsat imagery change detection methods}. Remote Sensing of Environment, 114, 1676-1686.
}
}
\seealso{ \code{\link{NASIS}}, \code{\link{NCSS}}, \code{\link{geochm}}}
\examples{
data(hoosier)
data(indiana)
border <- SpatialLines(list(Lines(list(Line(indiana)), "border")))
gridded(hoosier) <- ~s1+s2
\dontrun{
# overlay profiles and grids:
lst = c("C02SUR3","C04SUR3","C23SUR3")
col = grey(rev(seq(0.1,0.975,0.025)))
spplot(hoosier[lst], col.regions=col, sp.layout=list("sp.lines", border))
# write to a GIS format:
writeGDAL(hoosier["TD1MOD3"], "TD1MOD3.sdat", "SAGA")
}
}
\keyword{datasets}