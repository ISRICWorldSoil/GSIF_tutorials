\name{GSIF.env}
\alias{GSIF.env}
\alias{GSIF.opts}
\title{GSIF specific environmental variables / paths}
\description{Sets the environmental, package specific parameters and settings (URLs, names, default cell size and similar) that can be later on passed to other functions.}
\usage{GSIF.env(wps.server = "http://wps.worldgrids.org", 
    ref_CRS = "+proj=longlat +datum=WGS84",
    NAflag = -99999, 
    license_url = "http://creativecommons.org/licenses/by/3.0/", 
    project_url = "http://gsif.r-forge.r-project.org/", 
    stdepths = c(-2.5, -10, -22.5, -45, -80, -150)/100, 
    stsize = c(5, 10, 15, 30, 40, 100)/100, 
    cellsize = rev(c(6/120, 3/120, 1/120, 1/240, 1/600, 1/1200, 1/3600)),
    REST.server = 'http://rest.soilgrids.org/',
    attributes = c("ORCDRC","PHIHOX","SNDPPT","SLTPPT","CLYPPT",
                   "CFRVOL","CEC","BLD","TAXGWRB","TAXOUSDA"),
    TimeSpan = list(begin=as.POSIXct("1950-01-01"), end=as.POSIXct("2005-12-30")),
    show.env = TRUE)
}
\arguments{
  \item{wps.server}{character; location of the WPS server}
  \item{ref_CRS}{the referent CRS proj4string (\code{"+proj=longlat +datum=WGS84"})}
  \item{NAflag}{the default missing value flag (usually \code{"-99999"})}
  \item{license_url}{the default license URL}
  \item{project_url}{the default location of the package documentation}
  \item{stdepths}{numeric; standard depths}
  \item{stsize}{numeric; standard horizon thicknesses}
  \item{cellsize}{numeric; standard grid cell sizes on WGS84 geographical coordinates}
  \item{REST.server}{character; location of the SoilGrids REST service}
  \item{attributes }{character; default soil variables of interest}
  \item{TimeSpan}{list; default begin end times (temporal coverage of SoilGrids)}
  \item{show.env}{logical; specify whether to print all environmental parameters}
}
\author{Tomislav Hengl}
\note{To further customize the GSIF options, consider putting:\cr
 
\code{library(GSIF); GSIF.env(\dots, show.env = FALSE)}\cr 

in your \code{"/etc/Rprofile.site"}.
}
\examples{
# environmental variables:
GSIF.env()
get("cellsize", envir = GSIF.opts)
}
\keyword{options}