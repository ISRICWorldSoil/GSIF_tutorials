\name{WPS-class}
\docType{class}
\alias{WPS-class}
\alias{over,WPS,SpatialPoints-method}
\alias{describe,WPS-method}
\alias{getProcess,WPS-method}
\alias{getProcess}
\alias{describe}
\alias{show,WPS-method}
\alias{subset,WPS-method}
\title{A class for a Web Processing Service}
\description{A class for a Web Processing Service. Can be used to overlay points or fetch grid values for rasters located remotely on a server and specified via the \code{inRastername} slot.}
\section{Slots}{
  \describe{
    \item{\code{server}:}{object of class \code{"list"}; contains the location of the CGI script that executes WPS (\code{"URI"}); service name (\code{"service.name"}), version (\code{"version"}), request type (\code{"request"}), identifier (\code{"identifier"})}
    \item{\code{inRastername}:}{object of class \code{"character"}; name of the objects on the server}
  }
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(object = "WPS")}: gets the complete server capabilities }
	\item{getProcess}{\code{signature(x = "WPS")}: gets a list of processes available from a server }
	\item{describe}{\code{signature(x = "WPS")}: lists parameters specific to some service identifier}
  \item{over}{\code{signature(x = "WPS", y = "SpatialPoints")}: overlays spatial points and the target grids defined via the WPS-class (point-by-point)}
  \item{subset}{\code{signature(x = "WPS")}: subsets a grid (from server) and loads it to R; use \code{bbox} argument to specify the bounding box}
  }
}
\note{More examples of overlay, subset and aggregation functions are available via \href{http://worldgrids.org/doku.php?id=wiki:functions.r}{WorldGrids.org}. WPS WorldGrids.org uses the \href{http://pywps.wald.intevation.org}{PyWPS} module on a Debian system with Webserver, GDAL, Python and Scipy. The standard format for the gridded data on the WorldGrids.org repository is \code{"GeoTiff"}. Use of the \code{"bbox"} object to obtain grids that cover more than 30 percent of the global coverage is not recommended. Consider instead downloading the compressed images directly from \href{http://worldgrids.org}{WorldGrids.org}.}
\references{
\itemize{
\item PyWPS module (\url{http://pywps.wald.intevation.org})
\item WorldGrids.org (\url{http://worldgrids.org})
}
}
\author{ Tomislav Hengl & Hannes I. Reuter}
\seealso{ \code{\link{landmask}}
}
\examples{
\dontrun{
library(XML)
library(sp)
URI = "http://wps.worldgrids.org/pywps.cgi"
server <- list(URI=URI, request="execute", 
    version="version=1.0.0", service.name="service=wps", 
    identifier="identifier=sampler_local1pt_nogml")
glcesa3.wps <- new("WPS", server=server, inRastername="glcesa3a")
# show(biocl15.wps)
prl <- getProcess(glcesa3.wps)
prl[7]
describe(glcesa3.wps, identifier="overlay")
p1 <- data.frame(lon=15, lat=15)
coordinates(p1) <- ~lon+lat
proj4string(p1) <- CRS("+proj=longlat +datum=WGS84")
p1
over(glcesa3.wps, p1)
# fetch grids and load the to R:
glcesa3 <- subset(glcesa3.wps, bbox=matrix(c(20,40,22,42), nrow=2))
image(glcesa3)
}
}
\keyword{classes}
