\name{REST.SoilGrids-class}
\docType{class}
\alias{REST.SoilGrids-class}
\alias{REST.SoilGrids}
\alias{over,REST.SoilGrids,SpatialPoints-method}
\title{A class for SoilGrids REST API}
\description{A class for \href{http://rest.soilgrids.org/}{SoilGrids REST API} Service. Can be used to overlay points or fetch grid values from SoilGrids Soil Information System.}
\section{Slots}{
  \describe{
    \item{\code{server}:}{object of class \code{"character"}; contains the location of the server that executes REST.SoilGrids calls}
    \item{\code{query}:}{object of class \code{"list"}; contains parameters or REST.SoilGrids query}
    \item{\code{stream}:}{object of class \code{"character"}; contains parameters or REST.SoilGrids stream operation}
  }
}
\section{Methods}{
  \describe{
  \item{over}{\code{signature(x = "REST.SoilGrids", y = "SpatialPoints")}: overlays spatial points and the target grids defined via the REST.SoilGrids-class (point-by-point) and returns list of objects of \code{"SpatialPixelsDataFrame"}-class}
  }
}
\note{More examples of overlay and download functions are available via \url{http://rest.soilgrids.org/}. \code{over} method is not recommended for large point data sets.}
\references{
\itemize{
\item SoilGrids --- a system for automated soil mapping (\url{http://www.soilgrids.org})
}
}
\author{ Tomislav Hengl & Jorge S. Mendes de Jesus}
\seealso{ \code{\link{SoilGrids-class}}, \code{\link{WPS-class}}
}
\examples{
\dontrun{
library(rjson)
library(sp)
## 2 points:
pnts <- data.frame(lon=c(10.65,5.36), lat=c(51.81,51.48), id=c("p1","p2"))
coordinates(pnts) <- ~lon+lat
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
pnts
## REST example:
soilgrids.r <- REST.SoilGrids(c("ORCDRC","PHIHOX"))
ov <- over(soilgrids.r, pnts)
ORCDRC.pnt1 <- data.frame(
  top=unlist(ov[1,grep("depthCodesMeters", names(ov))])*-100, 
  M=unlist(ov[1,grep("ORCDRC.M", names(ov))]), 
  L=unlist(ov[1,grep("ORCDRC.L", names(ov))]), 
  U=unlist(ov[1,grep("ORCDRC.U", names(ov))]))
ORCDRC.pnt1$variable <- "ORCDRC"
## plot the result:
library(lattice)
library(aqp)
data(soil.legends)
## Soil organic carbon:
ORCDRC.range = range(soil.legends[["ORCDRC"]]$MIN, soil.legends[["ORCDRC"]]$MAX)
dev.new(width=5, height=6)
xyplot(top ~ M | variable, data=ORCDRC.pnt1, ylab='Depth in cm',
  xlab='5th and 95th percentiles', xlim=ORCDRC.range,
  lower=ORCDRC.pnt1$L, upper=ORCDRC.pnt1$U, ylim=c(150,0),
  panel=panel.depth_function,
  alpha=0.25, sync.colors=TRUE,
  par.settings=list(superpose.line=list(col='RoyalBlue', lwd=3)),
  strip=strip.custom(bg=grey(0.8))
)

## Soil pH:
PHIHOX.range = range(soil.legends[["PHIHOX"]]$MIN, soil.legends[["PHIHOX"]]$MAX)
for(i in 1:nrow(ov)){
  PHIHOX.pnt <- data.frame(
    top=unlist(ov[i,grep("depthCodesMeters", names(ov))])*-100, 
    M=unlist(ov[i,grep("PHIHOX.M", names(ov))]), 
    L=unlist(ov[i,grep("PHIHOX.L", names(ov))]), 
    U=unlist(ov[i,grep("PHIHOX.U", names(ov))]))
  PHIHOX.pnt$variable <- "PHIHOX"
  png(paste("PHIHOX_depth_", i, ".png", sep=""), width=300, height=6/5*300)
  p <- xyplot(top ~ M/10 | variable, data=PHIHOX.pnt, ylab='Depth in cm',
    xlab='5th and 95th percentiles', xlim=PHIHOX.range/10,
    lower=PHIHOX.pnt$L/10, upper=PHIHOX.pnt$U/10, ylim=c(150,0),
    panel=panel.depth_function,
    alpha=0.25, sync.colors=TRUE,
    par.settings=list(superpose.line=list(col='Red', lwd=3)),
    strip=strip.custom(bg=grey(0.8))
  )
  print(p)
  graphics.off()
}
## plot in Google Earth:
library(plotKML)
kml(pnts, colour=id, file="PHIHOX_depth.kml", 
  shape=paste("PHIHOX_depth_", 1:nrow(ov), ".png", sep=""), 
  size=6, points_names=pnts$id, 
  colour_scale=rep("#FFFFFF", 2))
}
}
\keyword{classes}
