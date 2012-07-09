\name{GlobalSoilMap-class}
\docType{class}
\alias{GlobalSoilMap-class}
\title{A class for GlobalSoilMap soil property maps}
\description{A class containing multiple realizations of soil properties following the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}, i.e. predictions of some of the target global soil properties at six depths (\code{sd1} = 2.5 cm, \code{sd2} = 10 cm, \code{sd3} = 22.5 cm, \code{sd4} = 45 cm, \code{sd5} = 80 cm, \code{sd6} = 150 cm).}
\section{Slots}{
  \describe{    
    \item{\code{varname}:}{object of class \code{"character"}; abbreviated variable name registered in the Global Soil Data registry}
    \item{\code{sd1}:}{object of class \code{"SpatialPixelsDataFrame"}; number of realizations of the target variable at depth 2.5 cm}
    \item{\code{sd2}:}{object of class \code{"SpatialPixelsDataFrame"}; number of realizations of the target variable at depth 10 cm}
    \item{\code{sd3}:}{object of class \code{"SpatialPixelsDataFrame"}; number of realizations of the target variable at depth 22.5 cm}
    \item{\code{sd4}:}{object of class \code{"SpatialPixelsDataFrame"}; number of realizations of the target variable at depth 45 cm}
    \item{\code{sd5}:}{object of class \code{"SpatialPixelsDataFrame"}; number of realizations of the target variable at depth 80 cm}
    \item{\code{sd6}:}{object of class \code{"SpatialPixelsDataFrame"}; number of realizations of the target variable at depth 150 cm}
    \item{\code{validation}:}{object of class \code{"SpatialPointsDataFrame"}; cross-validation data set with values of the target variable at six depths}
  }
}
\section{Methods}{
  \describe{
	\item{summary}{\code{signature(x = "GlobalSoilMap")}: generates summary statistics for the object}
  }
}
\author{ Tomislav Hengl and Robert A. MacMillan }
\seealso{ \code{\link{SpatialComponents-class}}
}
\examples{
\dontrun{# load soil samples from the plotKML package: 
library(plotKML)
data(eberg)
# 30 percent subsample:
eberg <- eberg[runif(nrow(eberg)) < .3,]
# sites table:
s.lst <- c("ID", "soiltype", "TAXGRSC", "X", "Y")
sites <- eberg[,s.lst]
# get horizons table:
horizons <- getHorizons(eberg, idcol="ID", sel=c("UHDICM","LHDICM","SNDMHT","SLTMHT","CLYMHT"))
# create object of type "SoilProfileCollection"
eberg.spc <- join(horizons, sites, type='inner')
depths(eberg.spc) <- ID ~ UHDICM + LHDICM
site(eberg.spc) <- as.formula(paste("~", paste(sel.site[-1], collapse="+"), sep=""))
coordinates(eberg.spc) <- ~X+Y
proj4string(eberg.spc) <- CRS("+init=epsg:31467")
# convert to geosamples:
eberg.geo <- as.geosamples(eberg.spc)
# prepare covariates:
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
eberg_spc <- spc(eberg_grid, as.formula("~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6"))
# fit gstat Model:
formulaString = as.formula(paste("observedValue ~ ", paste(names(eberg_spc@sp), collapse="+"), "+ ns(altitude, df=5)"))
SNDMHT.m <- fit.gstatModel(eberg.geo, formulaString, covariates=eberg_spc, methodid="SNDMHT")
summary(SNDMHT.m@regModel)
}
}
\keyword{classes}