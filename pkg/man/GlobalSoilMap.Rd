\name{GlobalSoilMap-class}
\docType{class}
\encoding{latin1}
\alias{GlobalSoilMap-class}
\alias{GlobalSoilMap}
\title{A class for GlobalSoilMap soil property maps}
\description{A class containing predictions and prediction error (or multiple realizations) of some of the target global soil property at six standard depths following the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}: \code{sd1} = 2.5 cm (0--5), \code{sd2} = 10 cm (5--15), \code{sd3} = 22.5 cm (15--30), \code{sd4} = 45 cm (30--60), \code{sd5} = 80 cm (60--100), \code{sd6} = 150 cm (100--200).}
\section{Slots}{
  \describe{    
    \item{\code{varname}:}{object of class \code{"character"}; abbreviated variable name registered in the Global Soil Data registry}
    \item{\code{TimeSpan.begin}:}{object of class \code{"POSIXct"}; begin of the sampling period}
    \item{\code{TimeSpan.end}:}{object of class \code{"POSIXct"}; end of the sampling period}
    \item{\code{sd1}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 2.5 cm (0--5)}
    \item{\code{sd2}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 10 cm (5--15)}
    \item{\code{sd3}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 22.5 cm (15--30)}
    \item{\code{sd4}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 45 cm (30--60)}
    \item{\code{sd5}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 80 cm (60--100)}
    \item{\code{sd6}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 150 cm (100--200)}
  }
Gridded data submitted to \code{sd*} slots of the \code{"GlobalSoilMap"} class must satisfy all of the following requirements (class validity):
  \itemize{
   \item All grids submitted must have the same grid topology (identical \code{grid} slot in the object of class \code{"SpatialPixelsDataFrame"});
   \item All grids must be projected in the referent coordinate system \href{http://spatialreference.org/ref/epsg/4326/}{WGS84} (geographical coordinates), with 3D dimension (altitude) expressed as distance from the land surface in meters (e.g. altitude of \code{-.025} corresponds to the 2.5 cm depth);
   \item The grid cell size must correspond to some standard resolution e.g. 0.0008333333 (1/1200 or about 100 m), 0.0016666667 (1/600 or about 250 m) etc;
   \item Only standard abbreviated names registered in the Global Soil Data registry can be used in the \code{varname} slot; 
  } 
}
\section{Methods}{
  \describe{
	\item{summary}{\code{signature(x = "GlobalSoilMap")}: generates summary statistics for the object}
  }
}
\references{
\itemize{
  \item Sanchez, P. A., S. Ahamed, F. Carre, A. E. Hartemink, J. Hempel, J. Huising, P. Lagacherie, A. B. McBratney, N. J. McKenzie, M L. de \enc{Mendonça}{Mendonca}-Santos, et al., (2009) \href{http://dx.doi.org/10.1126/science.1175084}{Digital Soil Map of the World}. Science, 325(5941): 680-681.
}
}
\author{ Tomislav Hengl and Robert A. MacMillan }
\seealso{ \code{\link{SpatialComponents-class}}, \code{\link{geosamples-class}} 
}
\examples{
# load soil samples from the plotKML package: 
library(plotKML)
library(aqp)
library(plyr)
library(splines)
library(rgdal)
library(raster)

data(eberg)
# subset data to 10%:
eberg <- eberg[runif(nrow(eberg)) < .1,]
# sites table:
s.lst <- c("ID", "soiltype", "TAXGRSC", "X", "Y")
h.lst <- c("UHDICM","LHDICM","SNDMHT","SLTMHT","CLYMHT")
sites <- eberg[,s.lst]
# get horizons table:
horizons <- getHorizons(eberg, idcol="ID", sel=h.lst)
# create object of type "SoilProfileCollection"
eberg.spc <- join(horizons, sites, type='inner')
depths(eberg.spc) <- ID ~ UHDICM + LHDICM
site(eberg.spc) <- as.formula(paste("~", paste(s.lst[-1], collapse="+"), sep=""))
coordinates(eberg.spc) <- ~X+Y
proj4string(eberg.spc) <- CRS("+init=epsg:31467")
# convert to logits:
eberg.spc@horizons$SNDMHT.t <- log((eberg.spc@horizons$SNDMHT/100)/
    (1-eberg.spc@horizons$SNDMHT/100))
# convert to geosamples:
eberg.geo <- as.geosamples(eberg.spc)
# load gridded data:
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
# derive spc's:
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
# build a 3D "gstatModel": 
glm.formulaString = as.formula(paste("observedValue ~ ", 
     paste(names(eberg_spc@predicted), collapse="+"), "+ ns(altitude, df=4)"))
\dontrun{
SNDMHT.m <- fit.gstatModel(observations=eberg.geo, glm.formulaString, 
     covariates=eberg_spc@predicted, methodid="SNDMHT.t")
summary(SNDMHT.m@regModel)
SNDMHT.m@vgmModel
# prepare new locations (6 standard depths): 
new3D <- sp3D(eberg_spc@predicted)
# Make predictions at six depths:
sd.l <- lapply(new3D, FUN=function(x){predict(SNDMHT.m, predictionLocations=x, nfold=0)})
# back-transform values from logits:
for(j in 1:length(sd.l)){ 
    sd.l[[j]]@predicted$observedValue <- exp(sd.l[[j]]@predicted$observedValue)/
      (1+exp(sd.l[[j]]@predicted$observedValue))*100 
}
# reproject to WGS84 system (100 m resolution):
p = get("cellsize", envir = GSIF.opts)[1]
s = get("stdepths", envir = GSIF.opts)
sd.ll <- sapply(1:length(sd.l), FUN=function(x){make.3Dgrid(sd.l[[x]]@predicted[3:4],
     pixsize=p, stdepths=s[x])})
# save to a "GlobalSoilMap" object:
SNDMHT.gsm <- GlobalSoilMap(varname="SNDMHT", sd.ll, period=c("1999-02-01", "2001-07-01"))
str(SNDMHT.gsm, max.level=2)
# visualize all maps in Google Earth:
data(R_pal)
z0 = mean(eberg_grid$DEMSRT6, na.rm=TRUE)
# export grids:
for(j in 1:length(sd.ll)){
  kml(slot(SNDMHT.gsm, paste("sd", j, sep="")), folder.name = paste("eberg_sd", j, sep=""),
     file = paste("SNDMHT_sd", j, ".kml", sep=""), colour = observedValue, z.lim=c(10,85),
     raster_name = paste("SNDMHT_sd", j, ".png", sep=""), altitude = z0+5000+(s[j]*2500))
}
}
}
\keyword{classes}