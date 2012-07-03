\name{SpatialComponents-class}
\docType{class}
\alias{SpatialComponents-class}
\alias{extract,SpatialComponents}
\title{A class for gridded components derived using the \code{spc} method}
\description{A class containing a list of gridded components and results of principal component analysis.}
\section{Slots}{
  \describe{    
    \item{\code{predicted}:}{object of class \code{"SpatialPixelsDataFrame"}; predicted values for components}
    \item{\code{pca}:}{object of class \code{"list"}; output objects from the \code{stats::prcomp} process --- contains objects: 'stdev', 'rotation', 'center' and 'scale'}
  }
}
\section{Methods}{
  \describe{
	\item{extract}{\code{signature(x = "SpatialComponents", y = "SpatialPoints", \dots)}: extracts data from a list of components and creates a regression matrix (\code{"data.frame"})}
  }
}
\author{ Tomislav Hengl }
\seealso{ \code{raster::extract}
}
\examples{ 
\dontrun{# load soil samples from the plotKML package:
library(plotKML)
data(eberg)
# 20 percent subsample:
eberg <- eberg[runif(nrow(eberg)) < .2,]
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
# overlay with SPCs (create the regression matrix):
SNDMHT.eberg <- extract(eberg_spc, eberg.geo, methodid="SNDMHT")
# fit a glm:
SNDMHT.glm <- glm(observedValue~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+ns(altitude, df=4), SNDMHT.eberg, family=gaussian())
summary(SNDMHT.glm)
}
}
\keyword{classes}
