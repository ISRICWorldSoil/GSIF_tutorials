\name{makeGstatCmd}
\alias{makeGstatCmd}
\title{Make a gstat command script}
\description{Generates a command script based on the regression model and variogram. This can then be used to run predictions/simulations by using the pre-compiled binary \code{gstat.exe}.}
\usage{
makeGstatCmd(formString, vgmModel, outfile, easfile, 
        nsim = 0, nmin = 20, nmax = 40, radius, zmap = 0, 
        predictions = "var1.pred.hdr", variances = "var1.svar.hdr", 
        xcol = 1, ycol = 2, zcol = 3, vcol = 4, Xcols)
}
\arguments{
  \item{formString}{object of class \code{"formula"} --- regression model} 
  \item{vgmModel}{object of class \code{"vgmmodel"} or \code{"data.frame"}} 
  \item{outfile}{character; output file for the command script}
  \item{easfile}{character; file name for the GeoEAS file with observed values}
  \item{nsim}{integer; number of simulations}
  \item{nmin}{integer; smallest number of points in the search radius (see gstat user's manual)}
  \item{nmax}{integer; largest number of points in the search radius (see gstat user's manual)}
  \item{radius}{numeric; search radius (see gstat user's manual)}
  \item{zmap}{numeric; fixed value for the 3D dimension in the case of 3D kriging}
  \item{predictions}{character; output file name for predictions}
  \item{variances}{character; output file name for kriging variances}
  \item{xcol}{integer; position of the x column in the GeoEAS file}
  \item{ycol}{integer; position of the y column in the GeoEAS file}
  \item{zcol}{integer; position of the z column in the GeoEAS file}
  \item{vcol}{integer; position of the target variable column in the GeoEAS file}
  \item{Xcols}{integer; column numbers for the list of covariates}
}
\details{To run the script under Windows OS you need to obtain the pre-compiled \code{gstat.exe} program from the www.gstat.org website, and put it in some directory e.g. \code{c:/gstat/}. Then add the program to your path (see environmental variable under Windows > Control panel > System > Advanced > Environmental variables), or copy the exe program directly to some windows system directory.}
\note{The advantage of using \code{gstat.exe} is that it loads large grids much faster to memory than if you use gstat in R, hence it is potentially more suited for computing with large grids. The draw back is that you can only pass simple linear regression models to gstat.exe. The stand-alone gstat is not maintained by the author of gstat any more.}
\author{ Tomislav Hengl }
\references{
\itemize{
\item Bivand, R.S., Pebesma, E.J., and \enc{Gómez}{Gomez}-Rubio, V., (2008) \href{http://www.asdar-book.org/}{Applied Spatial Data Analysis with R}. Springer, 378 p. 
\item Pebesma, E., (2003) \href{http://www.gstat.org/gstat.pdf}{Gstat user's manual}. Dept. of Physical Geography, Utrecht University, p. 100, www.gstat.org 
}
}
\seealso{ \code{\link{write.data}}, \code{\link{fit.gstatModel}}, \code{gstat::krige} }
\examples{
# universal kriging with Meuse data:
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
\dontrun{# fit a model:
library(GSIF)
omm <- fit.gstatModel(observations = meuse, formulaString = om~dist, 
  family = gaussian(log), covariates = meuse.grid)
str(omm@vgmModel)
# write the regression matrix to GeoEAS:
meuse$log_om <- log1p(meuse$om)
write.data(obj=meuse, covariates=meuse.grid["dist"], outfile="meuse.eas", methodid="log_om")
writeGDAL(meuse.grid["dist"], "dist.rst", drivername="RST", mvFlag="-99999")
makeGstatCmd(log_om~dist, vgmModel=omm@vgmModel, outfile="meuse_om_sims.cmd", easfile="meuse.eas", nsim=50, nmin=20, nmax=40, radius=1500)
# compare the processing times:
system.time(system("gstat meuse_om_sims.cmd"))
vgmModel = omm@vgmModel
class(vgmModel) <- c("variogramModel", "data.frame")
system.time(om.rk <- krige(log_om~dist, meuse[!is.na(meuse$log_om),], meuse.grid, nmin=20, nmax=40, model=vgmModel, nsim=50))
## 3D kriging:
data(eberg)
# subset data:
eberg <- eberg[runif(nrow(eberg)) < .3,]
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
eberg.spc@horizons$SNDMHT.t <- log((eberg.spc@horizons$SNDMHT/100)/(1-eberg.spc@horizons$SNDMHT/100))
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
glm.formulaString = as.formula(paste("observedValue ~ ", paste(names(eberg_spc@predicted), collapse="+"), "+ ns(altitude, df=4)"))
SNDMHT.m <- fit.gstatModel(observations=eberg.geo, glm.formulaString, covariates=eberg_spc, methodid="SNDMHT.t")
# prepare new locations and make predictions: 
new3D <- sp3D(eberg_spc@predicted, stdepths=c(-.05, -.20))
# make simulations (skip cross-validation):
system.time(SNDMHT.sd1 <- predict(SNDMHT.m, predictionLocations = new3D[[1]], nfold = 0, nmin=40, nmax=50))
# same thing via the gstat.exe:
masks = all.vars(SNDMHT.m@regModel$formula)[-1]
for(i in masks){ writeGDAL(eberg_spc@predicted[i], paste(i, ".rst", sep=""), drivername="RST", mvFlag="-99999") }
# altitude does not exist as a map!
formString = SNDMHT.t ~ PC1+PC2+PC4+PC5+PC6+PC7+PC8+PC9+PC10
masks = all.vars(formString)[-1]
write.data(obj=eberg.geo, covariates=eberg_spc@predicted[masks], outfile="eberg.eas", methodid="SNDMHT.t")
makeGstatCmd(formString, vgmModel=SNDMHT.m@vgmModel, outfile="eberg_SNDMHT_rk.cmd", easfile="eberg.eas", nmin=40, nmax=50, radius=2500, zmap=-.05)
# compare the processing times:
system.time(system("gstat eberg_SNDMHT_rk.cmd"))
## Conclusion: gstat.exe takes more time with smaller maps and can result in artifacts;
## Note also that we can not pass the "+ ns(altitude, df=4)" model to gstat.exe
}
}
\keyword{methods}
