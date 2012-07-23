# title         : tutorial_eberg.R
# purpose       : Pedometric mapping using the Ebergotzen data set;
# reference     : [http://gsif.r-forge.r-project.org/tutorial_eberg.php]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Jul 2012.
# inputs        : Ebergotzen data set [http://plotkml.r-forge.r-project.org/eberg.html]; 3670 observations of soil classes and textures; 100 m and 25 resolution grids (covariates)
# outputs       : 3D predictions of soil properties and classes;


##-----------------------------------
## Load packages / data
##-----------------------------------

library(plotKML)
library(GSIF)

# load data:
data(eberg)
data(eberg_grid)
data(eberg_grid25)

str(eberg)
summary(eberg$SNDMHT_A)
library(StatDA)
par(mar=c(2.5,2.5,0.5,0.5), oma=c(0,0,0,0))
edaplot(eberg$SNDMHT_A[!is.na(eberg$SNDMHT_A)], H.freq=TRUE, box=FALSE, S.pch=3, S.cex=0.5, D.lwd=1.5, P.ylab="", P.log=FALSE, P.logfine=c(5,10), P.main="", P.xlab="", B.pch=3, B.cex=0.5)
## http://gsif.r-forge.r-project.org/Fig_eberg_hist_SNDMHT.png

# prepare data for spatial analysis:
eberg.xy <- eberg[runif(nrow(eberg)) < .3,]
coordinates(eberg.xy) <- ~X+Y
proj4string(eberg.xy) <- CRS("+init=epsg:31467")

# format gridded data:
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")

# Point pattern statistics:
library(spatstat)
mg_owin <- as.owin(eberg_grid[1])
eberg.ppp <- ppp(x=coordinates(eberg.xy)[,1], y=coordinates(eberg.xy)[,2], window=mg_owin)
summary(nndist(eberg.ppp))
# Complete Spatial Randomness:
env.eberg.xy <- envelope(eberg.ppp, fun=Gest)
par(mar=c(4.5,4.5,0.5,0.5), oma=c(0,0,0,0))
plot(env.eberg.xy, lwd=list(3,1,1,1), main="")
## http://gsif.r-forge.r-project.org/Fig_eberg_CRS_test.png

# MaxEnt analysis:
me.eberg <- MaxEnt(occurrences=eberg.ppp, covariates=eberg_grid)
par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,0.5), oma=c(0,0,0,0))
image(as(me.eberg@predicted, "SpatialPixelsDataFrame"), col=rev(heat.colors(25)), xlab="", ylab="")
points(me.eberg@occurrences, pch="+", cex=.7)
image(me.eberg@sp.domain, col="grey", xlab="", ylab="")
## http://gsif.r-forge.r-project.org/Fig_eberg_MaxEnt_test.png
# plot(me.eberg@maxent)

## Convert to SoilProfileCollection:
# list columns of interest:
s.lst <- c("ID", "soiltype", "TAXGRSC", "X", "Y")
h.lst <- c("UHDICM","LHDICM","SNDMHT","SLTMHT","CLYMHT")
# get sites table:
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
str(eberg.geo)
levels(eberg.geo@data$methodid)

# Derive SPCs:
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
eberg_spc@pca$rotation
pal = rev(rainbow(65)[1:48])
rd = range(eberg_spc@predicted@data[,1], na.rm=TRUE)
spplot(eberg_spc@predicted[1:4], at=seq(rd[1], rd[2], length.out=48), col.regions=pal)
## http://gsif.r-forge.r-project.org/Fig_eberg_SPCs1_4.png

## Build a 3D "gstatModel" 
glm.formulaString = as.formula(paste("observedValue ~ ", paste(names(eberg_spc@predicted), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString
SNDMHT.m <- fit.gstatModel(observations=eberg.geo, glm.formulaString, covariates=eberg_spc@predicted, methodid="SNDMHT.t")
summary(SNDMHT.m@regModel)
SNDMHT.m@vgmModel

## Prepare prediction locations:
new3D <- sp3D(eberg_spc@predicted)
str(new3D[[1]]@grid)
# Make predictions at six depths:
sd.l <- lapply(new3D, FUN=function(x){predict(SNDMHT.m, predictionLocations=x, nfold=0)})
# back-transform values from logits:
for(j in 1:length(sd.l)){ sd.l[[j]]@predicted$observedValue <- exp(sd.l[[j]]@predicted$observedValue)/(1+exp(sd.l[[j]]@predicted$observedValue))*100 }
# reproject to WGS84 system (100 m resolution):
p = get("cellsize", envir = GSIF.opts)[2]
s = get("stdepths", envir = GSIF.opts)
s
sd.ll <- sapply(1:length(sd.l), FUN=function(x){make.3Dgrid(sd.l[[x]]@predicted[3:4], pixelsize=p, stdepths=s[x])})
# save to a "GlobalSoilMap" object:
SNDMHT.gsm <- GlobalSoilMap(varname="SNDMHT", sd.ll, period=c("1999-02-01", "2001-07-01"))
save(SNDMHT.gsm, file="SNDMHT.rda", compress="xz") 

# visualize all maps in Google Earth:
z0 = mean(eberg_grid$DEMSRT6, na.rm=TRUE)
# export grids:
for(j in 1:length(sd.ll)){
  kml(slot(SNDMHT.gsm, paste("sd", j, sep="")), folder.name=paste("eberg_sd", j, sep=""), file=paste("SNDMHT_sd", j, ".kml", sep=""), colour=observedValue, zlim=c(10,85), raster_name=paste("SNDMHT_sd", j, ".png", sep=""), altitude=z0+5000+(s[j]*2500))
}
# export points:
SNDMHT.geo <- subset(eberg.geo, method="SNDMHT")
SNDMHT.geo <- SNDMHT.geo[SNDMHT.geo$latitude>51.57&SNDMHT.geo$latitude<51.59,]
SNDMHT.geo$observedValue <- as.numeric(SNDMHT.geo$observedValue)
coordinates(SNDMHT.geo) <- ~ longitude + latitude + altitude
proj4string(SNDMHT.geo) <- CRS("+proj=longlat +datum=WGS84")
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(SNDMHT.geo, shape=shape, colour=observedValue, zlim=c(10,85), file="SNDMHT_eberg.kml", altitude=z0+5000+(SNDMHT.geo@coords[,3]*2500), balloon=FALSE, labels="", extrude=FALSE, altitudeMode="absolute", size=.3)

## Predicting soil types:
formulaString = soiltype ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
eberg_sm <- spfkm(formulaString, eberg.xy, eberg_spc@predicted)
eberg_sm@class.c
# plot memberships:
pal = seq(0, 1, 1/50)
spplot(eberg_sm@mu, at=pal, col.regions=rev(grey(pal)))
# classes predicted:
eberg_sm@predicted$soiltype <- as.factor(eberg_sm@predicted$soiltype)
Ls = length(levels(eberg_sm@predicted$soiltype))
pnts = list("sp.points", eberg.xy, pch="+", cex=.6, col="black")
spplot(eberg_sm@predicted, col.regions=rainbow(Ls)[rank(runif(Ls))], sp.layout=pnts)
## http://gsif.r-forge.r-project.org/Fig_eberg_Soiltypes_spfkm.png

# predict using soil types:
glm.formulaString2 = as.formula(paste("SNDMHT_A ~ ", paste(names(eberg_sm@mu), collapse="+"), "-1"))
glm.formulaString2
SNDMHT.m2 <- fit.gstatModel(observations=eberg.xy, glm.formulaString2, covariates=eberg_sm@mu)
summary(SNDMHT.m2@regModel)

## Predicting with multiscale data:


## Predicting with multisource data:

# end of script;