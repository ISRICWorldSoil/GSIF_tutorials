# title         : tutorial_eberg.R
# purpose       : Pedometric mapping using the Ebergotzen data set;
# reference     : [http://gsif.r-forge.r-project.org/tutorial_eberg.php]
# producer      : Prepared by T. Hengl, Bas Kempen, Gerard Heuvelink
# last update   : In Wageningen, NL, Oct 2012.
# inputs        : Ebergotzen data set [http://plotkml.r-forge.r-project.org/eberg.html]; 3670 observations of soil classes and textures; 100 m and 25 resolution grids (covariates)
# outputs       : 3D predictions of soil properties and classes;


##-----------------------------------
## Load packages / data
##-----------------------------------

library(plotKML)
library(GSIF)
library(raster)

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
sel <- runif(nrow(eberg)) < .3
eberg.xy <- eberg[sel,]
coordinates(eberg.xy) <- ~X+Y
proj4string(eberg.xy) <- CRS("+init=epsg:31467")

# format gridded data:
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")

## Point pattern statistics:
library(spatstat)
# mg_owin <- as.owin(eberg_grid[1])
mg_owin <- as.owin(data.frame(x = data.frame(eberg_grid)[,"x"], y = data.frame(eberg_grid)[,"y"], window = TRUE))
eberg.ppp <- ppp(x=coordinates(eberg.xy)[,1], y=coordinates(eberg.xy)[,2], window=mg_owin)
summary(nndist(eberg.ppp))
# Complete Spatial Randomness:
env.eberg.xy <- envelope(eberg.ppp, fun=Gest)
par(mar=c(4.5,4.5,0.5,0.5), oma=c(0,0,0,0))
plot(env.eberg.xy, lwd=list(3,1,1,1), main="")
## http://gsif.r-forge.r-project.org/Fig_eberg_CRS_test.png

## MaxEnt analysis:
#jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
#if (file.exists(jar)) {
#me.eberg <- MaxEnt(occurrences=eberg.ppp, covariates=eberg_grid)
#par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,0.5), oma=c(0,0,0,0))
#image(as(me.eberg@predicted, "SpatialPixelsDataFrame"), col=rev(heat.colors(25)), xlab="", ylab="")
#points(me.eberg@occurrences, pch="+", cex=.7)
#image(me.eberg@sp.domain, col="grey", xlab="", ylab="")
#}
## http://gsif.r-forge.r-project.org/Fig_eberg_MaxEnt_test.png
# plot(me.eberg@maxent)

## Convert to SoilProfileCollection:
# list columns of interest:
s.lst <- c("ID", "soiltype", "TAXGRSC", "X", "Y")
h.lst <- c("UHDICM","LHDICM","SNDMHT","SLTMHT","CLYMHT")
# get sites table:
sites <- eberg[sel,s.lst]
# get horizons table:
horizons <- getHorizons(eberg[sel,], idcol="ID", sel=h.lst)
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
# str(eberg.geo)
levels(eberg.geo@data$methodid)
# the observationid is missing, so we make our own:
eberg.geo@data$observationid <- paste("eberg", 1:length(eberg.geo@data$observationid), sep="")

# Derive SPCs:
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
# eberg_spc@pca$rotation
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
# run a proper cross-validation:
rk.cv <- validate(SNDMHT.m)

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
SNDMHT.gsm <- GlobalSoilMap(sd.ll, varname="SNDMHT", period=c("1999-02-01", "2001-07-01"))
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
kml(SNDMHT.geo, shape=shape, colour=observedValue, zlim=c(10,85), file.name="SNDMHT_eberg.kml", altitude=z0+5000+(SNDMHT.geo@coords[,3]*2500), balloon=FALSE, labels="", extrude=FALSE, altitudeMode="relativeToGround", size=.3)

## Uncertainty at two arbitrary locations:
loc <- eberg_spc@predicted[1200:1201,]
new3D.loc <- sp3D(loc)
str(new3D.loc[[1]])
sd.loc <- predict(SNDMHT.m, predictionLocations=new3D.loc[[1]], nfold=0)
## 95% interval:
int95 <- sd.loc@predicted$var1.pred[1] + c(-1.96, 1.96)*sqrt(sd.loc@predicted$var1.var[1])
exp(int95)/(1+exp(int95))*100
new3D.loc[[1]]@coords[1,]
## Subset / stripe:
SNDMHT.xy <- spTransform(SNDMHT.geo, CRS("+init=epsg:31467"))
sel.stripe <- eberg_spc@predicted@coords[,2] > min(SNDMHT.xy@coords[,2])  # 2400 locations
loc <- eberg_spc@predicted[sel.stripe,]
new3D.loc <- sp3D(loc)
sd.loc <- predict(SNDMHT.m, predictionLocations=new3D.loc[[1]], nsim=10, block=c(0,0,0))
## TH: geostat simulations on a block support in gstat is very time-consuming;  -> not a good idea!
## back-transform the values:
sd.loc@realizations <- calc(sd.loc@realizations, fun=function(x){exp(x)/(1+exp(x))*100})
str(sd.loc, max.level=2)
plotKML(sd.loc, file.name="SNDMHT_sims.kml", zlim=c(10,85))
## http://gsif.r-forge.r-project.org/Fig_eberg_sims_cross_section.png


## Predicting soil types:
formulaString = soiltype ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
eberg_sm <- spmultinom(formulaString, eberg.xy, eberg_spc@predicted)
eberg_sm@class.c
# plot memberships:
pal = seq(0, 1, 1/50)
spplot(eberg_sm@mu, at=pal, col.regions=rev(grey(pal)))
# classes predicted:
Ls = length(levels(eberg_sm@predicted$soiltype))
pnts = list("sp.points", eberg.xy, pch="+", cex=.6, col="black")
spplot(eberg_sm@predicted, col.regions=rainbow(Ls)[rank(runif(Ls))], sp.layout=pnts)
## http://gsif.r-forge.r-project.org/Fig_eberg_Soiltypes_spfkm.png

# predict using soil types, let's try using the superfied fuzzy k-means:
eberg_sm <- spfkm(formulaString, eberg.xy, eberg_spc@predicted)
glm.formulaString2 = as.formula(paste("SNDMHT_A ~ ", paste(names(eberg_sm@mu), collapse="+"), "-1"))
glm.formulaString2
SNDMHT.m2 <- fit.gstatModel(observations=eberg.xy, glm.formulaString2, covariates=eberg_sm@mu)
summary(SNDMHT.m2@regModel)

## Predicting with multiscale data:
eberg_grids <- list(eberg_grid, eberg_grid25)
unlist(sapply(eberg_grids, names))
glm.formulaString3 = observedValue ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6+LNCCOR6+TWITOPx+NVILANx+ns(altitude, df=4)
SNDMHT.m3 <- fit.gstatModel(observations=eberg.geo, glm.formulaString3, covariates=eberg_grids, methodid="SNDMHT.t")  # this takes slightly more time...
summary(SNDMHT.m3@regModel)
# new3D2s <- sp3D(eberg_grids, stdepths=-0.025)
# sd.l2 <- predict(SNDMHT.m3, predictionLocations=new3D2s[[1]], nfold=0)
## using original covariates will often lead to problems as not all classes have been sampled; 

## instead, we can downscale the 100 m resolution grids:
eberg_grid25p <- eberg_grid25
eberg_grid25p@data <- cbind(eberg_grid25@data, gdalwarp(eberg_grid, pixsize=eberg_grid25@grid@cellsize[1], GridTopology=eberg_grid25@grid, resampling_method="cubicspline")@data)
# and now we can again run spc:
formulaString2 <- ~ TWITOPx+NVILANx+PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc25 <- spc(eberg_grid25p, formulaString2)
# rd2 = range(eberg_spc25@predicted@data[,1], na.rm=TRUE)
# spplot(eberg_spc25@predicted[1:4], at=seq(rd2[1], rd2[2], length.out=48), col.regions=pal)
## http://gsif.r-forge.r-project.org/Fig_eberg_SPCs2_4.png

## fit the model:
glm.formulaString3 = as.formula(paste("observedValue ~ ", paste(names(eberg_spc25@predicted), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString3
SNDMHT.m3 <- fit.gstatModel(observations=eberg.geo, glm.formulaString3, covariates=eberg_spc25@predicted, methodid="SNDMHT.t")
summary(SNDMHT.m3@regModel)
## Prepare prediction locations (this requires downscaling!):
new3D2s <- sp3D(eberg_spc25@predicted, stdepths=-0.025)
sd.l2 <- predict(SNDMHT.m3, predictionLocations=new3D2s[[1]], nfold=0) ## this can take 3-5 mins!
sd.l2@predicted$observedValue <- exp(sd.l2@predicted$observedValue)/(1+exp(sd.l2@predicted$observedValue))*100
## compare predictions at two scales:
data(SAGA_pal)
rg <- range(sd.l2@predicted$observedValue, na.rm=TRUE)
rx <- rev(as.character(round(c(round(rg[1], 0), NA, round(mean(rg), 0), NA, round(rg[2], 0)), 2))) 
par(mfrow=c(1,2), mar=c(.5,.5,3.5,0.5), oma=c(0,0,0,0))
image(raster(sd.l2@predicted["observedValue"]), col=SAGA_pal[[1]], main="25 m", axes = FALSE, xlab="", ylab="", zlim=rg, asp=1)
points(eberg.xy, pch="+", cex=.5)
image(raster(sd.l[[1]]@predicted["observedValue"]), col=SAGA_pal[[1]], main="100 m", axes = FALSE, xlab="", ylab="", zlim=rg, asp=1)
points(eberg.xy, pch="+", cex=.5)
## http://gsif.r-forge.r-project.org/Fig_eberg_comparison_25_100_m.png
plotKML(sd.l2@predicted["observedValue"], file.name="SNDMHT_25m.kml", zlim=c(10,85))

## Predicting with multisource data:
formulaString.l <- list(~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6, ~ DEMTOPx+TWITOPx+NVILANx)
eberg_grids_spc <- spc(eberg_grids, formulaString.l)
## fit a list of models:
glm.formulaString.l <- lapply(eberg_grids_spc, FUN=function(x){as.formula(paste("observedValue ~ ", paste(names(x@predicted), collapse="+"), "+ ns(altitude, df=4)"))})
glm.formulaString.l
## focus on the regression model (this returns a list of models):
SNDMHT.ml <- fit.gstatModel(observations=eberg.geo, glm.formulaString.l, lapply(eberg_grids_spc, slot, "predicted"), methodid="SNDMHT.t")
# summary(SNDMHT.ml[[1]]@regModel)
new3D.ml <- sapply(eberg_grids_spc, FUN=function(x){sp3D(x@predicted, stdepths=-0.025)})
## create predictions using a list of gstatModels:
sd.ml <- predict(SNDMHT.ml, predictionLocations=new3D.ml, nfold=2, verbose=TRUE, mask.extra=FALSE)
sd.mlc <- merge(sd.ml[[1]], sd.ml[[2]], silent=FALSE)
sd.mlc@data[,1] <- exp(sd.mlc@data[,1])/(1+exp(sd.mlc@data[,1]))*100
par(mar=c(.5,.5,3.5,0.5), oma=c(0,0,0,0))
image(raster(sd.mlc["observedValue"]), col=SAGA_pal[[1]], main="25 m + 100 m", axes = FALSE, xlab="", ylab="", zlim=rg, asp=1)
points(eberg.xy, pch="+", cex=.5)
## http://gsif.r-forge.r-project.org/Fig_eberg_merge_25_100_m.png


# end of script;