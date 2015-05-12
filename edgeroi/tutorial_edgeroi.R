# title         : tutorial_edgeroi.R
# purpose       : Mapping soil organic carbon for the Edgeroi area using 3D regression-kriging;
# reference     : [https://code.google.com/p/gsif/source/browse/trunk/edgeroi/]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : The Edgeroi [http://plotkml.r-forge.r-project.org/edgeroi.html] is one of the standard soil data sets used to test soil mapping methods in Australia. It contains 359 soil profiles with soil observations in sites and horizons tables;
# outputs       : 3D predictions of soil properties and classes;

#install.packages("GSIF", repos="http://R-Forge.R-project.org", type="source")
#install.packages("plotKML", repos="http://R-Forge.R-project.org", type="source")
library(aqp)
library(splines)
library(GSIF)
library(plyr)
library(splines)
library(plotKML)
library(rgdal)
library(spatstat)
library(maptools)

## load the data:
data(edgeroi)

## simple example with Clay content (2D):
edgeroi.spc <- join(edgeroi$sites, edgeroi$horizons, type='inner')
h1 <- edgeroi.spc[edgeroi.spc$LSQINT==1,c("SOURCEID","LONGDA94","LATGDA94","ORCDRC","PHIHO5","CLYPPT")]
coordinates(h1) <- ~ LONGDA94 + LATGDA94
proj4string(h1) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
h1.xy <- spTransform(h1, CRS("+init=epsg:28355"))
## load the 250 m grids:
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids.rda")
load(con)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
m <- fit.gstatModel(h1.xy, CLYPPT~DEMSRT5+TWISRT5+PMTGEO5+EV1MOD5+EV2MOD5+EV3MOD5, edgeroi.grids)
m@vgmModel
rk <- predict(m, edgeroi.grids)
plot(rk)
## or in one line...
rk0 <- autopredict(h1.xy["CLYPPT"], edgeroi.grids)
plot(rk0)

## 3D soil mapping:
data(edgeroi)
edgeroi.spc <- join(edgeroi$horizons, edgeroi$sites, type='inner')
depths(edgeroi.spc) <- SOURCEID ~ UHDICM + LHDICM
site(edgeroi.spc) <- ~ LONGDA94 + LATGDA94 + TAXGAUC + NOTEOBS
coordinates(edgeroi.spc) <- ~LONGDA94+LATGDA94
proj4string(edgeroi.spc) <- CRS("+proj=longlat +ellps=GRS80 +datum=WGS84 +no_defs")
edgeroi.geo <- as.geosamples(edgeroi.spc)

## Q1: linear regression ORCDRC ~ altitude:
ORCDRC.geo <- subset(edgeroi.geo, "ORCDRC")
ORCDRC.geo$observedValue <- as.numeric(ORCDRC.geo$observedValue)
summary(lm(log1p(observedValue)~altitude, ORCDRC.geo))
plot(log1p(observedValue)~altitude, ORCDRC.geo)
## A1: 52% in log-scale

## fit model and make predictions:
glm.formulaString <- as.formula(paste(paste("log1p(ORCDRC) ~ "), paste(names(edgeroi.grids), collapse="+"), paste("+ ns(altitude, df=4)")))
glm.formulaString
ORCDRC.m <- fit.gstatModel(observations=edgeroi.geo, glm.formulaString, edgeroi.grids)

## Q2: Which covariate is the best predictor of organic carbon content?
summary(ORCDRC.m@regModel)
## A2: "altitude", then "PMTGEO5Qrt/Tv", "PMTGEO5Tv"

## Q3: Amoung of the original variance in soil organic carbon content explained by this geostatical model?
cv_glm <- validate(ORCDRC.m, nfold=2)
tvar <- 1-var(cv_glm[[1]]$residual, na.rm=T)/var(cv_glm[[1]]$observed, na.rm=T)
signif(tvar*100, 3)
## A3: 75% of variation in log-scale based on cross-validation

## Q4: GLM vs randomForest:
ORCDRC.m2 <- fit.gstatModel(observations=edgeroi.geo, log1p(ORCDRC) ~ DEMSRT5 + TWISRT5 + PMTGEO5 + EV1MOD5 + EV2MOD5 + EV3MOD5 + altitude, edgeroi.grids, method="rpart")
var(residuals(ORCDRC.m2@regModel))
var(residuals(ORCDRC.m@regModel))
## A4: GLM is better

## save model:
save(ORCDRC.m, file="ORCDRC.m.rda")

edgeroi.grids$PMTGEO5[edgeroi.grids$PMTGEO5 == "Ts"] <- "Qrs"
new3D <- sp3D(edgeroi.grids)
sd.l <- lapply(new3D, FUN=function(x){predict(ORCDRC.m, predictionLocations=x, nfold=0)})

## plot results in GE:
s = get("stdepths", envir = GSIF.opts)
z0 = 300
zlim = range(sapply(sd.l, function(x){range(x@predicted@data[,"ORCDRC"], na.rm=TRUE)}))
for(j in 1:length(s)){
kml(sd.l[[j]]@predicted, folder.name = paste("edgeroi_sd", j, sep=""),
 file.name = paste("ORCDRC_sd", j, ".kml", sep=""), colour = ORCDRC, z.lim=zlim,
 colour_scale = SAGA_pal[[1]], raster_name = paste("ORCDRC_sd", j, ".png", sep=""), 
 altitude = z0+5000+(s[j]*2500),
 png.width=edgeroi.grids@grid@cells.dim[1]*5, png.height=edgeroi.grids@grid@cells.dim[2]*5)
}

## sampled values:
ORCDRC.geo <- subset(edgeroi.geo, method="ORCDRC")
ORCDRC.geo$observedValue <- as.numeric(ORCDRC.geo$observedValue)
coordinates(ORCDRC.geo) <- ~ longitude + latitude + altitude
proj4string(ORCDRC.geo) <- CRS("+proj=longlat +datum=WGS84")
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(ORCDRC.geo, shape=shape, colour=log1p(observedValue), colour_scale=SAGA_pal[[1]], z.lim=zlim, file.name="ORCDRC_edgeroi.kml", altitude=z0+5000+(ORCDRC.geo@coords[,3]*2500), balloon=FALSE, labels="", extrude=FALSE, altitudeMode="relativeToGround", size=.3)

## plot as polygons (this takes time!!)
#x <- grid2poly(sd.l[[1]]@predicted["ORCDRC"])
#kml(x, colour=ORCDRC, colour_scale = SAGA_pal[[1]])

## end of script;