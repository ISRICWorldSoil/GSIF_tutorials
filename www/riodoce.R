# title         : riodoce.R
# purpose       : soil mapping and analysis for the Rio Doce case study;
# reference     : [http://gsif.r-forge.r-project.org/riodoce.html];
# producer      : Prepared by Eliana De Souza, T. Hengl & Bas Kempen
# last update   : In Wageningen, NL, June 2012.
# inputs        : to acces the riodoce data first install GSIF package [http://gsif.r-forge.r-project.org/]
# outputs       : maps and visualizations;
# remarks 1     : ;

library(GSIF)
library(plotKML)
library(aqp)
utm.csy <- "+proj=utm +zone=23 +south +ellps=GRS67 +units=m +no_defs"

## Load data
data(riodoce)
## Promote to SoilProfileCollection
riodoce.spc <- join(riodoce$horizons, riodoce$sites, type='inner')
## specify point IDs and depths:
depths(riodoce.spc) <- SOURCEID ~ UHDICM + LHDICM
## extract site data
site(riodoce.spc) <- ~ LONWGS84 + LATWGS84 + TAXBRC  
## generate SpatialPoints
coordinates(riodoce.spc) <- ~ LONWGS84 + LATWGS84
proj4string(riodoce.spc) <- "+proj=latlong +datum=WGS84"
# str(riodoce)
## convert to geosamples:
riodoce.geo <- as.geosamples(riodoce.spc)
str(riodoce.geo)
## available variables
summary(riodoce.geo@data$methodid)

# load the gridded data:
download.file("http://worldgrids.org/rda/riodoce.grids.rda", "riodoce.grids.rda")
load("riodoce.grids.rda")
gridded(riodoce.grids) <- ~x+y
proj4string(riodoce.grids) <- CRS(utm.csy)
## convert to SPCs:
formulaString = ~ TAXBRC3+DEMSRT3+SLPSRT3+TWISRT3+EV1MOD3+EV2MOD3+TD1MOD3+TN1MOD3
riodoce_spc <- spc(riodoce.grids, formulaString)
# riodoce_spc@pca$rotation
# spplot(riodoce_spc@predicted[1:3])
# export all produced maps to SAGA:
# for(j in names(riodoce_spc@predicted)){ writeGDAL(riodoce_spc@predicted[j], paste(j, ".sdat", sep=""), "SAGA", mvFlag=-99999) }
## Ignore soil polygon map:
formulaString2 = ~ DEMSRT3+SLPSRT3+TWISRT3+EV1MOD3+EV2MOD3+TD1MOD3+TN1MOD3
riodoce_spc2 <- spc(riodoce.grids, formulaString2)


#########################################
# Modeling soil types (memberships)
#########################################

TAXBRC.xy <- riodoce$sites
coordinates(TAXBRC.xy) <- ~LONWGS84+LATWGS84
proj4string(TAXBRC.xy) <- CRS("+proj=latlong +datum=WGS84")
TAXBRC.xy <- spTransform(TAXBRC.xy["TAXBRC"], CRS(utm.csy))
formulaString = as.formula(paste(paste("TAXBRC ~"), paste(names(riodoce_spc@predicted), collapse="+")))
riodoce_sm <- spfkm(formulaString, TAXBRC.xy, riodoce_spc@predicted)
## this can take >2-3 minutes!!
## prediction error 0.693
## groups ‘M’ ‘RR’ are empty
## plot memberships:
pal = seq(0, 1, 1/50)
pnts = list("sp.points", TAXBRC.xy, pch="+", col="black")
# classes predicted:
riodoce_sm@predicted$TAXBRC <- as.factor(riodoce_sm@predicted$TAXBRC)
Ls = length(levels(riodoce_sm@predicted$TAXBRC))
cols = rainbow(Ls)[rank(runif(Ls))]
spplot(riodoce_sm@mu, at=pal, col.regions=rev(grey(pal)))
names(riodoce_sm@mu)
## write predictions to GIS:
# for(j in names(riodoce_sm@mu)){ writeGDAL(riodoce_sm@mu[j], paste("mu_", j, ".sdat", sep=""), "SAGA", mvFlag=-99999) }
## plot the final classification result (Area class map):
spplot(riodoce_sm@predicted, col.regions=rainbow(Ls)[rank(runif(Ls))], sp.layout=pnts)
# image(riodoce_sm@mu[2], col=rev(grey(pal)))
# text(TAXBRC.xy@coords, paste(TAXBRC.xy$TAXBRC), cex=.6, col="black")
# confusion matrix:
# riodoce_sm@confusion
# library(vcd)
# agreementplot(riodoce_sm@confusion[-20,])

#########################################
# (WA) soil polygon map only
#########################################

## fit splines to estimate values of ORC at fixed depths?

## fit a GLM-variogram:
glm.formulaString.WA = observedValue ~ TAXBRC3 -1 + ns(altitude, df=4)
glm.formulaString.WA
ORCDRC.WA <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.WA, covariates=riodoce.grids, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.WA@regModel)
# ORCDRC.WA@vgmModel
## force a zero-nugget variogram:
ORCDRC.WA@vgmModel$psill[2] = 0
## prepare new locations and make predictions: 
new3D.WA <- sp3D(riodoce.grids["TAXBRC3"], stdepths=c(-.05, -.20))
ORCDRC.WA.sd1 <- predict(ORCDRC.WA, predictionLocations = new3D.WA[[1]])
writeGDAL(ORCDRC.WA.sd1@predicted["observedValue"], "ORCDRC_sd1_WA.sdat", "SAGA", mvFlag=-99999)
# image(ORCDRC.WA.sd1@predicted["observedValue"])
# ORCDRC.WA.sd1sim <- predict(ORCDRC.WA, predictionLocations = new3D.WA[[1]], method="RK", nsim=5)
# writeRaster(ORCDRC.WA.sd1sim@realizations[[1]], "ORCDRC_sd1sim_WA.tif", NAflag=-99999, overwrite=TRUE)

#########################################
# (SMT) using soil class memberships
#########################################

glm.formulaString.SMT = as.formula(paste("observedValue ~ ", paste(names(riodoce_sm@mu), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString.SMT
ORCDRC.SMT <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.SMT, covariates=riodoce_sm@mu, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.SMT@regModel)
## prepare new locations and make predictions: 
new3D.SMT <- sp3D(riodoce_sm@mu, stdepths=c(-.05, -.20))
ORCDRC.SMT.sd1 <- predict(ORCDRC.SMT, predictionLocations = new3D.SMT[[1]])
writeGDAL(ORCDRC.SMT.sd1@predicted["observedValue"], "ORCDRC_sd1_SMT.sdat", "SAGA", mvFlag=-99999)

#########################################
# (DMS1) Using soil mapping units
#########################################

glm.formulaString.DMS1 = as.formula(paste("observedValue ~ ", paste(names(riodoce_spc@predicted), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString.DMS1
ORCDRC.DMS1 <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.DMS1, covariates=riodoce_spc@predicted, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.DMS1@regModel)  # 39%?
# ORCDRC.m@vgmModel
# prepare new locations and make predictions: 
new3D.DMS1 <- sp3D(riodoce_spc@predicted, stdepths=c(-.05, -.20))
ORCDRC.DMS1.sd1 <- predict(ORCDRC.DMS1, predictionLocations = new3D.DMS1[[1]])
writeGDAL(ORCDRC.DMS1.sd1@predicted["observedValue"], "ORCDRC_sd1_DMS1.sdat", "SAGA", mvFlag=-99999)


#########################################
# (DMS2) ignoring the soil mapping units
#########################################

glm.formulaString.DMS2 = as.formula(paste("observedValue ~ ", paste(names(riodoce_spc2@predicted), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString.DMS2
ORCDRC.DMS2 <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.DMS2, covariates=riodoce_spc2@predicted, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.DMS2@regModel)  # %?
# prepare new locations and make predictions: 
new3D.DMS2 <- sp3D(riodoce_spc2@predicted, stdepths=c(-.05, -.20))
ORCDRC.DMS2.sd1 <- predict(ORCDRC.DMS2, predictionLocations = new3D.DMS2[[1]])
writeGDAL(ORCDRC.DMS2.sd1@predicted["observedValue"], "ORCDRC_sd1_DMS2.sdat", "SAGA", mvFlag=-99999)

#########################################
# Cross-validation results (4 methods)
#########################################

data(SAGA_pal)
r <- range(ORCDRC.DMS1.sd1@predicted$observedValue, na.rm=TRUE)
rx <- rev(as.character(round(c(round(r[1], 0), NA, round(mean(r), 0), NA, round(r[2], 0)), 2))) 
par(mfrow=c(2,2), mar=c(.5,.5,3.5,0.5), oma=c(0,0,0,0))
image(log1p(raster(ORCDRC.WA.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='WA', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8, bty="n")
image(log1p(raster(ORCDRC.SMT.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='SMT', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8, bty="n")
image(log1p(raster(ORCDRC.DMS1.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='DMS1', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8, bty="n")
image(log1p(raster(ORCDRC.DMS2.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='DMS2', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8, bty="n")

## cross-validation for 4 methods
summary(ORCDRC.DMS2.sd1)

# end of script;
