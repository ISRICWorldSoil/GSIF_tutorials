# title         : demo_riodoce.R
# purpose       : soil mapping and analysis for the Rio Doce case study;
# reference     : [http://gsif.r-forge.r-project.org/riodoce.html];
# producer      : Prepared by Eliana De Souza & T. Hengl
# last update   : In Wageningen, NL, June 2012.
# inputs        : to acces the riodoce data first install GSIF package [http://gsif.r-forge.r-project.org/]
# outputs       : maps and visualizations;
# remarks 1     : ;

# basic definitions:
library(GSIF)
library(plotKML)
library(raster)
library(aqp)
utm.csy <- "+proj=utm +zone=23 +south +ellps=GRS67 +units=m +no_defs"

## Load data
data(riodoce)
## Promote to SoilProfileCollection
riodoce.spc <- join(riodoce$horizons, riodoce$sites, type='inner')
## specify point IDs and depths:
depths(riodoce.spc) <- SOURCEID ~ UHDICM + LHDICM
## extract site data
site(riodoce.spc) <- ~ LONWGS84 + LATWGS84 + TAXBRC + legend  
## generate SpatialPoints
coordinates(riodoce.spc) <- ~ LONWGS84 + LATWGS84
proj4string(riodoce.spc) <- "+proj=latlong +datum=WGS84"
# str(riodoce)

## convert to geosamples:
riodoce.geo <- as.geosamples(riodoce.spc)
str(riodoce.geo)
## available variables
summary(riodoce.geo@data$methodid)
## write to a shape file:
riodoce.shp <- as.data.frame(riodoce.spc)
coordinates(riodoce.shp) <- ~ LONWGS84 + LATWGS84
proj4string(riodoce.shp) <- "+proj=latlong +datum=WGS84"
writePointsShape(spTransform(riodoce.shp, CRS(utm.csy)), "riodoce.shp")
# riodoce_PVe <- riodoce.shp[riodoce.shp$legend=="PVe",c("legend", "TAXBRC")]
# writePointsShape(spTransform(riodoce_PVe, CRS(utm.csy)), "riodoce_PVe.shp")

## load the gridded data:
download.file("http://worldgrids.org/rda/riodoce.grids.rda", "riodoce.grids.rda")
load("riodoce.grids.rda")
gridded(riodoce.grids) <- ~x+y
proj4string(riodoce.grids) <- CRS(utm.csy)
## derive Soil Predictive Components:
formulaString = ~ TAXBRC3+DEMSRT3+SLPSRT3+TWISRT3+EV1MOD3+EV2MOD3+TD1MOD3+TN1MOD3
riodoce_spc <- spc(riodoce.grids, formulaString)
# riodoce_spc@pca$rotation
# spplot(riodoce_spc@predicted[1:3])
# export all produced maps to SAGA:
# for(j in names(riodoce_spc@predicted)){ writeGDAL(riodoce_spc@predicted[j], paste(j, ".sdat", sep=""), "SAGA", mvFlag=-99999) }
## SPCs without the soil polygon map:
formulaString2 = ~ DEMSRT3+SLPSRT3+TWISRT3+EV1MOD3+EV2MOD3+TD1MOD3+TN1MOD3
riodoce_spc2 <- spc(riodoce.grids, formulaString2)

writeGDAL(riodoce.grids["mask"], "riodoce_mask.sdat", "SAGA", mvFlag=0)
rsaga.geoprocessor(lib="shapes_grid", module=6, param=list(GRID="riodoce_mask.sgrd", POLYGONS="riodoce_mask.shp", CLASS_ALL=0))
## borders of RioDoce:
riodoce.pol <- as(readShapePoly("riodoce_mask.shp"), "SpatialLines")

#########################################
# Modeling soil types (memberships)
#########################################

TAXBRC.xy <- riodoce$sites
coordinates(TAXBRC.xy) <- ~LONWGS84+LATWGS84
proj4string(TAXBRC.xy) <- CRS("+proj=latlong +datum=WGS84")
TAXBRC.xy <- spTransform(TAXBRC.xy["legend"], CRS(utm.csy))
## plot the study area:
alt <- raster(riodoce.grids["DEMSRT3"])
slope <- terrain(alt, opt='slope')
aspect <- terrain(alt, opt='aspect')
hill <- hillShade(slope, aspect, 40, 270)
par(mar=c(.5,.5,.5,.5), oma=c(0,0,0,0))
image(hill, col=grey(0:100/100), main='', axes=FALSE, asp=1)
lines(riodoce.pol)
points(TAXBRC.xy, pch=21, bg="white", cex=1, col="black")
dev.off()
## Fig_riodoce_study_area.pdf

formulaString.legend = as.formula(paste(paste("legend ~"), paste(names(riodoce_spc@predicted), collapse="+")))
riodoce_sm <- spmultinom(formulaString.legend, TAXBRC.xy, riodoce_spc@predicted)
## this can take >2-3 minutes!!
## Cohen weighted Kappa: 43%
## some groups are empty
pal = seq(0, 1, 1/50)
pnts = list("sp.points", TAXBRC.xy, pch="+", col="black")
# classes predicted:
riodoce_sm@predicted$legend <- as.factor(riodoce_sm@predicted$legend)
Ls = length(levels(riodoce_sm@predicted$legend))
cols = rainbow(Ls)[rank(runif(Ls))]
## plot memberships:
spplot(riodoce_sm@mu, at=pal, col.regions=rev(grey(pal)))
## Fig_RioDoce_memberships.png
names(riodoce_sm@mu)
## plot the final classification result (Area class map):
spplot(riodoce_sm@predicted, col.regions=rainbow(Ls)[rank(runif(Ls))], sp.layout=pnts)
## Fig_RioDoce_best_class.png

## confusion matrix:
# riodoce_sm@confusion
# library(vcd)
# agreementplot(riodoce_sm@confusion[-20,])
## Export to a GIS:
# riodoce_sm@predicted$legend.img <- as.integer(riodoce_sm@predicted$legend)
# riodoce_sm@predicted$legend.img <- ifelse(riodoce.grids$mask==0, NA, riodoce_sm@predicted$legend.img)
# writeGDAL(riodoce_sm@predicted["legend.img"], "mu_predicted.sdat", "SAGA", mvFlag=-99999)


#########################################
# (WA) soil polygon map only
#########################################

## ORCDRC distribution:
# xyplot(cbind(UHDICM,LHDICM) ~ ORCDRC, data=riodoce.spc@horizons, id=riodoce.spc@horizons$SOURCEID, panel=panel.depth_function, ylim=c(250,-10), scales=list(y=list(tick.number=10)), xlab='Soil organic carbon', ylab='Depth (cm)', main='')
plot(y=(riodoce.spc@horizons$LHDICM+riodoce.spc@horizons$UHDICM)/2, x=riodoce.spc@horizons$ORCDRC, type="p", ylim=c(250,-10), ylab='depth (cm)', xlab='', col=grey(.2), main='Organic carbon (promilles)', log=c("x", "y")) 
## Fig_log_log_ORCDRC.pdf
## fit splines to estimate values of ORC at fixed depths?

## fit a GLM-variogram:
glm.formulaString.WA = observedValue ~ TAXBRC3 -1 + log(altitude+100)
glm.formulaString.WA
ORCDRC.WA <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.WA, covariates=riodoce.grids, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.WA@regModel)
# ORCDRC.WA@vgmModel
## force a zero-nugget variogram:
ORCDRC.WA@vgmModel$psill[2] = 0
## prepare new locations and make predictions: 
new3D.WA <- sp3D(riodoce.grids[riodoce.grids$mask==1,"TAXBRC3"], stdepths=c(-.05, -.20))
ORCDRC.WA.sd1 <- predict(ORCDRC.WA, predictionLocations = new3D.WA[[1]], method="RK", nfold=0)
writeGDAL(ORCDRC.WA.sd1@predicted["observedValue"], "ORCDRC_sd1_WA.sdat", "SAGA", mvFlag=-99999)
# ORCDRC.WA.sd1sim <- predict(ORCDRC.WA, predictionLocations = new3D.WA[[1]], method="RK", nsim=5)
# writeRaster(ORCDRC.WA.sd1sim@realizations[[1]], "ORCDRC_sd1sim_WA.tif", NAflag=-99999, overwrite=TRUE)
ORCDRC.WA.cv <- validate(ORCDRC.WA)

#########################################
# (SMT) using soil class memberships
#########################################

glm.formulaString.SMT = as.formula(paste("observedValue ~ ", paste(names(riodoce_sm@mu), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString.SMT
ORCDRC.SMT <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.SMT, covariates=riodoce_sm@mu, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.SMT@regModel)
## prepare new locations and make predictions: 
new3D.SMT <- sp3D(riodoce_sm@mu[riodoce.grids$mask==1,], stdepths=c(-.05, -.20))
ORCDRC.SMT.sd1 <- predict(ORCDRC.SMT, predictionLocations = new3D.SMT[[1]])
writeGDAL(ORCDRC.SMT.sd1@predicted["observedValue"], "ORCDRC_sd1_SMT.sdat", "SAGA", mvFlag=-99999)
ORCDRC.SMT.cv <- validate(ORCDRC.SMT)

## predict for a single profile:
pnt <- spTransform(riodoce.shp[riodoce.shp$SOURCEID=="UFV277",], CRS(utm.csy))
ov.pnt <- overlay(riodoce.grids, pnt)
new3D.SMT.pnt <- sp3D(riodoce_sm@mu[c(ov.pnt, ov.pnt+1),], stdepths=seq(from=0,to=-2, by=-0.01), stsize = rep(0.01, 201))
ORCDRC.SMT.pnt <- lapply(new3D.SMT.pnt, FUN=function(x){predict(ORCDRC.SMT, predictionLocations = x, nfold=0)})
## plot measured vs predicted:
x.l <- sapply(ORCDRC.SMT.pnt, FUN=function(x){x@predicted@data$observedValue[1]})
plot(y=seq(from=0,to=2, by=0.01)*100, x=x.l, type="l", ylim=c(200,-10), ylab='depth (cm)', xlab='', col=grey(.2), main='Predicted OC (promilles)', log=c("x", "y"), xlim=c(1,80))
points(y=(riodoce.spc@horizons[riodoce.spc@horizons$SOURCEID=="UFV277",]$LHDICM+riodoce.spc@horizons[riodoce.spc@horizons$SOURCEID=="UFV277",]$UHDICM)/2, x=riodoce.spc@horizons[riodoce.spc@horizons$SOURCEID=="UFV277",]$ORCDRC)

#########################################
# (DMS1) Using soil mapping units
#########################################

glm.formulaString.DMS1 = as.formula(paste("observedValue ~ ", paste(names(riodoce_spc@predicted), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString.DMS1
ORCDRC.DMS1 <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.DMS1, covariates=riodoce_spc@predicted, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.DMS1@regModel)  # 39%?
# ORCDRC.m@vgmModel
# prepare new locations and make predictions: 
new3D.DMS1 <- sp3D(riodoce_spc@predicted[riodoce.grids$mask==1,], stdepths=c(-.05, -.20))
ORCDRC.DMS1.sd1 <- predict(ORCDRC.DMS1, predictionLocations = new3D.DMS1[[1]])
writeGDAL(ORCDRC.DMS1.sd1@predicted["observedValue"], "ORCDRC_sd1_DMS1.sdat", "SAGA", mvFlag=-99999)
ORCDRC.DMS1.cv <- validate(ORCDRC.DMS1)


#########################################
# (DMS2) ignoring the soil mapping units
#########################################

glm.formulaString.DMS2 = as.formula(paste("observedValue ~ ", paste(names(riodoce_spc2@predicted), collapse="+"), "+ ns(altitude, df=4)"))
glm.formulaString.DMS2
ORCDRC.DMS2 <- fit.gstatModel(observations=riodoce.geo, glm.formulaString.DMS2, covariates=riodoce_spc2@predicted, methodid="ORCDRC", family=gaussian(log))
summary(ORCDRC.DMS2@regModel)  # %?
# prepare new locations and make predictions: 
new3D.DMS2 <- sp3D(riodoce_spc2@predicted[riodoce.grids$mask==1,], stdepths=c(-.05, -.20))
ORCDRC.DMS2.sd1 <- predict(ORCDRC.DMS2, predictionLocations = new3D.DMS2[[1]])
writeGDAL(ORCDRC.DMS2.sd1@predicted["observedValue"], "ORCDRC_sd1_DMS2.sdat", "SAGA", mvFlag=-99999)
ORCDRC.DMS2.cv <- validate(ORCDRC.DMS2)


#########################################
# Summary statistics / visualization (4 methods)
#########################################

data(SAGA_pal)
r <- range(ORCDRC.DMS1.sd1@predicted$observedValue, na.rm=TRUE)
rx <- rev(as.character(round(c(round(r[1], 0), NA, round(mean(r), 0), NA, round(r[2], 0)), 2))) 
par(mfrow=c(2,2), mar=c(.5,.5,3.5,0.5), oma=c(0,0,0,0))
image(log1p(raster(ORCDRC.WA.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='WA', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.6, bty="n")
image(log1p(raster(ORCDRC.SMT.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='SMT', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.6, bty="n")
image(log1p(raster(ORCDRC.DMS1.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='DMS1', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.6, bty="n")
image(log1p(raster(ORCDRC.DMS2.sd1@predicted["observedValue"])), col=SAGA_pal[[1]], main='DMS2', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch="+", cex=.5)
legend("topleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.6, bty="n")
dev.off()
## Fig_RioDoce_organic_carbon.png


## cross-validation for 4 methods
1-var(ORCDRC.WA.cv[[1]]$residual)/var(ORCDRC.WA.cv[[1]]$observed); mean(ORCDRC.WA.cv[[1]]$zscore)
1-var(ORCDRC.SMT.cv[[1]]$residual)/var(ORCDRC.SMT.cv[[1]]$observed); mean(ORCDRC.SMT.cv[[1]]$zscore)
1-var(ORCDRC.DMS1.cv[[1]]$residual)/var(ORCDRC.DMS1.cv[[1]]$observed); mean(ORCDRC.DMS1.cv[[1]]$zscore)
1-var(ORCDRC.DMS2.cv[[1]]$residual)/var(ORCDRC.DMS2.cv[[1]]$observed); mean(ORCDRC.DMS2.cv[[1]]$zscore)

## Zoom into some smaller area:
sel <- ORCDRC.DMS1.sd1@predicted@coords[,1] < 768000 & ORCDRC.DMS1.sd1@predicted@coords[,1] > 708000 & ORCDRC.DMS1.sd1@predicted@coords[,2] < 7749000 & ORCDRC.DMS1.sd1@predicted@coords[,2] > 7689000
sel.z <- riodoce.grids@coords[,1] < 768000 & riodoce.grids@coords[,1] > 708000 & riodoce.grids@coords[,2] < 7749000 & riodoce.grids@coords[,2] > 7689000

# subset:
grd.sub0 <- ORCDRC.SMT.sd1@predicted[which(sel),"observedValue"]
grd.sub0 <- data.frame(grd.sub0)
gridded(grd.sub0) <- ~longitude + latitude
grd.sub1 <- ORCDRC.DMS1.sd1@predicted[which(sel),"observedValue"]
grd.sub1 <- data.frame(grd.sub1)
gridded(grd.sub1) <- ~longitude + latitude
grd.sub2 <- ORCDRC.DMS2.sd1@predicted[which(sel),"observedValue"]
grd.sub2 <- data.frame(grd.sub2)
gridded(grd.sub2) <- ~longitude + latitude
## Comparison using 60 x 60 pixels only:
r.s <- range(grd.sub0@data, na.rm=TRUE)
par(mfrow=c(1,3), mar=c(.5,.5,3.5,.5), oma=c(0,0,0,0))
image(log1p(raster(grd.sub0)), col=SAGA_pal[[1]], main='SMT', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch=21, bg="white", cex=1.5, lwd=1)
image(log1p(raster(grd.sub1)), col=SAGA_pal[[1]], main='DMS1', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch=21, bg="white", cex=1.5, lwd=1)
image(log1p(raster(grd.sub2)), col=SAGA_pal[[1]], main='DMS2', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
points(TAXBRC.xy, pch=21, bg="white", cex=1.5, lwd=1)
##  Fig_RioDoce_organic_carbon_subset.png

## Perspective view plot:
library(fields)
library(colorspace)
z.sub <- riodoce.grids[which(sel.z),"DEMSRT3"]
z.sub <- data.frame(z.sub)
gridded(z.sub) <- ~x+y
DMS1.val <- as.matrix(grd.sub[1])[,grd.sub@grid@cells.dim[2]:1]
z.val <- as.matrix(z.sub)[,z.sub@grid@cells.dim[2]:1]
x=seq(from=grd.sub@bbox[1,1], length.out=nrow(z.val), by=1000)
y=seq(from=grd.sub@bbox[2,1], length.out=ncol(z.val), by=1000)
zcol <- drape.color(log1p(DMS1.val), col=SAGA_pal[[1]], zlim=c(2,5.5))
z.3d <- persp(x, y, z=z.val, theta=-90, col=zcol$color.index, phi=45, scale=FALSE, ltheta=-165, shade=0.75, axes=TRUE, box=TRUE, main="Organic carbon", expand=5)

# end of script;
