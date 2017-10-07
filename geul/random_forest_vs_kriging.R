## Comparison RF vs kriging, see slides at: https://github.com/ISRICWorldSoil/GSIF_tutorials/blob/master/geul/5H_Hengl.pdf
## tom.hengl@isric.org

list.of.packages <- c("plyr", "parallel", "randomForest", "quantregForest", "plotKML", "GSIF", "ranger", "RCurl", "raster", "rgdal", "geoR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

setwd("~/git/GSIF_tutorials/geul")
load(".RData")
library(GSIF)
library(rgdal)
library(raster)
library(gstat)
library(randomForest)
library(quantregForest)
library(plotKML)
library(scales)
library(ranger)
library(RCurl)
library(parallel)
library(geoR)
library(geostatsp)
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")
## Load the Meuse data set:
demo(meuse, echo=FALSE)

## Compare GLM vs RF ----
m <- glm(zinc~log1p(dist)+ffreq, meuse, family=gaussian(link=log))
plot(m$fitted.values~m$y, asp=1)
abline(0,1)
rf <- quantregForest(x=meuse@data[,c("dist","ffreq")], y=meuse$zinc)
plot(rf$predicted~rf$y, asp=1)
abline(0,1)
meuse.grid$glm.zinc <- predict(m, meuse.grid@data, type="response")
meuse.grid$rf.zinc <- predict(rf, meuse.grid@data)[,2]

## Plot predictions next to each other:
meuse.grid$glm.zinc = ifelse(meuse.grid$glm.zinc<expm1(4.8), expm1(4.8), meuse.grid$glm.zinc)
png(file = "Fig_comparison_GLM_RF_zinc_meuse.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(meuse.grid["glm.zinc"])), col=leg, zlim=c(4.8,7.4), main="GLM")
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["rf.zinc"])), col=leg, zlim=c(4.8,7.4), main="Random Forest")
points(meuse, pch="+")
dev.off()
## TH: Very similar

## Zinc predicted using ordinary kriging ----
zinc.geo <- as.geodata(meuse["zinc"])
#plot(variog4(zinc.geo, lambda=0, max.dist=1500, messages=FALSE), lwd=2)
zinc.vgm <- likfit(zinc.geo, lambda=0, messages=FALSE, ini=c(var(log1p(zinc.geo$data)),500), cov.model="exponential")
locs = meuse.grid@coords
zinc.ok <- krige.conv(zinc.geo, locations=locs, krige=krige.control(obj.model=zinc.vgm))
meuse.grid$zinc_ok = zinc.ok$predict
meuse.grid$zinc_ok_var = zinc.ok$krige.var

## Zinc predicted using RF and buffer distances only ----
grid.dist0 <- buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
dn0 <- paste(names(grid.dist0), collapse="+")
fm0 <- as.formula(paste("zinc ~ ", dn0))
ov.zinc <- over(meuse["zinc"], grid.dist0)
m.zinc <- ranger(fm0, cbind(meuse@data["zinc"], ov.zinc), keep.inbag = TRUE)
zinc.rfd <- predict(m.zinc, grid.dist0@data, type = "se")
meuse.grid$zinc_rfd = zinc.rfd$predictions
meuse.grid$zinc_rfd_var = zinc.rfd$se

## Plot predictions next to each other:
#png(file = "Fig_comparison_OK_RF_zinc_meuse.png", res = 150, width = 1750, height = 1200)
var.max = max(c(meuse.grid$zinc_rfd_var, sqrt(meuse.grid$zinc_ok_var)))
axis.ls = list(at=c(4.8,5.7,6.5,7.4), labels=round(expm1(c(4.8,5.7,6.5,7.4))))
pdf(file = "Fig_comparison_OK_RF_zinc_meuse.pdf", width=9, height=9)
par(mfrow=c(2,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_ok"])), col=leg, zlim=c(4.8,7.4), main="Ordinary Kriging (OK)", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd"])), col=leg, zlim=c(4.8,7.4), main="Random Forest (RF)", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(sqrt(raster(meuse.grid["zinc_ok_var"])), col=rev(bpy.colors()), zlim=c(0,var.max), main="OK prediction error", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfd_var"]), col=rev(bpy.colors()), zlim=c(0,var.max), main="RF prediction error", axes=FALSE, box=FALSE)
points(meuse, pch="+")
dev.off()
## TH: RF smooths somewhat more than OK

## cross-validation:

## RF with combined covariates ----
meuse.grid$SW_occurrence = readGDAL("Meuse_GlobalSurfaceWater_occurrence.tif")$band1[meuse.grid@grid.index]
meuse.grid$AHN = readGDAL("ahn.asc")$band1[meuse.grid@grid.index]
meuse.grid$LGN5 = as.factor(readGDAL("lgn5.asc")$band1[meuse.grid@grid.index])
grids.spc = spc(meuse.grid, as.formula("~ SW_occurrence + AHN + ffreq + dist"))
## fit hybrid RF model:
fm1 <- as.formula(paste("zinc ~ ", dn0, " + ", paste(names(grids.spc@predicted), collapse = "+")))
ov.zinc1 <- over(meuse["zinc"], grids.spc@predicted)
m1.zinc <- ranger(fm1, do.call(cbind, list(meuse@data["zinc"], ov.zinc, ov.zinc1)), keep.inbag = TRUE, importance = "impurity")
m1.zinc
zinc.rfd1 <- predict(m1.zinc, cbind(grid.dist0@data, grids.spc@predicted@data), type = "se")
meuse.grid$zinc_rfd1 = zinc.rfd1$predictions
meuse.grid$zinc_rfd1_var = zinc.rfd1$se
xl <- as.list(ranger::importance(m1.zinc))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
m2.zinc <- ranger(paste("zinc ~ ", paste(names(grids.spc@predicted), collapse = "+")), do.call(cbind, list(meuse@data["zinc"], ov.zinc1)))
m2.zinc
meuse.grid$zinc_rfd2 = predict(m2.zinc, grids.spc@predicted@data)$predictions

pdf(file = "Fig_RF_covs_bufferdist_zinc_meuse.pdf", width=9, height=5)
par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_rfd2"])), col=leg, zlim=c(4.8,7.4), main="Random Forest (RF) covs only", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd1"])), col=leg, zlim=c(4.8,7.4), main="Random Forest (RF) covs + buffer dist.", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
dev.off()

## SIC 1997 data set ----
## measurements made in Switzerland on the 8th of May 1986
sic97.sp = readRDS("sic97.rds")
swiss1km = readRDS("swiss1km.rds")
ov = over(y=swiss1km, x=sic97.sp)
#plot(stack(swiss1km[1:2]))
## linear geostatistical model from: https://www.jstatsoft.org/article/view/v063i12/v63i12.pdf 
swissFit <- lgm(rainfall ~ CHELSA_rainfall + DEM, sic97.sp, grid = raster(swiss1km["border"]), covariates = stack(swiss1km[c("CHELSA_rainfall","DEM")]), aniso = TRUE)
## Example from the paper not working for some reason?!
#swissFit <- lgm(rainfall ~ CHELSA_rainfall + DEM, sic97.sp, grid = raster(swiss1km["border"]), covariates = stack(swiss1km[c("CHELSA_rainfall","DEM")]),, shape = 1, fixShape = TRUE, boxcox = 0.5, fixBoxcox = TRUE, aniso = TRUE)
# Error in (function (formula, data, paramToEstimate = c("range", "nugget"),  : 
# L-BFGS-B needs finite values of 'fn'
swiss1km$rainfall_UK = as(swissFit$predict[["predict"]], "SpatialGridDataFrame")@data[swiss1km@grid.index,1]
swiss1km$rainfall_krigeSd = as(swissFit$predict[["krigeSd"]], "SpatialGridDataFrame")@data[swiss1km@grid.index,1]

## Random Forest example:
swiss.dist0 <- buffer.dist(sic97.sp["rainfall"], swiss1km[1], as.factor(1:nrow(sic97.sp))) ## takes 2-3 mins
ov.swiss = over(sic97.sp["rainfall"], swiss.dist0)
sw.dn0 <- paste(names(swiss.dist0), collapse="+")
sw.fm1 <- as.formula(paste("rainfall ~ ", sw.dn0, " + CHELSA_rainfall + DEM"))
ov.rain <- over(sic97.sp["rainfall"], swiss1km[1:2])
sw.rm = do.call(cbind, list(sic97.sp@data["rainfall"], ov.rain, ov.swiss))
m1.rain <- ranger(sw.fm1, sw.rm[complete.cases(sw.rm),], keep.inbag = TRUE, importance = "impurity")
m1.rain
rain.rfd1 <- predict(m1.rain, cbind(swiss.dist0@data, swiss1km@data), type = "se")
swiss1km$rainfall_rfd1 = rain.rfd1$predictions
swiss1km$rainfall_rfd1_var = rain.rfd1$se
xl1 <- as.list(ranger::importance(m1.rain))
print(t(data.frame(xl1[order(unlist(xl1), decreasing=TRUE)[1:15]])))

rain.max = max(swiss1km$rainfall_rfd1, na.rm = TRUE)
swiss1km$rainfall_UK = ifelse(swiss1km$rainfall_UK<0, 0, ifelse(swiss1km$rainfall_UK>rain.max, rain.max, swiss1km$rainfall_UK))
## Plot predictions next to each other:
pdf(file = "Fig_Swiss_rainfall_UK_vs_RF.pdf", width=12, height=8)
par(mfrow=c(2,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
plot(raster(swiss1km["rainfall_UK"]), col=leg, main="Universal kriging (UK)", axes=FALSE, box=FALSE, zlim=c(0, rain.max))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1"]), col=leg, main="Random Forest (RF)", axes=FALSE, box=FALSE, zlim=c(0, rain.max))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_krigeSd"]), col=rev(bpy.colors()), main="Universal kriging (UK) prediction error", axes=FALSE, box=FALSE, zlim=c(0,140))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1_var"]), col=rev(bpy.colors()), main="Random Forest (RF) prediction error", axes=FALSE, box=FALSE, zlim=c(0,140))
points(sic97.sp, pch="+")
dev.off()

nl.rd <- getURL("http://spatialreference.org/ref/sr-org/6781/proj4/")
## Geul data ----
geul <- read.table("geul.dat", header = TRUE, as.is = TRUE)
geul$pb = as.numeric(geul$pb)
geul = geul[!is.na(geul$pb),]
coordinates(geul) <- ~x+y
proj4string(geul) <- CRS(nl.rd) 
grd25 <- readGDAL("dem25.txt")
grd25 <- as(grd25, "SpatialPixelsDataFrame")
proj4string(grd25) = proj4string(geul) 

## Pb predicted using OK ----
pb.geo <- as.geodata(geul["pb"])
pb.vgm <- likfit(pb.geo, lambda=0, messages=FALSE, ini=c(var(log1p(pb.geo$data)),500), cov.model="exponential")
locs2 = grd25@coords
pb.ok <- krige.conv(pb.geo, locations=locs2, krige=krige.control(obj.m=pb.vgm))
grd25$pb_ok = pb.ok$predict
## Pb predicted using RF only ----
ov.geul = over(geul["pb"], grd25)
summary(ov.geul$band1)
geul.s = geul[!is.na(ov.geul$band1),"pb"]
grid.dist1 <- buffer.dist(geul.s, grd25[1], as.factor(1:nrow(geul.s)))
dn1 <- paste(names(grid.dist1), collapse="+")
fm1 <- as.formula(paste("pb ~", dn1))
m1 <- fit.gstatModel(geul.s, fm1, grid.dist1, method="ranger", rvgm=NULL)
rk.m1 <- predict(m1, grid.dist1)

## Plot predictions next to each other:
grd25$pb_ok = ifelse(grd25$pb_ok<expm1(4.2), expm1(4.2), grd25$pb_ok)
png(file = "Fig_comparison_OK_RF_Pb_Geul.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(grd25["pb_ok"])), col=leg, , zlim=c(4.2,6.6), main="geoR (krige.conv)")
points(geul.s, pch="+")
plot(log1p(raster(rk.m1@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest")
points(geul.s, pch="+")
dev.off()

## RF with both buffer dist and covariates ----
grd25$swi <- readGDAL("swi.sdat")$band1[grd25@grid.index]
grd25$dis <- readGDAL("riverdist.txt")$band1[grd25@grid.index]
plot(stack(grd25))
grd25T <- grd25[c("band1","swi","dis")]
grd25T@data <- cbind(grd25T@data, grid.dist1@data)
## Run principal component analysis:
grd25.spc <- spc(grd25T, as.formula(paste("~", paste(names(grd25T), collapse = "+"))))
plot(stack(grd25.spc@predicted[1:6]))
fm2 <- as.formula(paste("pb ~", paste(names(grd25.spc@predicted), collapse = "+")))
m2 <- fit.gstatModel(geul.s, fm2, grd25.spc@predicted, method="quantregForest", rvgm=NULL)
plot(m2)
dev.off()
rk.m2 <- predict(m2, grd25.spc@predicted)
#plot(rk.m2, col=leg)
varImpPlot(m1@regModel)

rk.m2@predicted$pb = ifelse(rk.m2@predicted$pb<expm1(4.2), expm1(4.2), rk.m2@predicted$pb)
png(file = "Fig_comparison_RF_covariates_Pb_Geul.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(rk.m2@predicted[2])), col=leg, , zlim=c(4.2,6.6), main="Random Forest + covs")
points(geul.s, pch="+")
plot(log1p(raster(rk.m1@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest")
points(geul.s, pch="+")
dev.off()

## Factor type variable (Ebergotzen data set) ----
data(eberg)
eberg <- eberg[runif(nrow(eberg))<.3,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
eberg = eberg[!is.na(eberg$soiltype),]
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
## predict soil types:
soiltype <- autopredict(eberg["soiltype"], eberg_grid, auto.plot=FALSE)
plot(stack(soiltype$predicted), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,100))
## Plot soil type "G"
r.G = soiltype$predicted["G"]
r.G$G = ifelse(r.G$G>40, 40, r.G$G)
plot(raster(r.G), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40))
points(eberg[eberg$soiltype=="G",], pch=19)
points(eberg[!eberg$soiltype=="G",], pch="+")
## Plot soil type "D"
r.D = soiltype$predicted["D"]
r.D$D = ifelse(r.D$D>40, 40, r.D$D)
plot(raster(r.D), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40))
points(eberg[eberg$soiltype=="D",], pch=19)
points(eberg[!eberg$soiltype=="D",], pch="+")
## Conclusion: looks like regression-kriging on class probs