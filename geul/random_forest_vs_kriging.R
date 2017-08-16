## Comparison RF vs kriging
## tom.hengl@isric.org

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
#library(geoR)
#.onLoad failed in loadNamespace() for 'tcltk', details
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")
## Load the Meuse data set:
demo(meuse, echo=FALSE)

## compare GLM vs RF
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

## Zinc predicted using ordinary kriging:
zinc.geo <- as.geodata(meuse["zinc"])
#plot(variog4(zinc.geo, lambda=0, max.dist=1500, messages=FALSE), lwd=2)
zinc.vgm <- likfit(zinc.geo, lambda=0, messages=FALSE, ini=c(var(log1p(zinc.geo$data)),500), cov.model="exponential")
locs = meuse.grid@coords
zinc.ok <- krige.conv(zinc.geo, locations=locs, krige=krige.control(obj.m=zinc.vgm))
meuse.grid$zinc_ok = zinc.ok$predict

## Zinc predicted using RF and buffer distances only:
grid.dist0 <- buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
dn0 <- paste(names(grid.dist0), collapse="+")
fm0 <- as.formula(paste("zinc ~", dn0))
m0 <- fit.gstatModel(meuse, fm0, grid.dist0, method="ranger", rvgm=NULL)
rk.m0 <- predict(m0, grid.dist0)
plot(rk.m0)

## Plot predictions next to each other:
png(file = "Fig_comparison_OK_RF_zinc_meuse.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(meuse.grid["zinc_ok"])), col=leg, zlim=c(4.8,7.4), main="geoR (krige.conv)")
points(meuse, pch="+")
plot(log1p(raster(rk.m0@predicted[2])), col=leg, zlim=c(4.8,7.4), main="Random Forest")
points(meuse, pch="+")
dev.off()
## TH: RF smooths somewhat more than OK (TO-DO: test bias and accuracy using cross-validation)

nl.rd <- getURL("http://spatialreference.org/ref/sr-org/6781/proj4/")
## Geul data:
geul <- read.table("geul.dat", header = TRUE, as.is = TRUE)
geul$pb = as.numeric(geul$pb)
geul = geul[!is.na(geul$pb),]
coordinates(geul) <- ~x+y
proj4string(geul) <- CRS(nl.rd) 
grd25 <- readGDAL("dem25.txt")
grd25 <- as(grd25, "SpatialPixelsDataFrame")
proj4string(grd25) = proj4string(geul) 

## Pb predicted using OK:
pb.geo <- as.geodata(geul["pb"])
pb.vgm <- likfit(pb.geo, lambda=0, messages=FALSE, ini=c(var(log1p(pb.geo$data)),500), cov.model="exponential")
locs2 = grd25@coords
pb.ok <- krige.conv(pb.geo, locations=locs2, krige=krige.control(obj.m=pb.vgm))
grd25$pb_ok = pb.ok$predict
## RF
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

## RF with both buffer dist and covariates:
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

## Ebergotzen data set:
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