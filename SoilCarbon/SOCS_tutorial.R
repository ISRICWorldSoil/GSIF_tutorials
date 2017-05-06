## Spatial prediction and derivation of Soil Organic Carbon Stock (OCS)
## prepared by: tom.hengl@isric.org
## http://gsif.isric.org/doku.php/wiki:soil_organic_carbon

list.of.packages = c("GSIF", "plotKML", "nnet", "aqp", "plyr", "ROCR", "randomForest", "parallel", "psych", "mda", "dismo", "grDevices", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret", "mboost")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#library(geospt)
library(aqp)
library(plyr)
library(sp)
library(scales)
library(rgdal)
library(raster)
library(caret)
library(ranger)
library(randomForest)
library(xgboost)
library(plotKML)
library(GSIF)
## GDAL (https://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries):
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate = "gdal_translate"
  gdalwarp = "gdalwarp"
}
load(".RData")
source('SOCS_functions.R')

## 0. Fitting splines for profiles ----

## http://www.asris.csiro.au/mapping/hyperdocs/NatSoil/399%5EEDGEROI%5Eed079.pdf
lon = 149.73; lat = -30.09; id = "399_EDGEROI_ed079"; TIMESTRR = "1987-01-05"
top = c(0, 10, 20, 55, 90) 
bottom = c(10, 20, 55, 90, 116)
ORC = c(8.2, 7.5, 6.1, 3.3, 1.6)
BLD = c(1340, 1367, 1382, 1433, 1465)
CRF = c(6, 6, 7, 8, 8)
#OCS = OCSKGM(ORC, BLD, CRF, HSIZE=bottom-top)
prof1 <- join(data.frame(id, top, bottom, ORC, BLD, CRF), 
              data.frame(id, lon, lat, TIMESTRR), type='inner')
depths(prof1) <- id ~ top + bottom
site(prof1) <- ~ lon + lat + TIMESTRR
coordinates(prof1) <- ~ lon + lat
proj4string(prof1) <- CRS("+proj=longlat +datum=WGS84")
## fit a spline:
ORC.s <- mpspline(prof1, var.name="ORC", d=t(c(0,30,100,200)), vhigh = 2200)
BLD.s <- mpspline(prof1, var.name="BLD", d=t(c(0,30,100,200)), vhigh = 2200)
CRF.s <- mpspline(prof1, var.name="CRF", d=t(c(0,30,100,200)), vhigh = 2200)
## derive OCS for standard depths:
OCSKGM(ORC.s$var.std$`0-30 cm`, BLD.s$var.std$`0-30 cm`, CRF.s$var.std$`0-30 cm`, HSIZE=30)
OCSKGM(ORC.s$var.std$`30-100 cm`, BLD.s$var.std$`30-100 cm`, CRF.s$var.std$`30-100 cm`, HSIZE=70)

## Organic soil from Canada:
lon = -97.08639; lat = 51.12972; id = "CanSIS:3010"; TIMESTRR = "1966-06-05"
top = c(0, 31, 61, 91, 122) 
bottom = c(31, 61, 91, 122, 130)
ORC = c(472, 492, 487, 502, 59)
BLD = c(179, 166, 169, 160, 787)
CRF = c(5, 6, 6, 6, 6)
#OCS = OCSKGM(ORC, BLD, CRF, HSIZE=bottom-top)
prof2 <- join(data.frame(id, top, bottom, ORC, BLD, CRF), 
              data.frame(id, lon, lat, TIMESTRR), type='inner')
depths(prof2) <- id ~ top + bottom
site(prof2) <- ~ lon + lat + TIMESTRR
coordinates(prof2) <- ~ lon + lat
proj4string(prof2) <- CRS("+proj=longlat +datum=WGS84")
## fit a spline:
ORC.s <- mpspline(prof2, var.name="ORC", d=t(c(0,30,100,200)), vhigh = 2200)
BLD.s <- mpspline(prof2, var.name="BLD", d=t(c(0,30,100,200)), vhigh = 2200)
CRF.s <- mpspline(prof2, var.name="CRF", d=t(c(0,30,100,200)), vhigh = 2200)
## derive OCS for standard depths:
OCSKGM(ORC.s$var.std$`0-30 cm`, BLD.s$var.std$`0-30 cm`, CRF.s$var.std$`0-30 cm`, HSIZE=30)
OCSKGM(ORC.s$var.std$`30-100 cm`, BLD.s$var.std$`30-100 cm`, CRF.s$var.std$`30-100 cm`, HSIZE=70)

## 1. Model performance ----
fitControl <- trainControl(method="repeatedcv", number=2, repeats=2)
demo(meuse, echo=FALSE)
meuse.ov <- cbind(over(meuse, meuse.grid), meuse@data)
mFit0 <- train(om~1, data=meuse.ov, method="glm", family=gaussian(link=log), trControl=fitControl, na.action=na.omit)
mFit1 <- train(om~soil, data=meuse.ov, method="glm", family=gaussian(link=log), trControl=fitControl, na.action=na.omit)
mFit2 <- train(om~dist+soil+ffreq, data=meuse.ov, method="glm", family=gaussian(link=log), trControl=fitControl, na.action=na.omit)
mFit3 <- train(om~dist+soil+ffreq, data=meuse.ov, method="ranger", trControl=fitControl, na.action=na.omit)
## compare performance:
resamps <- resamples(list(Mean=mFit0, Soilmap=mFit1, GLM=mFit2, RF=mFit3))
summary(resamps)
bwplot(resamps, layout = c(3, 1))
## Improvement in efficiency = 32%:
round((1-min(mFit3$results$RMSE)/min(mFit0$results$RMSE))*100)

## 2. SOCS points from La Libertad Research Center (Colombia) ----
## https://github.com/cran/geospt/tree/master/data
load("COSha10.rda")
load("COSha30.rda")
str(COSha30)
load("COSha30map.rda")
proj4string(COSha30map) = "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
plotKML(COSha30map, color=var1.pred, colour_scale=SAGA_pal[[1]], png.type = "cairo")
coordinates(COSha30) = ~ x+y
proj4string(COSha30) = proj4string(COSha30map)
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(COSha30, colour=COSha30, shape=shape, labels=COSha30$COSha30, colour_scale=SAGA_pal[[1]])
writeOGR(COSha30, "COSha30.shp", "COSha30", "ESRI Shapefile")

## Download extra covariates:
#system(paste(gdalwarp, ' /home/tom/Downloads/N004W074/AVERAGE/N004W074_AVE_DSM.tif DEM30m.sdat -of \"SAGA\" -t_srs \"', proj4string(COSha30map), '\" -tr 30 30 -te ', paste(as.vector(COSha30map@bbox), collapse=" ")))
if(!file.exists("SRTMDEM30m.sdat")){
  system(paste(gdalwarp, ' /home/tom/Downloads/n04_w074_1arc_v3.tif SRTMDEM30m.sdat -of \"SAGA\" -t_srs \"', proj4string(COSha30map), '\" -tr 30 30 -te ', paste(as.vector(COSha30map@bbox), collapse=" ")))
}
plot(raster("SRTMDEM30m.sdat"), col=SAGA_pal[[1]])
## DEM derivatives:
saga_DEM_derivatives("SRTMDEM30m.sgrd")
plot(raster("SRTMDEM30m_vbf.sdat"), col=SAGA_pal[[1]])
## Landsat images:
for(i in c("B3","B4","B5","B6","B7","B10")){
  if(!file.exists(paste0('L8',i,'.tif'))){
    system(paste0(gdalwarp, ' /home/tom/Downloads/LC80070572014090LGN00/LC80070572014090LGN00_',i,'.TIF L8',i,'.tif -co \"COMPRESS=DEFLATE\" -t_srs \"', proj4string(COSha30map), '\" -tr 30 30 -te ', paste(as.vector(COSha30map@bbox), collapse=" ")))
  }
}

## Overlay points / covariates:
covs30m = stack(c(list.files(pattern=glob2rx("SRTMDEM30m*.sdat$")),list.files(pattern=glob2rx("L8B*.tif$"))))
covs30m = as(as(covs30m, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
plot(stack(covs30m), col=SAGA_pal[[1]])
classes = cut(COSha30$COSha30, breaks=seq(0, 100, length=10))
covs30mdist = buffer.dist(COSha30["COSha30"], covs30m[1], classes)
#plot(stack(covs30mdist))
covs30m@data = cbind(covs30m@data, covs30mdist@data)
fm.spc = as.formula(paste(" ~ ", paste(names(covs30m), collapse = "+")))
fm.spc
proj4string(covs30m) = proj4string(COSha30)
covs30m.spc = spc(covs30m, fm.spc)
ov.COSha30 = cbind(as.data.frame(COSha30), over(COSha30, covs30m.spc@predicted))
#plot(stack(covs30m.spc@predicted))

## Model building:
fm.COSha30 = as.formula(paste("COSha30 ~ ", paste(names(covs30m.spc@predicted), collapse = "+")))
fm.COSha30
fitControl <- trainControl(method="repeatedcv", number=3, repeats=3)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
mFit1 <- train(fm.COSha30, data=ov.COSha30, method="ranger", trControl=fitControl)
mFit1
mFit2 <- train(fm.COSha30, data=ov.COSha30, method="xgbTree", trControl=fitControl, tuneGrid=gb.tuneGrid)
mFit3 <- train(fm.COSha30, data=ov.COSha30, method="gamboost", trControl=fitControl)
w1 = max(mFit1$results$Rsquared) #24%
w2 = max(mFit2$results$Rsquared) #25%
w3 = max(mFit3$results$Rsquared) #17%
COSha30.pr = covs30m[1]
names(COSha30.pr) = "SOCS30cm_pred"
COSha30.pr@data[,1] = rowSums(cbind(w1/(w1+w2)*predict(mFit1, covs30m.spc@predicted@data), w2/(w1+w2)*predict(mFit2, covs30m.spc@predicted@data)))
plot(raster(COSha30.pr), col=SAGA_pal[[1]])
points(COSha30, pch="+")
plotKML(COSha30.pr, colour_scale=SAGA_pal[[1]], png.type = "cairo")
writeGDAL(COSha30.pr, "SOCS30cm_pred.tif", options="COMPRESS=DEFLATE")
## Average and total SOCS:
mean(COSha30.pr$SOCS30cm_pred); mean(COSha30$COSha30, na.rm=TRUE)
## 50,1 tonnes / ha in average
sum(COSha30.pr$SOCS30cm_pred*30^2/1e4)
## 106,809 tonnes

## Plot fitting success:
yr = range(COSha30.pr$SOCS30cm_pred, na.rm=TRUE)
panel_line = function(x,y,...) {
  panel.xyplot(x, y)
  panel.abline(0,1, ...)
}
xyplot(mFit1$finalModel$predictions~COSha30$COSha30, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE, cex=1.8), y=list(log=TRUE, equispaced.log=FALSE, cex=1.8)), xlab=list(label="measured (t / ha)", cex=1.8), ylab=list(label="predicted with machine learning (t / ha)", cex=1.8), xlim=yr, ylim=yr, panel=panel_line)
## Uncertainty stripped off:
yrl = range(log1p(COSha30$COSha30), na.rm=TRUE)
## RMSE in log space:
byl = sqrt((1-w1)*var(log1p(COSha30$COSha30), na.rm=TRUE))/2
brks1 = seq(yrl[1], yrl[2], by=byl)
labl = round(expm1(rowMeans(cbind(brks1[-length(brks1)], brks1[-1]))))
COSha30.pr$SOCS30cm_predC = cut(log1p(COSha30.pr$SOCS30cm_pred), breaks=brks1, labels=labl)
spplot(COSha30.pr["SOCS30cm_predC"], col.regions=SAGA_pal[[1]], sp.layout=list(list("sp.points", COSha30, pch="+", col="black")))

