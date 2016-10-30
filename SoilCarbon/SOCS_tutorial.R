## Derivation of Soil Organic Carbon Stocks (SOCS)
## tom.hengl@isric.org

list.of.packages = c("GSIF", "plotKML", "nnet", "plyr", "ROCR", "randomForest", "plyr", "parallel", "psych", "mda", "h2o", "dismo", "grDevices", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret", "mboost")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#library(geospt)
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

## 1. SOCS points from La Libertad Research Center (Colombia)
## https://github.com/cran/geospt/tree/master/data
load("COSha10.rda")
load("COSha30.rda")
load("COSha30map.rda")
proj4string(COSha30map) = "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
plotKML(COSha30map, color=var1.pred, colour_scale=SAGA_pal[[1]])
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
proj4string(covs30m) = proj4string(COSha30)
covs30m.spc = spc(covs30m, fm.spc)
ov.COSha30 = cbind(as.data.frame(COSha30), over(COSha30, covs30m.spc@predicted))
#plot(stack(covs30m.spc@predicted))

## Model building:
fm.COSha30 = as.formula(paste("COSha30 ~ ", paste(names(covs30m.spc@predicted), collapse = "+")))
fitControl <- trainControl(method="repeatedcv", number=3, repeats=3)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
mFit1 <- train(fm.COSha30, data=ov.COSha30, method="ranger", trControl=fitControl)
mFit1
mFit2 <- train(fm.COSha30, data=ov.COSha30, method="xgbTree", trControl=fitControl, tuneGrid=gb.tuneGrid)
mFit3 <- train(fm.COSha30, data=ov.COSha30, method="gamboost", trControl=fitControl)
w1 = max(mFit1$results$Rsquared) #26%
w2 = max(mFit2$results$Rsquared) #18%
w3 = max(mFit3$results$Rsquared) #14%
COSha30.pr = covs30m[1]
names(COSha30.pr) = "SOCS30cm_pred"
COSha30.pr@data[,1] = rowSums(cbind(w1/(w1+w2)*predict(mFit1, covs30m.spc@predicted@data), w2/(w1+w2)*predict(mFit2, covs30m.spc@predicted@data)))
plotKML(COSha30.pr, colour_scale=SAGA_pal[[1]])
writeGDAL(COSha30.pr, "SOCS30cm_pred.tif", options="COMPRESS=DEFLATE")
## Average and total SOCS:
mean(COSha30.pr$SOCS30cm_pred); mean(COSha30$COSha30, na.rm=TRUE)
## 50,1 tonnes / ha
sum(COSha30.pr$SOCS30cm_pred*30^2/1e4)
## 105,346 tonnes

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

## 2. SOCS using legacy soil profiles from Ethiopia
## http://www.isric.org/data/cascape-ethiopia-woredas-soil-landscape-resources
tig_profs = readOGR("eth_profs_tigray.shp", "eth_profs_tigray")
tig_admin = readOGR("eth_tigray.shp", "eth_tigray")
et.prj = "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
tig_profs = spTransform(tig_profs, CRS(et.prj))
tig_admin = spTransform(tig_admin, CRS(et.prj))

## Covariates at 250 m resolution
if(!file.exists("SRTMDEM250m.sdat")){
  system(paste(gdalwarp, ' /home/tom/Downloads/EarthEnv-DEM90_N10E035/EarthEnv-DEM90_N10E035.bil SRTMDEM250m.sdat -of \"SAGA\" -t_srs \"', proj4string(tig_admin), '\" -tr 250 250 -te ', paste(as.vector(tig_admin@bbox), collapse=" ")))
}
plot(raster("SRTMDEM250m.sdat"), col=SAGA_pal[[1]])
## DEM derivatives:
saga_DEM_derivatives("SRTMDEM250m.sgrd")
plot(raster("SRTMDEM250m_vbf.sdat"), col=SAGA_pal[[1]])
## Landsat images:
for(i in c("B3","B4","B5","B6","B7","B10")){
  if(!file.exists(paste0('L8_250m_',i,'.tif'))){
    system(paste0(gdalwarp, ' /home/tom/Downloads/LC81680512014322LGN00/LC81680512014322LGN00_',i,'.TIF L8_250m_',i,'.tif -co \"COMPRESS=DEFLATE\" -t_srs \"', proj4string(tig_admin), '\" -tr 250 250 -r \"average\" -te ', paste(as.vector(tig_admin@bbox), collapse=" ")))
  }
}
## Land cover from 2010:
if(!file.exists("LC_250m.tif")){
  system(paste(gdalwarp, ' /home/tom/Envirometrix/LDN/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif LC_250m.tif -t_srs \"', proj4string(tig_admin), '\" -co \"COMPRESS=DEFLATE\" -r \"near\" -tr 250 250 -te ', paste(as.vector(tig_admin@bbox), collapse=" ")))
}

## SoilGrids ORC and BLD maps:
for(i in 1:4){
  if(!file.exists(paste0('ORCDRC_sl',i,'.tif'))){
    system(paste0(gdalwarp, ' /home/tom/Downloads/geonode-orcdrc_m_sl',i,'_250m.tif ORCDRC_sl',i,'.tif -co \"COMPRESS=DEFLATE\" -t_srs \"', proj4string(tig_admin), '\" -tr 250 250 -te ', paste(as.vector(tig_admin@bbox), collapse=" ")))
  }
  if(!file.exists(paste0('BLDFIE_sl',i,'.tif'))){
    system(paste0(gdalwarp, ' /home/tom/Downloads/geonode-bldfie_m_sl',i,'_250m.tif BLDFIE_sl',i,'.tif -co \"COMPRESS=DEFLATE\" -t_srs \"', proj4string(tig_admin), '\" -tr 250 250 -te ', paste(as.vector(tig_admin@bbox), collapse=" ")))
  }
}

## Study area of interest:
et_mask = rasterize(tig_admin, raster("LC_250m.tif"))

## Overlay points / covariates:
covs250m = stack(c(list.files(pattern=glob2rx("SRTMDEM250m*.sdat$")),list.files(pattern=glob2rx("L8_250m_*.tif$"))))
covs250m = as(as(covs250m, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
plot(stack(covs250m), col=SAGA_pal[[1]])
classes = cut(COSha30$COSha30, breaks=seq(0, 100, length=10))
covs250mdist = buffer.dist(COSha30["COSha30"], covs250m[1], classes)
covs250m@data = cbind(covs250m@data, covs250mdist@data)
fm.spc = as.formula(paste(" ~ ", paste(names(covs250m), collapse = "+")))
proj4string(covs250m) = proj4string(COSha30)
covs250m.spc = spc(covs250m, fm.spc)
ov.COSha30 = cbind(as.data.frame(COSha30), over(COSha30, covs250m.spc@predicted))

## Model building:
fm.COSha30 = as.formula(paste("COSha30 ~ ", paste(names(covs250m.spc@predicted), collapse = "+")))
fitControl <- trainControl(method="repeatedcv", number=3, repeats=3)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
mFit1 <- train(fm.COSha30, data=ov.COSha30, method="ranger", trControl=fitControl)
mFit1
mFit2 <- train(fm.COSha30, data=ov.COSha30, method="xgbTree", trControl=fitControl, tuneGrid=gb.tuneGrid)
mFit3 <- train(fm.COSha30, data=ov.COSha30, method="gamboost", trControl=fitControl)
w1 = max(mFit1$results$Rsquared) #26%
w2 = max(mFit2$results$Rsquared) #18%
w3 = max(mFit3$results$Rsquared) #14%
COSha30.pr = covs250m[1]
names(COSha30.pr) = "SOCS30cm_pred"
COSha30.pr@data[,1] = rowSums(cbind(w1/(w1+w2)*predict(mFit1, covs250m.spc@predicted@data), w2/(w1+w2)*predict(mFit2, covs250m.spc@predicted@data)))
plotKML(COSha30.pr, colour_scale=SAGA_pal[[1]])
writeGDAL(COSha30.pr, "SOCS30cm_pred.tif", options="COMPRESS=DEFLATE")
## Average and total SOCS:
mean(COSha30.pr$SOCS30cm_pred); mean(COSha30$COSha30, na.rm=TRUE)
## 50,1 tonnes / ha
sum(COSha30.pr$SOCS30cm_pred*30^2/1e4)
## 105,346 tonnes