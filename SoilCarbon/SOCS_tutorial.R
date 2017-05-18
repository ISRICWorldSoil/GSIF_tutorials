## Spatial prediction and derivation of Soil Organic Carbon Stock (OCS)
## prepared by: tom.hengl@isric.org
## http://gsif.isric.org/doku.php/wiki:soil_organic_carbon

list.of.packages = c("GSIF", "plotKML", "nnet", "aqp", "plyr", "ROCR", "randomForest", "parallel", "psych", "mda", "dismo", "grDevices", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret", "mboost")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

load(".RData")
#library(geospt)
library(aqp)
library(plyr)
library(sp)
library(scales)
library(rgdal)
library(raster)
library(caret)
library(utils)
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
source('SOCS_functions.R')
md = getwd()

## 0. Fitting splines for profiles ----

## Mineral soil profile
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

## 2. Global PTF for estimating bulk density from profile data ----
if(!file.exists("wosis_tbl.rds")){
  download.file("http://gsif.isric.org/lib/exe/fetch.php?media=wosis_profiles.rda", "wosis_profiles.rda")
  load("wosis_profiles.rda"); unlink("wosis_profiles.rda")
  ## Get USDA taxa:
  wosis_pnts = wosis_profiles$sites[,c("profile_id","longitude","latitude")]
  coordinates(wosis_pnts) = ~ longitude + latitude
  proj4string(wosis_pnts) = "+proj=longlat +datum=WGS84"
  ov.tax = raster::extract(raster("/mnt/cartman/ftp.soilgrids.org/data/recent/TAXOUSDA_250m_ll.tif"), wosis_pnts)
  ## Bulk density fine earth:
  summary(wosis_profiles$horizons$bd_fe)
  wosis_profiles$horizons$bd_fe = ifelse(wosis_profiles$horizons$bd_fe < .085, NA, wosis_profiles$horizons$bd_fe)
  wosis_tbl = plyr::join(wosis_profiles$horizons[,c("profile_id", "upper_depth", "lower_depth", "bd_fe", "carbon_o", "ph_h2o", "sand_t", "clay_t")], data.frame(profile_id=wosis_pnts$profile_id, TAXOUSDA=as.factor(ov.tax)))
  wosis_tbl = wosis_tbl[complete.cases(wosis_tbl),]
  #'data.frame':	101570 obs. of  9 variables:
  summary(wosis_tbl$TAXOUSDA)
  wosis_tbl$depth = (wosis_tbl$lower_depth - wosis_tbl$upper_depth)/2 + wosis_tbl$upper_depth
  wosis_tbl$bd = wosis_tbl$bd_fe * 1000
  ind.tax = data.frame(model.matrix(~TAXOUSDA-1, wosis_tbl))
  wosis_tbl = cbind(wosis_tbl, ind.tax)
  #str(wosis_tbl)  
  saveRDS(wosis_tbl, "wosis_tbl.rds")
} else {
  wosis_tbl = readRDS("wosis_tbl.rds")
}
library(ranger)
hist(log1p(wosis_tbl$carbon_o), breaks=45)
fm = as.formula(paste0('bd ~ carbon_o + ph_h2o + sand_t + clay_t + depth +', paste(names(ind.tax), collapse="+")))
rf_BD <- ranger(fm, data=wosis_tbl, num.trees = 85, importance='impurity')
rf_BD
library(scales)
library(hexbin)
library(lattice)
library(plotKML)
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)
  panel.loess(x, y, ..., col = "black",lty=1,lw=2,span=1/18)
}
pal = R_pal[["bpy_colors"]][3:18]
hexbinplot(wosis_tbl$bd~wosis_tbl$carbon_o, colramp=colorRampPalette(pal), xlab="Organic carbon (permiles)", ylab="Bulk density (kg/cubic-m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, asp=1, xbins=40, ybins=40, xlim=c(0,400), panel=pfun, colorcut=c(0,0.008,0.01,0.018,0.028,0.07,0.15,0.25,0.5,0.75,1))

## 3. SOCS points from La Libertad Research Center (Colombia) ----
## https://github.com/cran/geospt/tree/master/data
setwd("./LRC")
load("COSha10.rda")
load("COSha30.rda")
str(COSha30)
load("COSha30map.rda")
proj4string(COSha30map) = "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#plotKML(COSha30map, color=var1.pred, colour_scale=SAGA_pal[[1]], png.type = "cairo")
coordinates(COSha30) = ~ x+y
proj4string(COSha30) = proj4string(COSha30map)
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(COSha30, colour=COSha30, shape=shape, labels=round(COSha30$COSha30), colour_scale=SAGA_pal[[1]])
writeOGR(COSha30, "COSha30.shp", "COSha30", "ESRI Shapefile")
## resample to 30 by 30 m:
writeGDAL(COSha30map["var1.pred"], "COSha30map_var1pred.tif", options=c("COMPRESS=DEFLATE"))
system('gdalwarp COSha30map_var1pred.tif COSha30map_var1pred_.tif -tr 30 30')
rLRC = raster("COSha30map_var1pred_.tif")
unlink("COSha30map_var1pred.tif")

## Prepare covariates:
if(!file.exists("covs30m.rds")){
  get_30m_covariates(te=paste(as.vector(extent(rLRC))[c(1,3,2,4)], collapse=" "), tr=res(rLRC)[1], p4s=proj4string(rLRC))
  saga_DEM_derivatives("SRTMGL1_SRTMGL1.2.tif")
  plot(raster("SRTMGL1_SRTMGL1.2_vbf.sdat"), col=SAGA_pal[[1]])
  unlink("SRTMGL1_SRTMGL1.2.tif")
  
  ## Overlay points / covariates:
  cov.lst = c(list.files(pattern=glob2rx("*.sdat$")),list.files(pattern=glob2rx("*.tif$")))
  covs30m = stack(cov.lst)
  covs30m = as(as(covs30m, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
  saveRDS(covs30m, "covs30m.rds")
} else {
  covs30m = readRDS("covs30m.rds")  
}

#plot(stack(covs30m), col=SAGA_pal[[1]])
#plot(raster(covs30m["Landsat_bare2010"]), col=SAGA_pal[[1]])
#plot(raster(covs30m["GlobalSurfaceWater_occurrence"]), col=SAGA_pal[[1]])
classes = cut(COSha30$COSha30, breaks=seq(30, 95, length=10))
covs30mdist = buffer.dist(COSha30["COSha30"], covs30m[1], classes)
#plot(stack(covs30mdist))
covs30m@data = cbind(covs30m@data, covs30mdist@data)
sel.rm = c("GlobalSurfaceWater_occurrence", "GlobalSurfaceWater_extent", "Landsat_bare2010", "COSha30map_var1pred_")
fm.spc = as.formula(paste(" ~ ", paste(names(covs30m)[-which(names(covs30m@data) %in% sel.rm)], collapse = "+")))
fm.spc
proj4string(covs30m) = proj4string(COSha30)
covs30m.spc = spc(covs30m, fm.spc)
ov.COSha30 = cbind(as.data.frame(COSha30), over(COSha30, covs30m.spc@predicted))
#plot(stack(covs30m.spc@predicted))

## Model building:
fm.COSha30 = as.formula(paste("COSha30 ~ ", paste(names(covs30m.spc@predicted), collapse = "+")))
fm.COSha30
fitControl <- trainControl(method="repeatedcv", number=3, repeats=2)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample=1)
mFit1 <- train(fm.COSha30, data=ov.COSha30, method="ranger", trControl=fitControl, importance='impurity')
mFit1
xl <- as.list(ranger::importance(mFit1$finalModel))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
#spplot(covs30m.spc@predicted["PC7"])
mFit2 <- train(fm.COSha30, data=ov.COSha30, method="xgbTree", trControl=fitControl, tuneGrid=gb.tuneGrid)
mFit2
mFit3 <- train(fm.COSha30, data=ov.COSha30, method="gamboost", trControl=fitControl)
mFit3
w1 = max(mFit1$results$Rsquared) #13%
w2 = max(mFit2$results$Rsquared) #13%
w3 = max(mFit3$results$Rsquared) #15%
## Two best models gamboost
COSha30.pr = covs30m["COSha30map_var1pred_"]
COSha30.pr@data[,"COSha30map_MLA"] = rowSums(cbind(w1/(w1+w3)*predict(mFit1, covs30m.spc@predicted@data), w3/(w1+w3)*predict(mFit2, covs30m.spc@predicted@data)))
COSha30.pr$COSha30map_MLA = ifelse(is.na(COSha30.pr$COSha30map_var1pred_), NA, COSha30.pr$COSha30map_MLA)
spplot(COSha30.pr, col.regions=SAGA_pal[[1]], sp.layout = list("sp.points", COSha30, pch = "+", col="black", cex=1.5))
#plotKML(COSha30.pr["COSha30map_MLA"], colour_scale=SAGA_pal[[1]], png.type = "cairo")
#writeGDAL(COSha30.pr["COSha30map_MLA"], "SOCS30cm_pred.tif", options="COMPRESS=DEFLATE")
## Average and total SOCS:
mean(COSha30.pr$COSha30map_MLA, na.rm=TRUE); mean(COSha30$COSha30, na.rm=TRUE)
## 51 tonnes / ha in average
sum(COSha30.pr$COSha30map_MLA*30^2/1e4, na.rm=TRUE)
## 48,608 tonnes in total

## Plot fitting success:
yr = range(COSha30.pr$COSha30map_MLA, na.rm=TRUE)
panel_line = function(x,y,...) {
  panel.xyplot(x, y)
  panel.abline(0,1, ...)
}
xyplot(mFit1$finalModel$predictions~COSha30$COSha30, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE, cex=1.8), y=list(log=TRUE, equispaced.log=FALSE, cex=1.8)), xlab=list(label="measured (t / ha)", cex=1.8), ylab=list(label="predicted with machine learning (t / ha)", cex=1.8), xlim=yr, ylim=yr, panel=panel_line)
## Uncertainty 'stripped-off':
yrl = range(log1p(COSha30$COSha30), na.rm=TRUE)
## RMSE in log space:
byl = sqrt((1-w1)*var(log1p(COSha30$COSha30), na.rm=TRUE))/2
brks1 = seq(yrl[1], yrl[2], by=byl)
labl = round(expm1(rowMeans(cbind(brks1[-length(brks1)], brks1[-1]))))
COSha30.pr$SOCS30cm_predC = cut(log1p(COSha30.pr$COSha30map_MLA), breaks=brks1, labels=labl)
spplot(COSha30.pr["SOCS30cm_predC"], col.regions=SAGA_pal[[1]], sp.layout=list(list("sp.points", COSha30, pch="+", col="black")))
setwd(md)

## 4. Soil profiles case study Bor (Serbia) ----
## https://github.com/pejovic/sparsereg3D
setwd("./bor")
load("BorData.rda")
coordinates(bor) = ~ x + y
proj4string(bor) = "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"
rBor = raster("bor_soiltype_50m.tif")
bor = spTransform(bor, CRS(proj4string(rBor)))
bor$depth = (bor$Bottom - bor$Top)/2+bor$Top

if(!file.exists("covs50m.rds")){
  get_30m_covariates(te=paste(as.vector(extent(rBor))[c(1,3,2,4)], collapse=" "), tr=res(rBor)[1], p4s=proj4string(rBor))
  saga_DEM_derivatives("SRTMGL1_SRTMGL1.2.tif")
  cov2.lst = c(list.files(pattern=glob2rx("*.sdat$")),list.files(pattern=glob2rx("*.tif$")))
  covs50m = stack(cov2.lst)
  covs50m = as(as(covs50m, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
  #plot(raster(covs50m["bor_soiltype_50m"]))
  #plot(raster(covs50m["GlobalForestChange2000.2014_first_REDL00"]))
  saveRDS(covs50m, "covs50m.rds")
} else {
  covs30m = readRDS("covs50m.rds")  
}

## Convert soil organic matter to OC:
bor$ORC = bor$SOM*0.58*10
summary(bor$ORC)
classes2 = cut(bor$ORC, breaks=seq(0, 250, length=15))
covs50mdist = buffer.dist(bor["ORC"], covs50m["bor_soiltype_50m"], classes2)
#plot(stack(covs50mdist))
covs50m@data = cbind(covs50m@data, covs50mdist@data)
covs50m$bor_soiltype_50m = as.factor(covs50m$bor_soiltype_50m)
sel.rm = c("GlobalSurfaceWater_occurrence", "GlobalSurfaceWater_extent", "Landsat_bare2010", "COSha30map_var1pred_")
fm.spc2 = as.formula(paste(" ~ ", paste(names(covs50m)[-which(names(covs50m@data) %in% sel.rm)], collapse = "+")))
covs50m.spc = spc(covs50m, fm.spc2)
#plot(stack(covs50m.spc@predicted[1:6]), col=SAGA_pal[[1]])

## Fill in the bulk density values using SoilGrids:
bor.ll = spTransform(bor, CRS("+proj=longlat +datum=WGS84"))
bor.ll = bor.ll[!duplicated(bor.ll$ID),c("ID","Soil.Type")]
## 206 unique points
## REST example:
library(rjson)
library(sp)
soilgrids.r <- REST.SoilGrids(c("BLDFIE"))
ov.sg <- over(soilgrids.r, bor.ll)
## Match values of BLD using standard depths:
## 7 standard depths
breaks = c(0, rowMeans(data.frame(c(0,5,15,30,60,100), c(5,15,30,60,100,200))), 200, 4500)
bor$depth_c = cut(bor$depth, breaks, labels = paste0("sl", 1:8))
summary(bor$depth_c)
#sl1 sl2 sl3 sl4 sl5 sl6 sl7 sl8 
#33  90 145 119  69  15   0   0
## Estimate Organic carbon density:
bor$BLDFIE = NA
bor$OCD = bor$ORC/1000 * bor$BLDFIE

ov.bor = cbind(as.data.frame(bor[c("ID","depth","OCD")]), over(bor, covs50m.spc@predicted))

## Model building:
fm.COSha30 = as.formula(paste("COSha30 ~ ", paste(names(covs30m.spc@predicted), collapse = "+")))
fm.COSha30
fitControl <- trainControl(method="repeatedcv", number=3, repeats=2)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample=1)
mFit1 <- train(fm.COSha30, data=ov.COSha30, method="ranger", trControl=fitControl, importance='impurity')
mFit1

## 5. Spatiotemporal modeling of OCD (USA48) ----
## 250,428 observations of organic carbon density (kg/m3) together with some 80+ covariates
OCD_stN = readRDS("usa48.OCD_spacetime_matrix.rds")
pr.lst <- names(OCD_stN)[-which(names(OCD_stN) %in% c("SOURCEID", "DEPTH.f", "OCDENS",   "YEAR", "YEAR_c", "LONWGS84", "LATWGS84"))]
## 132 vars
fm0.st <- as.formula(paste('OCDENS ~ DEPTH.f + ', paste(pr.lst, collapse="+")))
sel0.m = complete.cases(OCD_stN[,all.vars(fm0.st)])
rf0.OCD_st <- ranger(fm0.st, data=OCD_stN[sel0.m,all.vars(fm0.st)], importance="impurity", write.forest=TRUE, num.trees=120)
rf0.OCD_st
## RMSE = 10 kg/m3 / 60%
xl <- as.list(ranger::importance(rf0.OCD_st))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))

