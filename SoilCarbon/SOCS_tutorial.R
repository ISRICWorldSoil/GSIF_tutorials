## Spatial prediction and derivation of Soil Organic Carbon Stock (OCS)
## prepared by: tom.hengl@isric.org and bas.kempen@wur.nl
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
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")


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

## 1. Improving model performance ----
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

## 2. Fitting a global PTF for estimating BD from profile data ----
library(ranger)
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
  saveRDS(ind.tax, "ov_wosis_taxousda.rds")
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
} else {
  ind.tax = readRDS("ov_taxousda.rds")
  dfs_tbl = readRDS("wosis_tbl.rds")
  fm.BLD = as.formula(paste("BLD ~ ORCDRC + CLYPPT + SNDPPT + PHIHOX + DEPTH.f +", paste(names(ind.tax), collapse="+")))
  m.BLD_PTF <- ranger(fm.BLD, dfs_tbl, num.trees = 85, importance='impurity')
  ## takes 30--60 secs
  m.BLD_PTF
  saveRDS(m.BLD_PTF, "m.BLD_PTF.rds")
}

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

## 4. Soil profiles case study Edgeroi ----
setwd("./edgeroi")
data(edgeroi)
edgeroi.sp = edgeroi$sites
coordinates(edgeroi.sp) <- ~ LONGDA94 + LATGDA94
proj4string(edgeroi.sp) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
edgeroi.sp <- spTransform(edgeroi.sp, CRS("+init=epsg:28355"))
writeOGR(edgeroi.sp, "edgeroi.gpkg", "edgeroi", driver="GPKG")
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids.rda")
load(con)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
edgeroi.spc = spc(edgeroi.grids, ~DEMSRT5+TWISRT5+PMTGEO5+EV1MOD5+EV2MOD5+EV3MOD5)
#plot(stack(edgeroi.grids))
## Impute missing bulk density using SoilGrids250m:
summary(edgeroi$horizons[,c("ORCDRC","CLYPPT")])
ov.edgeroi.BLD = raster::extract(stack(paste0("/mnt/cartman/ftp.soilgrids.org/data/recent/BLDFIE_M_sl",1:7,"_250m_ll.tif")), spTransform(edgeroi.sp, CRS("+proj=longlat +datum=WGS84")))
## Derive averaged estimates for standard depth intervals:
ov.edgeroi.BLDm = data.frame(BLD.f=as.vector(sapply(2:ncol(ov.edgeroi.BLD), function(i){rowMeans(ov.edgeroi.BLD[,c(i-1,i)])})), DEPTH.c=as.vector(sapply(1:6, function(i){rep(paste0("sd",i), nrow(edgeroi$sites))})), SOURCEID=rep(edgeroi$sites$SOURCEID, 6))
str(ov.edgeroi.BLDm)
## Match BLD values with actual horizons:
edgeroi$horizons$DEPTH = edgeroi$horizons$UHDICM + (edgeroi$horizons$LHDICM - edgeroi$horizons$UHDICM)/2
edgeroi$horizons$DEPTH.c = cut(edgeroi$horizons$DEPTH, include.lowest=TRUE, breaks=c(0,5,15,30,60,100,1000), labels=paste0("sd",1:6))
summary(edgeroi$horizons$DEPTH.c)
edgeroi$horizons$BLD.f = plyr::join(edgeroi$horizons[,c("SOURCEID","DEPTH.c")], ov.edgeroi.BLDm)$BLD.f
## Derive Organic Carbon Density in x10kg/m3:
edgeroi$horizons$OCD = edgeroi$horizons$ORCDRC/1000 * edgeroi$horizons$BLD.f
summary(edgeroi$horizons$OCD)
## Prepare regression matrix:
ov2 <- over(edgeroi.sp, edgeroi.spc@predicted)
ov2$SOURCEID = edgeroi.sp$SOURCEID
h2 = hor2xyd(edgeroi$horizons)
## regression matrix:
m2 <- plyr::join_all(dfs = list(edgeroi$sites, h2, ov2))
## Fit a RF model:
fm.OCD = as.formula(paste0("OCD ~ DEPTH + ", paste(names(edgeroi.spc@predicted), collapse = "+")))
fm.OCD
m.OCD <- ranger(fm.OCD, m2[complete.cases(m2[,all.vars(fm.OCD)]),], keep.inbag = TRUE, importance = "impurity")
m.OCD
xl <- as.list(ranger::importance(m.OCD))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
## Predict OCD (with uncertainty):
for(i in c(0,30)){
  edgeroi.spc@predicted$DEPTH = i
  OCD.rf <- predict(m.OCD, edgeroi.spc@predicted@data)
  #system.time(OCD.rf <- predict(m.OCD, edgeroi.spc@predicted@data, type = "se"))
  edgeroi.grids@data[,paste0("OCD.", i, "cm")] = OCD.rf$predictions
  #edgeroi.grids@data[,paste0("OCD.", i, "cm_se")] = OCD.rf$se
}
## Deriving the error map is VERY computationally intensive (>5minutes), especially if the number of covariates is high
## Derive Organic carbon stocks in t/ha:
edgeroi.grids$OCS.30cm = rowMeans(edgeroi.grids@data[,paste0("OCD.", c(0,30), "cm")]) * 0.3 * 10
edgeroi.grids$OCS.30cm.f = ifelse(edgeroi.grids$OCS.30cm>76, 76, ifelse(edgeroi.grids$OCS.30cm<28, 28, edgeroi.grids$OCS.30cm))
## plot OCS 0-30 cm and the error map:
png(file = "Fig_RF_organic_carbon_stock_Edgeroi.png", width = 900, height = 1100, res=120)
par(mfrow=c(2,1), oma=c(0,0,0,1), mar=c(0,0,3.5,1.5))
plot(raster(edgeroi.grids["OCS.30cm.f"]), col=leg, main=paste0("Organic carbon stock 0", "\U2012", "30 cm (t/ha)"), axes=FALSE, box=FALSE, zlim=c(28,76))
points(edgeroi.sp, pch=21, bg="white", cex=.8)
plot(raster(edgeroi.grids["OCD.30cm_se"])*0.3*10, col=rev(bpy.colors()), main="Standard prediction error (t/ha)", axes=FALSE, box=FALSE)
points(edgeroi.sp, pch=21, bg="white", cex=.8)
dev.off()
## Total OCS per land use (http://data.environment.nsw.gov.au/dataset/nsw-landuseac11c):
edgeroi.grids$LandUse = readGDAL("edgeroi_LandUse.sdat")$band1
lu.leg = read.csv("LandUse.csv")
edgeroi.grids$LandUseClass = paste(join(data.frame(LandUse=edgeroi.grids$LandUse), lu.leg, match="first")$LU_NSWDeta)
OCS_agg.lu <- plyr::ddply(edgeroi.grids@data, .(LandUseClass), summarize, Total_OCS_kt=round(sum(OCS.30cm*250^2/1e4, na.rm=TRUE)/1e3), Area_km2=round(sum(!is.na(OCS.30cm))*250^2/1e6))
OCS_agg.lu$LandUseClass.f = strtrim(OCS_agg.lu$LandUseClass, 34)
OCS_agg.lu$OCH_t_ha_M = round(OCS_agg.lu$Total_OCS_kt*1000/(OCS_agg.lu$Area_km2*100))
OCS_agg.lu[OCS_agg.lu$Area_km2>5,c("LandUseClass.f","Total_OCS_kt","Area_km2","OCH_t_ha_M")]
sum(OCS_agg.lu$Total_OCS_kt, na.rm = TRUE)
setwd(md)

## 5. Spatiotemporal modeling of OCD (USA48) ----
## 250,428 observations of organic carbon density (kg/m3) together with some 80+ covariates
setwd("./USA48")
OCD_stN = readRDS("usa48.OCD_spacetime_matrix.rds")
pr.lst <- names(OCD_stN)[-which(names(OCD_stN) %in% c("SOURCEID", "DEPTH.f", "OCDENS",   "YEAR", "YEAR_c", "LONWGS84", "LATWGS84"))]
## 132 vars
fm0.st <- as.formula(paste('OCDENS ~ DEPTH.f + ', paste(pr.lst, collapse="+")))
sel0.m = complete.cases(OCD_stN[,all.vars(fm0.st)])
rf0.OCD_st <- ranger(fm0.st, data=OCD_stN[sel0.m,all.vars(fm0.st)], importance="impurity", write.forest=TRUE, num.trees=120)
# Growing trees.. Progress: 31%. Estimated remaining time: 1 minute, 11 seconds.
# Growing trees.. Progress: 71%. Estimated remaining time: 27 seconds.
rf0.OCD_st
## RMSE = 10 kg/m3 / 60%
## Most important variables:
xl <- as.list(ranger::importance(rf0.OCD_st))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))

## Map the decrease of OCD:
library(greenbrown)
library(raster)
## import and re-organize predictions:
tif.lst = list.files(pattern="_10km.tif")
g10km = as(readGDAL(tif.lst[1]), "SpatialPixelsDataFrame")
for(i in 2:length(tif.lst)){ g10km@data[,i] = readGDAL(tif.lst[i], silent=TRUE)$band1[g10km@grid.index] }
names(g10km) = basename(tif.lst)
g10km = as.data.frame(g10km)
gridded(g10km) = ~x+y
proj4string(g10km) = "+proj=longlat +datum=WGS84"
## Focus on Texas
library(maps)
library(maptools)
states <- map('state', plot=FALSE, fill=TRUE)
states = SpatialPolygonsDataFrame(map2SpatialPolygons(states, IDs=1:length(states$names)), data.frame(names=states$names))
proj4string(states) = "+proj=longlat +datum=WGS84"
ov.g10km = over(y=states, x=g10km)
txg10km = g10km[which(ov.g10km$names=="texas"),]
txg10km = as.data.frame(txg10km)
gridded(txg10km) = ~x+y
proj4string(txg10km) = "+proj=longlat +datum=WGS84"
spplot(log1p(stack(txg10km)), col.regions=SAGA_pal[[1]])
g10km.b = raster::brick(txg10km)
#names(g10km.b)
trendmap <- TrendRaster(g10km.b, start=c(1935, 1), freq=1, breaks=1) ## can be computationally intensive
#plot(trendmap, col=brgr.colors(20), legend.width=2)
plot(trendmap[["SlopeSEG1"]], col=rev(SAGA_pal[["SG_COLORS_GREEN_GREY_RED"]]), zlim=c(-1.5,1.5), main="Slope SEG1")
trendclassmap <- TrendClassification(trendmap, min.length=8, max.pval=0.05)
plot(trendclassmap)
setwd(md)
