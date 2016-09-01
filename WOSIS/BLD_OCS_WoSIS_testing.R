## Testing of PedoTransfer Functions for Bulk density based on WoSIS data
## By Tom.Hengl@isric.org and Amanda Rachmaran
## Pedo-Transfer Function for BLD (http://gsif.isric.org/doku.php/wiki:soil_data#using_machine_learning_to_build_pedo-transfer-functions)

library(R.utils)
library(raster)
library(rgdal)
library(snowfall)
library(utils)
library(tools)
library(plyr)
library(maps)
library(maptools)
library(ranger)
library(randomForestSRC)
library(parallel)
options(rf.cores=detectCores(), mc.cores=detectCores())

country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")

## Point data:
load("/data/models/SPROPS/ovA.rda")
## Add Bulk density from SoilGrids
gm <- readRDS(paste0("/data/models/SPROPS/mrf.BLD.rds"))
sel.na <- complete.cases(ovA[,gm$forest$independent.variable.names])
ovA$BLD.f = NA
ovA[sel.na,"BLD.f"] <- predict(gm, ovA[sel.na,])$predictions

## Only soil properties:
rfsrc_BD0 <- rfsrc(BLD ~ ORCDRC + PHIHOX + SNDPPT + CLYPPT + CRFVOL + DEPTH.f, data=ovA)
rfsrc_BD0
## Soil properties and SoilGrids BLD:
rfsrc_BD1 <- rfsrc(BLD ~ ORCDRC + PHIHOX + SNDPPT + CLYPPT + CRFVOL + DEPTH.f + BLD.f, data=ovA)
rfsrc_BD1

tvars = c("SOURCEID","BLD","ORCDRC", "PHIHOX", "SNDPPT", "CLYPPT", "CRFVOL", "DEPTH.f", "BLD.f","LONWGS84","LATWGS84")
sel.t = complete.cases(ovA[,tvars])
str(ovA[sel.t,tvars])
## 86,119 obs. of  11 variables
saveRDS(ovA[sel.t,tvars], "global_profiles_sept_2016.rds")

soilgrids_BD_xy = ovA[!duplicated(ovA$LOC_ID)&!is.na(ovA$BLD),c("LOC_ID","LONWGS84","LATWGS84")]
soilgrids_BD_xy = soilgrids_BD_xy[!is.na(soilgrids_BD_xy$LONWGS84),]
coordinates(soilgrids_BD_xy) = ~ LONWGS84+LATWGS84
proj4string(soilgrids_BD_xy) = CRS("+proj=longlat +datum=WGS84")
par(par(mar=c(.0,.0,.0,.0)))
plot(country, col="lightgrey", ylim=c(-61,82))
points(soilgrids_BD_xy, pch="+", cex=.5)

## Conclusion: global PTF using lab data is possible and will have an error of about +/- 160 kg / cubic-m
## If SoilGrids BLD is used as covariate, then the accuracy is even higher 

## Downalod WoSIS points:
download.file("http://www.isric.org/sites/default/files/datasets/WoSIS_2016_July.zip", "WoSIS_2016_July.zip")
system("7za e WoSIS_2016_July.zip")
lst = list.files(pattern=glob2rx("*_july.csv"))
wosis_07_2016 = lapply(list.files(pattern=glob2rx("*_july.csv")), read.csv, stringsAsFactors=FALSE)
names(wosis_07_2016) = file_path_sans_ext(lst)
#saveRDS(wosis_07_2016, "wosis_07_2016.rds")

wosis_BD = wosis_07_2016[["bulk_density_fine_earth_2016_july"]][,c("profile_layer_id","dataset_id","profile_code","date","top","bottom","method","latitude","longitude")]
wosis_BD$BLDFIE = wosis_07_2016[["bulk_density_fine_earth_2016_july"]]$value*1000
str(wosis_BD)
summary(duplicated(wosis_BD$profile_layer_id))
## many many duplicates
wosis_BD = wosis_BD[!duplicated(wosis_BD$profile_layer_id),]
wosis_BD$ORCDRC = plyr::join(wosis_BD, wosis_07_2016[["carbon_organic_2016_july"]], by="profile_layer_id", type="left", match = "first")$value
wosis_BD$PHIHOX = plyr::join(wosis_BD, wosis_07_2016[["ph_h2o_2016_july"]], by="profile_layer_id", type="left", match = "first")$value
wosis_BD$SNDPPT = plyr::join(wosis_BD, wosis_07_2016[["texture_sand_total_2016_july"]], by="profile_layer_id", type="left", match = "first")$value
wosis_BD$CLYPPT = plyr::join(wosis_BD, wosis_07_2016[["texture_clay_total_2016_july"]], by="profile_layer_id", type="left", match = "first")$value
wosis_BD$CRFVOL = plyr::join(wosis_BD, wosis_07_2016[["coarse_fragments_volumetric_total_2016_july"]], by="profile_layer_id", type="left", match = "first")$value
wosis_BD = wosis_BD[complete.cases(wosis_BD[,c("BLDFIE","ORCDRC","PHIHOX","SNDPPT","CLYPPT","CRFVOL")]),]
str(wosis_BD)
## 51,763 obs. of  14 variables
wosis_BD$DEPTH = wosis_BD$top + (wosis_BD$bottom - wosis_BD$top)/2
wosis_BD$LOC_ID = paste("ID", wosis_BD$longitude, wosis_BD$latitude, sep="_")
wosis_BD_xy = wosis_BD[!duplicated(wosis_BD$LOC_ID),c("LOC_ID","longitude","latitude")]
coordinates(wosis_BD_xy) = ~ longitude + latitude
proj4string(wosis_BD_xy) = CRS("+proj=longlat +datum=WGS84")
par(par(mar=c(.0,.0,.0,.0)))
plot(country, col="lightgrey", ylim=c(-61,82))
points(wosis_BD_xy, pch="+", cex=.5)
## 10,249 points only ??

## SoilGrids BLD:
sg.tifs = paste0("/data/GEOG/BLDFIE_M_sl", 1:7, "_250m_ll.tif")
sfInit(parallel=TRUE, cpus=length(sg.tifs))
sfExport("wosis_BD_xy", "sg.tifs")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(sg.tifs, function(i){try( raster::extract(raster(i), wosis_BD_xy) )}) 
snowfall::sfStop()
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = paste0("BLDFIE.sg_", 1:7)
## trapesoidal rule:

saveRDS(wosis_BD, "wosis_BD.rds")

## Only soil properties:
rfsrc_BD0 <- rfsrc(BLDFIE ~ ORCDRC + PHIHOX + SNDPPT + CLYPPT + CRFVOL + DEPTH, data=wosis_BD)
rfsrc_BD0
## Soil properties and SoilGrids BLD:
rfsrc_BD1 <- rfsrc(BLDFIE ~ ORCDRC + PHIHOX + SNDPPT + CLYPPT + CRFVOL + DEPTH + BLDFIE.sg, data=wosis_BD)
rfsrc_BD1
