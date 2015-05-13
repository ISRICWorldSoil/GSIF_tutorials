## Model fitting:

library(aqp)
library(GSIF)
library(plotKML)
library(sp)
library(plyr)
library(randomForest)
library(gdalUtils)
library(raster)
# global
gdal.dir = shortPathName("C:\\ms4w")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)

dir <- "E:\\data\\soildata\\depth\\points\\models\\"
setwd(dir)
## point data:
load("../profs/NA.sp.rda")
## bounding box:
na.bbox <- matrix(c(-175,14,-58,70), nrow=2)

## covariates:
#cov.lst <- c("TWISRE3a.tif", "SLPSRT3a.tif", "DEMSRE3a.tif",
#            "L3POBI3b.tif", "GLTUHA3x.tif")
cov.lst <- c("na_TWISRE3a.tif", "na_SLPSRT3a.tif", "na_DEMSRE3a.tif",
            "na_L3POBI3b.tif", "na_GLTUHA3x.tif")
## extract values:
sel <- runif(length(NA.sp))<0.05
ov.NA <- GSIF::extract(NA.sp[sel,], cov.lst, 
        path="..\\covs")

rmat <- plyr::join(as.data.frame(NA.sp[sel,]), ov.NA, type="left", match="first")
rmat$na_L3POBI3b.tif <- as.factor(rmat$na_L3POBI3b.tif)
rmat$na_GLTUHA3x.tif <- as.factor(rmat$na_GLTUHA3x.tif)
fm <- as.formula(paste("BDRICM ~ ", paste(basename(cov.lst), collapse="+")))
rf.BDRICM <- randomForest(fm, data=rmat, na.action=na.omit, importance=TRUE)
rf.BDRICM
gc()
#19%, 20%, 21%, 22%, 25%, 27%,
varImpPlot(rf.BDRICM)

## prepare prediction locations:
for(j in 1:length(cov.lst)){
  if(!file.exists(paste0("../covs/", cov.lst[j]))){
    gdalwarp(paste0("G:/SoilGrids1km/1km/", cov.lst[j]), 
    paste0("G:/NA/covs/na_", cov.lst[j]), 
    te=as.vector(na.bbox), r="near")
  }
}
## Read to R:
na1km <- stack(paste0("../covs/", cov.lst))
## subset:
## extent(c(-85, -75, 40, 50); extent(c(-125, -115, 35, 45)
s.na1km <- crop(na1km, extent(c(-120, -115, 65, 70))) 
s.na1km <- as(s.na1km, "SpatialPixelsDataFrame")
sel.C <- complete.cases(s.na1km@data)
summary(sel.C)
s.na1km <- s.na1km[sel.C,]

names(s.na1km)[1:3] <- cov.lst[1:3]
## fix missing factors:
l.L3P <- data.frame(levs=levels(rmat$na_L3POBI3b.tif),
         na_L3POBI3b=as.numeric(levels(rmat$na_L3POBI3b.tif)))
l.GLT <- data.frame(levs=levels(rmat$na_GLTUHA3x.tif), 
        na_GLTUHA3x=as.numeric(levels(rmat$na_GLTUHA3x.tif)))
s.na1km$na_L3POBI3b.tif <- plyr::join(s.na1km@data["na_L3POBI3b"], 
        l.L3P, type="left", match = "first")$levs
s.na1km$na_GLTUHA3x.tif <- plyr::join(s.na1km@data["na_GLTUHA3x"], 
        l.GLT, type="left", match = "first")$levs
str(s.na1km@data)

xp <- predict(rf.BDRICM, s.na1km@data, na.action=na.pass)
s.na1km$BDRICM <- ifelse(xp<0, 0, xp)
writeGDAL(s.na1km["BDRICM"], "../out/na_BDRICM.tif", type="Int16")
BDRICM1km <- log1p(raster(s.na1km["BDRICM"]))
plotKML(BDRICM1km, colour_scale=SAGA_pal[[1]], 
    png.width = BDRICM1km@ncols*3, png.height = BDRICM1km@nrows*3)
