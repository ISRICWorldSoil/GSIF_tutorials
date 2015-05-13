## Deriving Available Water capacity using SoilGrids1km
## prepared by: T. Hengl (tom.hengl@wur.nl)

##-----------------------------------
## subset SoilGrids1km for Africa
##-----------------------------------

library(snowfall)
library(gdalUtils)
library(raster)
library(GSIF)
library(rgdal)
library(sp)

gdal.dir <- shortPathName("C:\\Program Files\\GDAL")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)

tifout_s.lst <- read.csv(file="tifout.lst.csv") ## list of files needed to run this analysis
indir <- "G:/SoilGrids1km/zipped"  ## download all files from "ftp://ftp.soilgrids.org/data/recent/"
outdir <- "G:/AFSIS/soilgrids"
cellsize = 1000
v = "02_apr_2014"
e = "1km"
## African projection system:
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## resample to 1 km Africa only:
resample.grid <- function(x){
  inname <- paste0(x, "_", v, ".tif")
  fname <- paste("af", x, v, e, sep="_")
  if(!file.exists(paste0(fname, ".tif"))){
    try( gunzip(filename=paste0(indir, "/", inname, ".gz"), destname=paste0(outdir, "/", inname), remove=FALSE) )
     if(x %in% c("TAXGWRB", "TAXOUSDA")){
      gdalwarp(inname, dstfile=paste0(fname, ".tif"), tr=c(cellsize,cellsize), te=c(-3977500,-4321500,3397500,3554500), r="near", t_srs=af.csy)
    } else {
      gdalwarp(inname, dstfile=paste0(fname, ".tif"), tr=c(cellsize,cellsize), te=c(-3977500,-4321500,3397500,3554500), r="bilinear", t_srs=af.csy)
    }
    unlink(paste0(outdir, "/", inname))
  }
}

sfInit(parallel=TRUE, cpus=8)
sfLibrary(gdalUtils)
sfLibrary(rgdal)
sfLibrary(R.utils)
sfExport("indir", "outdir", "cellsize", "v", "e", "tifout_s.lst", "af.csy")
t <- sfLapply(tifout_s.lst, resample.grid)
sfStop()

##-----------------------------------
## available Water capacity for Africa
##-----------------------------------

## derive total AWCh1, h2, h3 and WWP:
tv.lst <- c("SNDPPT", "SLTPPT", "CLYPPT", "ORCDRC", "BLD", "CEC", "PHIHOX")
tif.lst <- lapply(1:6, function(x){list.files(pattern=paste0("sd", x, "_M_02_apr"))})

## soil mask:
smask <- readGDAL("af_SMKISR3a.sdat")
smask <- as(smask, "SpatialPixelsDataFrame")
## TAKES >30 mins!
str(smask)
save(smask, file="smask.rda")
## 33 million pixels!

out <- smask
for(j in 1:length(tif.lst)){
  if(any(!file.exists(paste0("af_", c("AWCh1", "AWCh2", "AWCh3", "WWP"), "_sd", j, ".tif")))){
    r <- stack(tif.lst[[j]])
    r <- as(r, "SpatialGridDataFrame")
    names(r) <- c("BLD", "CEC", "CLYPPT", "ORCDRC", "PHIHOX", "SLTPPT", "SNDPPT")
    xd <- AWCPTF(r$SNDPPT[smask@grid.index], r$SLTPPT[smask@grid.index], r$CLYPPT[smask@grid.index], r$ORCDRC[smask@grid.index], r$BLD[smask@grid.index], r$CEC[smask@grid.index], r$PHIHOX[smask@grid.index]/10) ## TAKES >5 mins!
    out@data <- xd
    for(i in 1:ncol(xd)){
      out@data[,i] <- round(100*xd[,i])
      writeGDAL(out[i], paste0("af_", names(out)[i], "_sd", j, ".tif"), type="Byte", driver="GTiff", mvFlag=255)
    }
    rm(xd); rm(r)
  }
}

## Derive total AWC in mm for the whole profile:
awc.lst <- paste0("af_", c("AWCh1", "AWCh2", "AWCh3", "WWP")) 
sdepth <- c(0, 100*(-get("stdepths", envir=GSIF.opts)+get("stsize", envir=GSIF.opts)/2))

## read coarse fragments and depth to bedrock:
smask$af_BDRICM <- readGDAL(paste0("af_BDRICM_", v, "_1km.tif"))$band1[smask@grid.index]
for(j in 1:6){ smask@data[,paste0("af_CRFVOL_sd", j)] <- readGDAL(paste0("af_CRFVOL_sd", j, "_M_", v, "_1km.tif"))$band1[smask@grid.index] }
gc(); gc()

for(k in 1:length(awc.lst)){
  if(!file.exists(paste0(awc.lst[k], "_Total_mm_1km.tif"))){
    ## for each AWC type read all layers:
    for(j in 1:6){ smask@data[,paste0("WCT_sd", j)] <- readGDAL(paste0(awc.lst[k], "_sd", j, ".tif"))$band1[smask@grid.index] } 
    ## adjust available Water content for (1-CRF/100) and BDR:
    for(j in 1:6){
      svol <- (sdepth[j+1]-sdepth[j])*10
      ## ca. 30-40% of pixels have depth to bedrock within 0-100 cm!
      smask@data[,paste0("WCT_sdf",j)] <- ifelse(smask$af_BDRICM > sdepth[j+1], smask@data[,paste0("WCT_sd", j)]/100*svol*(1-smask@data[,paste0("af_CRFVOL_sd", j)]/100), ifelse(smask$af_BDRICM > sdepth[j], smask@data[,paste0("WCT_sd", j)]/100*svol*(smask$af_BDRICM-sdepth[j])/sdepth[j+1]*(1-smask@data[,paste0("af_CRFVOL_sd", j)]/100), 0))
    }
    ## sum up all layers
    smask@data[,"WCT_1km"] <- rowSums(smask@data[,paste0("WCT_sdf",1:6)], na.rm = TRUE) 
    if(is.numeric(smask@data[,"WCT_1km"])){
      writeGDAL(smask["WCT_1km"], paste0(awc.lst[k], "_Total_mm_1km.tif"), type="Int16", mvFlag=-9999)
    }
  }
}

## end of script;



