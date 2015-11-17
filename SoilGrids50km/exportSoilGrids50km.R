# title         : exportSoilGrids50km.R
# purpose       : Aggregating SoilGrids1km to 0.5 grids for IPCC;
# reference     : SoilGrids1km are available for download from [http://www.soilgrids.org/]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : SoilGrids from [ftp://ftp.soilgrids.org/data/recent/]
# outputs       : Grid stack at 0.5 decimal degrees resolution (land mask = 67,484 pixels);
# remarks 1     : all predictions are at point support, hence the confidence limits represent mean upper and lower values for predictions at point support;
# remarks 2     : This script takes > 2 hours to run!

library(RSAGA)
library(R.utils)
library(snowfall)
library(raster)
library(plotKML)
library(GSIF)
indir <- "G:/SoilGrids1km/zipped"  ## all files on "ftp://ftp.soilgrids.org/data/recent/"
outdir <- "G:/SoilGrids1km/50km"
TAXGWRB.leg <- read.table("G:/SoilGrids1km/1km/TAXGWRB.txt", header=TRUE)
TAXOUSDA.leg <- read.table("G:/SoilGrids1km/1km/TAXOUSDA.txt", header=TRUE)
tifout.lst <- read.csv(file="tifout.lst.csv")
sel <- c(grep(tifout.lst$x, pattern="_L"), grep(tifout.lst$x, pattern="_U"))
tifout_s.lst <- as.vector(tifout.lst[-sel,"x"])
tifout.lst <- as.vector(tifout.lst[,"x"])
cellsize = 0.5
v = "02_apr_2014"
e = "50km"

## Resampling using SAGA GIS:
rsaga.env()
rsaga.get.usage("grid_tools", 0)

resample.grid <- function(x){
  inname <- paste0(x, "_", v, ".tif")
  fname <- paste(x, v, e, sep="_")
  if(!file.exists(paste0(fname, ".sgrd"))){
    gunzip(filename=paste0(indir, "/", inname, ".gz"), destname=paste0(outdir, "/", inname), remove=FALSE)
    unlink("tmp.sgrd"); unlink("tmp.sdat")
    tmp.file <- set.file.extension(tempfile(), ".sgrd")
    rsaga.geoprocessor(lib="io_gdal", module=0, param = list(GRIDS=tmp.file, FILES=inname, INTERPOL=0, TRANSFORM=FALSE))
    if(x %in% c("TAXGWRB", "TAXOUSDA")){
      rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT=tmp.file, TARGET=0, SCALE_UP_METHOD=9, USER_XMIN=-180+cellsize/2, USER_XMAX=180-cellsize/2, USER_YMIN=-90+cellsize/2, USER_YMAX=90-cellsize/2, USER_SIZE=cellsize, USER_GRID=paste0(fname, ".sgrd")),  show.output.on.console = FALSE)
    } else {
      rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT=tmp.file, TARGET=0, SCALE_UP_METHOD=5, USER_XMIN=-180+cellsize/2, USER_XMAX=180-cellsize/2, USER_YMIN=-90+cellsize/2, USER_YMAX=90-cellsize/2, USER_SIZE=cellsize, USER_GRID=paste0(fname, ".sgrd")),  show.output.on.console = FALSE)
    }
    ## recycle:
    unlink(tmp.file)
    unlink(inname)
  }
}

sfInit(parallel=TRUE, cpus=8)
sfLibrary(RSAGA)
sfLibrary(R.utils)
sfExport("indir", "outdir", "cellsize", "v", "e", "tifout.lst")
t <- sfLapply(tifout.lst, resample.grid)
sfStop()
## TAKES CA 1-2 hours!!

## Soil mask map:
download.file("http://worldgrids.org/lib/exe/fetch.php?media=smkisr3a.tif.gz", "smkisr3a.tif.gz")
gunzip("smkisr3a.tif.gz")
rsaga.geoprocessor(lib="io_gdal", module=0, param = list(GRIDS="tmp.sgrd", FILES="SMKISR3a.tif", INTERPOL=0, TRANSFORM=FALSE))
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="tmp.sgrd", TARGET=0, SCALE_UP_METHOD=9, USER_XMIN=-180+cellsize/2, USER_XMAX=180-cellsize/2, USER_YMIN=-90+cellsize/2, USER_YMAX=90-cellsize/2, USER_SIZE=cellsize, USER_GRID="SMKISR3a_50km.sgrd"))
## soil mask groups:
levs <- data.frame(ID=1:3, code=c("soils with vegetation cover", "urban areas", "bare soil areas"))
smask <- readGDAL("SMKISR3a_50km.sdat")
summary(as.factor(smask$band1))
smask$deserts <- ifelse(smask$band1==3, 1, NA)
smask$fillmask <- ifelse(!is.na(smask$band1), 1, NA)
writeGDAL(smask["deserts"], "deserts.sdat", "SAGA", mvFlag=-99999)
writeGDAL(smask["fillmask"], "fillmask.sdat", "SAGA", mvFlag=-99999)
unlink("SMKISR3a.tif")
smask.pix <- as(smask["band1"], "SpatialPixelsDataFrame")

## read all rasters to R and filter out all missing values etc:
sgrd.lst <- list.files(pattern=paste0("*", e, ".sdat$"))
s50km <- raster::stack(sgrd.lst) 
grid50km <- as(s50km, "SpatialGridDataFrame")
## TAKES >5 mins
#str(names(grid50km))

sel.f <- c("SMKISR3a_50km", "TAXGWRB_02_apr_2014_50km", "TAXOUSDA_02_apr_2014_50km")
## maps that do not need to be filtered:
tv <- names(grid50km)[-which(names(grid50km) %in% sel.f)]
for(j in 1:length(tv)){ ## TAKES > 5 mins
  ## filter memberships (fill in missing pixels):
  if(length(unlist(sapply(c("TAXGWRB","TAXOUSDA"), grep, x=tv[j])))>0){
    grid50km@data[,tv[j]] <- round(ifelse(is.na(grid50km@data[,tv[j]])&!is.na(smask$band1), 0, grid50km@data[,tv[j]]))
  }
  ## filter ORCDRC, SLTPPT, CLYPPT using empirical values
  if(length(unlist(sapply(c("ORCDRC","SLTPPT","CLYPPT","OCSTHA"), grep, x=tv[j])))>0){
    grid50km@data[,tv[j]] <- round(ifelse(smask$band1==3&!is.na(smask$band1), 0, grid50km@data[,tv[j]]))
  }
  if(length(unlist(sapply(c("SNDPPT"), grep, x=tv[j])))>0){
    grid50km@data[,tv[j]] <- round(ifelse(smask$band1==3&!is.na(smask$band1), 100, grid50km@data[,tv[j]]))
  }
  ## Filter all missing pixels using values from the edge of the deserts...
  if(length(unlist(sapply(c("BDRICM","PHIHOX","BDRLOG","CRFVOL","BLD","CEC"), grep, x=tv[j])))>0){
    unlink("flt.sgrd")
    rsaga.geoprocessor(lib="grid_tools", module=29, param=list(INPUT=paste0(tv[j],".sgrd"), MASK="fillmask.sgrd", RESULT="flt.sgrd", INTERPOLATION=4, START=0, START_SIZE=1), show.output.on.console = FALSE)
    img <- round(readGDAL("flt.sdat", silent=TRUE)$band1)
    grid50km@data[,tv[j]] <- img
    ## fix values outside the range:
    if(length(grep(pattern="_M", x=tv[j]))>0){
      grid50km@data[,tv[j]] <- ifelse(img<0, 0, img)
    }
    if(length(unlist(sapply(c("BDRLOG","CEC","CRFVOL"), grep, x=tv[j])))>0 & length(grep(pattern="_L", x=tv[j]))>0){
      grid50km@data[,tv[j]] <- ifelse(img<0, 0, img)
    }
    if(length(grep(pattern="PHIHOX", x=tv[j]))>0){
      grid50km@data[,tv[j]] <- ifelse(img<min(soil.legends$PHIHOX$MIN), min(soil.legends$PHIHOX$MIN), ifelse(img>max(soil.legends$PHIHOX$MAX), max(soil.legends$PHIHOX$MAX), img))
    }
    if(length(grep(pattern="BLD", x=tv[j]))>0){
      grid50km@data[,tv[j]] <- ifelse(img<min(soil.legends$BLD$MIN), min(soil.legends$BLD$MIN), ifelse(img>max(soil.legends$BLD$MAX), max(soil.legends$BLD$MAX), img))
    }
    if(length(unlist(sapply(c("BDRLOG","CRFVOL"), grep, x=tv[j])))>0 & length(grep(pattern="_U", x=tv[j]))>0){
      grid50km@data[,tv[j]] <- ifelse(img>100, 100, img)
    }
    if(length(grep(pattern="BDRICM", x=tv[j]))>0){
      grid50km@data[,tv[j]] <- ifelse(img>240, 240, ifelse(img<0, 0, img))
    }
  }
}

## Test plotting maps:
#plot(raster(grid50km["BDRLOG_L_02_apr_2014_50km"]), col=SAGA_pal[[1]])
#plot(raster(grid50km["PHIHOX_sd1_M_02_apr_2014_50km"]), col=SAGA_pal[[1]])
#plot(raster(grid50km["CRFVOL_sd1_M_02_apr_2014_50km"]), col=SAGA_pal[[1]])
#plot(raster(grid50km["TAXOUSDA_Aqualfs_02_apr_2014_50km"]), col=SAGA_pal[[10]])
#plot(raster(grid50km["TAXOUSDA_02_apr_2014_50km"]), col=SAGA_pal[[12]])

## write to a csv file:
grid50km.csv <- as.data.frame(grid50km[smask.pix@grid.index,])
grid50km.csv$SMKISR3a_50km <- as.factor(as.character(grid50km.csv$SMKISR3a_50km))
levels(grid50km.csv$SMKISR3a_50km) = levs$code
grid50km.csv$TAXGWRB_02_apr_2014_50km <- as.factor(as.character(grid50km.csv$TAXGWRB_02_apr_2014_50km))
levels(grid50km.csv$TAXGWRB_02_apr_2014_50km) <- TAXGWRB.leg$NAME[match(levels(grid50km.csv$TAXGWRB_02_apr_2014_50km), as.character(TAXGWRB.leg$MINIMUM))]
grid50km.csv$TAXOUSDA_02_apr_2014_50km <- as.factor(as.character(grid50km.csv$TAXOUSDA_02_apr_2014_50km))
levels(grid50km.csv$TAXOUSDA_02_apr_2014_50km) <- TAXOUSDA.leg$NAME[match(levels(grid50km.csv$TAXOUSDA_02_apr_2014_50km), as.character(TAXOUSDA.leg$MINIMUM))]

## fix columns:
names(grid50km.csv)[which(names(grid50km.csv) %in% c("s1","s2"))] = c("LONWGS84", "LATWGS84")
write.csv(grid50km.csv, file="SoilGrids50km.csv", na="")
gzip("SoilGrids50km.csv")

#require(ncdf)
#rnc <- writeRaster(s50km, filename='SoilGrids50km.nc', format="CDF", overwrite=TRUE, varname=, varunit=, longname=, xname="LONWGS84", yname="LATWGS84")


## end of script;