# title         : example_largeRaster.R
# purpose       : Working with large Rasters;
# reference     : Lecture at the GEOSTAT Quebec City 2013 [http://GEOSTAT-course.org]
# producer      : Prepared by T. Hengl (tom.hengl@wur.nl) and M. Kilibarda (milan.kili11@gmail.com)
# address       : In Quebec City, CA.
# inputs        : WorldGrid map at 1 km
# outputs       : SAGA GIS grids;
# remarks 1     : covariate data for spatial interpolation;

library(rgdal)
library(plotKML)
library(RCurl)
library(raster)
data(SAGA_pal)
library(splines)
library(nlme)

## Locate FWTools on the system:
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdalwarp
## Dowload 7z:
if(is.na(file.info("7za.exe")$size)){
  download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
  unzip("7za920.zip")
  unlink("7za920.zip")
}
si <- Sys.info()
aea.csy = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

## Download WorldGrids:
w.l <- c("tdmmod3a.tif", "tdlmod3a.tif")
for(j in 1:length(w.l)){
  is.na(file.info(w.l[j])$size){
  download.file(paste("http://worldgrids.org/lib/exe/fetch.php?media=", w.l[j], ".gz", sep=""), paste(w.l[j], ".gz", sep=""))
}
system("7za e td*mod3a.tif.gz")
GDALinfo("TDMMOD3a.tif")
## Huge raster -> rows 21600, columns 43200
## Reading this raster using 'readGDAL' is not a good idea!

## Area of interest (Spatial Prediction Competition GEOSTAT 2013):
## verify the certificate:
curl <- getCurlHandle()
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)

## download the table data:
cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdGd0RTJ4SGJVRU1IRWViYXN3WkhiaUE&single=true&gid=2&output=csv"), file="pnts.csv")
pnts <- read.csv("pnts.csv")
plot(pnts, asp=1)
names(pnts) <- c("LON", "LAT")
coordinates(pnts) <- ~LON+LAT
proj4string(pnts) = CRS(aea.csy)
pnts@bbox
pnts.ll <- spTransform(pnts, CRS("+proj=longlat +datum=WGS84"))

## use the raster package to 'connect' to the raster of interest:
TDMMOD3a <- raster("TDMMOD3a.tif")
TDLMOD3a <- raster("TDLMOD3a.tif")
TDMMOD3a
ext <- extent(-72, -70, 46.5, 48.5)
x <- crop(TDMMOD3a, ext)
x
plot(x, col=SAGA_pal[[1]])
## to prepare the tiles, consider using GSIF package:
tl <- getSpatialTiles(as(x, "SpatialPixelsDataFrame"), block.x=1)
lines(as(tl, "SpatialLines"))

## Overlay using raster package:
ov <- extract(TDMMOD3a, pnts.ll)
ov2 <- extract(TDLMOD3a, pnts.ll)
## 'extract' function is very efficient!!
pnts.xy <- SpatialPointsDataFrame(pnts, data.frame(TDMMOD3a = ov, TDLMOD3a = ov2, LON = pnts.ll@coords[,1], LAT = pnts.ll@coords[,2]))

## Build a regression model:
cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdGd0RTJ4SGJVRU1IRWViYXN3WkhiaUE&single=true&gid=0&output=csv"), file="cal.csv")
cal <- read.csv("cal.csv")
coordinates(cal) <- ~LON+LAT
proj4string(cal) = CRS("+proj=longlat +datum=WGS84")
cal.xy <- spTransform(cal, CRS(aea.csy))
cal.xy$LON <- cal@coords[,1]
cal.xy$LAT <- cal@coords[,2]
sel = cal.xy@coords[,1] > pnts@bbox[1,1] & cal.xy@coords[,1] < pnts@bbox[1,2] & cal.xy@coords[,2] > pnts@bbox[2,1] & cal.xy@coords[,2] < pnts@bbox[2,2]
cal.xy <- cal.xy[sel,]
## Overlay again:
ov <- extract(TDMMOD3a, cal[sel,])
ov2 <- extract(TDLMOD3a, cal[sel,])
## subset to the study area of interest:
cal.xy$TDMMOD3a <- ov
cal.xy$TDLMOD3a <- ov2
cal.xy <- remove.duplicates(cal.xy)
## regression diagnostic:
pairs(TEMPC~TDMMOD3a+TDLMOD3a+LAT+LON, cal.xy)
m = TEMPC~TDMMOD3a*TDLMOD3a*LAT+ns(LON, df=5)
TEMPC.m <- lm(m, cal.xy)
summary(TEMPC.m)
## the model explains 72% of variability
## optional: fit a GLS model using the 'nlme' package ? Takes > 20 minutes!
#TEMPC.m0 <- gls(TEMPC~TDMMOD3a+LAT, cal.xy[-TEMPC.m$na.action,], corExp(form = ~LON+LAT, nugget=TRUE))

## fit Vgm for residuals:
svar <- variogram(m, cal.xy[-TEMPC.m$na.action,])
plot(svar)
TEMPC.vgm <- fit.variogram(svar, vgm(nugget=5, model="Lin"))
plot(svar, TEMPC.vgm)

## check predictions:
TEMPC.uk.cv <- krige.cv(m, locations=cal.xy[-TEMPC.m$na.action,], TEMPC.vgm, nfold=5)
spplot(TEMPC.uk.cv["zscore"], col.regions=SAGA_pal[[1]]) ## some points are very difficult to predict!
sqrt(mean((TEMPC.uk.cv$var1.pred - TEMPC.uk.cv$observed)^2, na.rm = TRUE))
## RMSE = +/- 3.0 C

## Run predictions:
TEMPC.uk <- krige(m, locations=cal.xy[-TEMPC.m$na.action,], newdata=pnts.xy, TEMPC.vgm)
#spplot(TEMPC.uk[1], col.regions=SAGA_pal[[1]], sp.layout=list("sp.points", cal.xy, pch="+", col="black"))
## compare with actual measured values:
cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdGd0RTJ4SGJVRU1IRWViYXN3WkhiaUE&single=true&gid=3&output=csv"), file="preds.csv")
preds <- read.csv("preds.csv")
preds$TomHengl <- TEMPC.uk$var1.pred
lapply(preds[,5:14], function(x){sqrt(mean((preds$Measured - x)^2, na.rm = TRUE))}) 

## subset using FWTools:
unlink("na_TDMMOD3a.tif")
system(paste(gdalwarp, 'TDMMOD3a.tif', 'na_TDMMOD3a.tif', '-te', pnts.ll@bbox[1,1], pnts.ll@bbox[2,1], pnts.ll@bbox[1,2], pnts.ll@bbox[2,2]))
## Still a large file!
## resample to 5 times coarser resolution:
unlink("na_TDMMOD3a_20km.sdat")
system(paste(gdalwarp, 'na_TDMMOD3a.tif', 'na_TDMMOD3a_20km.sdat', '-r near', '-srcnodata -32768', '-dstnodata -99999', '-of SAGA', '-tr', res(TDMMOD3a)[1]*20, res(TDMMOD3a)[2]*20))
na_TDMMOD3a_20km <- readGDAL("na_TDMMOD3a_20km.sdat")
spplot(na_TDMMOD3a_20km, col.regions=SAGA_pal[[1]])

## clean-up:
unlink("TDMMOD3a.tif")
unlink("TDLMOD3a.tif")

## end of script;