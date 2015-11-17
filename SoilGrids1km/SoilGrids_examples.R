## Example of how to aggregate values from SoilGrids to 0--20 cm depths
## prepared by: T. Hengl (tom.hengl@wur.nl)
## requested by: Badza, Taruvinga <taruvinga.badza@wur.nl> and Schut, Tom <tom.schut@wur.nl> 

library(rgdal)
library(XML)
library(sp)

## We focus on Rwanda
wg.url <- url("http://gsif.isric.org/lib/exe/fetch.php?media=admin.af.rda")
load(wg.url)
proj4string(admin.af) <- "+proj=longlat +datum=WGS84"
country.af <- as(admin.af, "SpatialLines")
## rwanda bounding box:
rwanda <- admin.af[admin.af$FORMAL_EN=="Republic of Rwanda",]
rwanda@bbox

## load soil Africa Soil Profile DB:
sites <- read.csv("truvinga_badza_XY_Rwanda.csv")
sites$ID <- paste0("ID", 1:nrow(sites)) 
coordinates(sites) <- ~ POINT_X + POINT_Y
proj4string(sites) <- "+proj=longlat +datum=WGS84"
names(sites)
## plot country and profiles
plot(rwanda, col="red", lwd=2, asp=1)
lines(country.af)
points(sites, pch="+")

## get only rwanda:
library(gdalUtils)
fw <- utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdal_setInstallation(search_path=paste(fw, "bin", sep="\\"), rescan=TRUE)
te = as.vector(sites@bbox)

## layers of interest:
var.name <- c(paste("ORCDRC_sd", c(1,2,3), "_M", sep=""), paste("CLYPPT_sd", c(1,2,3), "_M", sep=""))
## Bounding box of interest (row column coordinates)
bb <- matrix(nrow=2, c(-180,-90,180,90))
o.x = 43200 + round(43200*(te[1]-bb[1,2])/(bb[1,2]-bb[1,1]))
o.y = round(21600*(bb[2,2]-te[4])/(bb[2,2]-bb[2,1]))
d.y = round(21600*(te[4]-te[2])/(bb[2,2]-bb[2,1]))
d.x = round(43200*(te[3]-te[1])/(bb[1,2]-bb[1,1]))
o.x; o.y; d.x; d.y

## now run everything in a loop...
wcs = "http://wms3.isric.org/geoserver/soilgrids1km/wcs?"
for(i in 1:length(var.name)){
 t1 <- newXMLNode("WCS_GDAL")
 t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
 t1.l <- newXMLNode("CoverageName", var.name[i], parent=t1)
 xml.out <- paste(var.name[i], ".xml", sep="")
 saveXML(t1, file=xml.out)
 f.name <- paste(var.name[i], "_rwanda.tif", sep="")
 if(!file.exists(f.name)){
   gdal_translate(src_dataset=xml.out, dst_dataset=f.name, of="GTiff", srcwin=paste(c(o.x, o.y, d.x, d.y), collapse=" "))
 }
}

## read all layers:
library(raster)
rwanda1km <- stack(paste(var.name, "_rwanda.tif", sep=""))
rwanda1km <- as(rwanda1km, "SpatialGridDataFrame")
str(rwanda1km@data)

## overlay points:
proj4string(rwanda1km) <- "+proj=longlat +datum=WGS84"
ov <- over(sites, rwanda1km)
str(ov)
ov.xy <- cbind(data.frame(sites), ov)
## create SoilProfiles:
upper=c(rep(0,nrow(sites)),rep(5,nrow(sites)),rep(15,nrow(sites)))
lower=c(rep(5,nrow(sites)),rep(15,nrow(sites)),rep(30,nrow(sites)))
ORCDRC=c(ov$ORCDRC_sd1_M_rwanda, ov$ORCDRC_sd2_M_rwanda, ov$ORCDRC_sd3_M_rwanda)
CLYPPT=c(ov$CLYPPT_sd1_M_rwanda, ov$CLYPPT_sd2_M_rwanda, ov$CLYPPT_sd3_M_rwanda)

library(aqp)
library(plyr)
horizons <- data.frame(ID=paste0("ID", rep(1:nrow(sites), 3)), upper, lower, ORCDRC, CLYPPT)
spc <- merge(data.frame(sites), horizons, type='inner')
depths(spc) <- ID ~ upper + lower
site(spc) <- ~ POINT_X + POINT_Y
## fit splines:
library(GSIF)
ORC.s <- mpspline(spc, var.name="ORCDRC", d=t(c(0,5,15,30)), mxd=30)
CLY.s <- mpspline(spc, var.name="CLYPPT", d=t(c(0,5,15,30)), mxd=30)
sortl <- match(sites$ID, spc@site[,"ID"])
## aggregate:
ov.xy$ORC_0_20cm <- signif(colMeans(ORC.s$var.1cm[1:20,sortl], na.rm = TRUE), 3)
ov.xy$CLY_0_20cm <- signif(colMeans(CLY.s$var.1cm[1:20,sortl], na.rm = TRUE), 3)
write.csv(ov.xy, file="truvinga_badza_XY_Rwanda_filled.csv")

## end of script;