# title         : tutorial_SoilGrids.R
# purpose       : Accessing and using SoilGrids;
# reference     : [http://gsif.isric.org/doku.php?id=wiki:tutorial_soilgrids]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : SoilGrids GeoTiffs [http://www.isric.org/content/soilgrids]
# outputs       : visualization of original and derived soil properties;

##-----------------------------------
## Accessing data
##-----------------------------------

## (a) FTP download:
library(RCurl)
## location of soilgrids:
sg.ftp <- "ftp://soilgrids:soilgrids@ftp.soilgrids.org/data/recent/"
filenames = getURL(sg.ftp, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames = strsplit(filenames, "\r*\n")[[1]]
filenames[1:5]

## download to a local directory:
ORC.name <- filenames[grep(filenames, pattern="ORCDRC_sd1_M")]
ORC.name
try(download.file(paste(sg.ftp, ORC.name[1], sep=""), ORC.name))
library(R.utils)
gunzip(ORC.name[1])

## check that everything is OK:
ORC.tif <- strsplit(ORC.name[1], ".gz")[[1]][1]
library(rgdal)
GDALinfo(ORC.tif)

## We focus on Ghana
wg.url <- url("http://gsif.isric.org/lib/exe/fetch.php?media=admin.af.rda")
load(wg.url)
proj4string(admin.af) <- "+proj=longlat +datum=WGS84"
country.af <- as(admin.af, "SpatialLines")
## Ghana bounding box:
ghana <- admin.af[admin.af$FORMAL_EN=="Republic of Ghana",]
ghana@bbox

## load soil Africa Soil Profile DB:
library(GSIF)
data(afsp)
sites <- afsp$sites
coordinates(sites) <- ~ LONWGS84 + LATWGS84
proj4string(sites) <- "+proj=longlat +datum=WGS84"
#af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
#sites.af <- spTransform(sites, CRS(af.csy))

## plot country and profiles
plot(ghana, col="red", lwd=2, asp=1)
lines(country.af)
points(sites, pch="+")
## in local projection system:
#ghana.xy <- spTransform(ghana, CRS(af.csy))
#ghana.xy@bbox

## get only Ghana:
library(raster)
library(gdalUtils)
getOption("gdalUtils_gdalPath")
## path is unknown, we have to help gdalUtils locate the software:
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)
getOption("gdalUtils_gdalPath")
te = as.vector(ghana@bbox)
unlink("ORC_sd1_Ghana.tif")
gdalwarp(ORC.tif, dstfile="ORC_sd1_Ghana.tif", te=te)
ORCDRC_sd1_ghana <- readGDAL("ORC_sd1_Ghana.tif")
plot(raster(ORCDRC_sd1_ghana))

## (b) Web Coverage Service
library(XML)
## location of service:
wcs = "http://wms3.isric.org/geoserver/soilgrids1km/wcs?"
## create an XML file:
l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", "ORCDRC_sd1_M", parent=l1)
l1
xml.out = "ORCDRC_sd1_M.xml"
saveXML(l1, file=xml.out)

## check if the layer exists:
gdalinfo(xml.out)

## Alternative: calculate offset and region dims:
bb <- matrix(nrow=2, c(-180,-90,180,90))
o.x = 43200 + round(43200*(te[1]-bb[1,2])/(bb[1,2]-bb[1,1]))
o.y = round(21600*(bb[2,2]-te[4])/(bb[2,2]-bb[2,1]))
d.y = round(21600*(te[4]-te[2])/(bb[2,2]-bb[2,1]))
d.x = round(43200*(te[3]-te[1])/(bb[1,2]-bb[1,1]))
o.x; o.y; d.x; d.y
## read only the piece of TIF file to R:
#ORCDRC_sd1_ghana <- readGDAL(ORC.tif, offset=c(o.y, o.x), region.dim=c(d.y, d.x), silent=TRUE)
## getCoverage using GDAL translate:
gdal_translate(src_dataset=xml.out, dst_dataset="ORC_sd1_Ghana.tif", of="GTiff", srcwin=paste(c(o.x, o.y, d.x, d.y), collapse=" "))
GDALinfo("ORC_sd1_Ghana.tif")
ORCDRC_sd1_ghana <- readGDAL("ORC_sd1_Ghana.tif")

## NOT Run:
#ORCDRC_sd1_ghana <- readGDAL(xml.out, offset=c(o.y, o.x), region.dim=c(d.y, d.x), silent=TRUE)

## plot the map:
data(soil.legends)
class.labels = paste(soil.legends[["ORCDRC"]]$MIN, "\226", soil.legends[["ORCDRC"]]$MAX)
ORCDRC_sd1_ghana$val <- cut(ORCDRC_sd1_ghana$band1, breaks=c(soil.legends[["ORCDRC"]]$MIN[1], soil.legends[["ORCDRC"]]$MAX), labels = class.labels)
bnd <- list(list("sp.points", sites, pch="+", col="black"), list("sp.lines", country.af))
spplot(ORCDRC_sd1_ghana["val"], col.regions=soil.legends[["ORCDRC"]]$COLOR, sp.layout=bnd, scales=list(draw=TRUE), main="Organic carbon in permilles (0\2265 cm)")
## Visualization of maps in Google Earth:
library(plotKML)
plotKML(ORCDRC_sd1_ghana["val"], colour_scale=soil.legends[["ORCDRC"]]$COLOR)

## plot soil types:
t1 <- newXMLNode("WCS_GDAL")
t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
t1.l <- newXMLNode("CoverageName", "TAXOUSDA", parent=t1)
saveXML(t1, file="TAXOUSDA.xml")
gdal_translate(src_dataset="TAXOUSDA.xml", dst_dataset="TAXOUSDA_Ghana.tif", of="GTiff", srcwin=paste(c(o.x, o.y, d.x, d.y), collapse=" "))
TAXUSDA_ghana <- readGDAL("TAXOUSDA_Ghana.tif")
TAXUSDA_ghana$band1 <- ifelse(TAXUSDA_ghana$band1==255, NA, TAXUSDA_ghana$band1)
summary(as.factor(TAXUSDA_ghana$band1))
TAXUSDA_ghana$val <- cut(TAXUSDA_ghana$band1, breaks=c(-1, soil.legends[["TAXOUSDA"]]$Number), labels=soil.legends[["TAXOUSDA"]]$Group)
xs <- summary(TAXUSDA_ghana$val)
sort(xs[-length(xs)], decreasing=TRUE)[1:5]
spplot(TAXUSDA_ghana["val"], col.regions=soil.legends[["TAXOUSDA"]]$COLOR, sp.layout=bnd, scales=list(draw=TRUE), main="Soil suborders based on USDA Taxonomy")

##-----------------------------------
## Deriving Soil Organic Carbon Stock
##-----------------------------------

## layers of interest:
var.name <- c(paste("ORCDRC_sd", c(1,2,3), "_M", sep=""), paste("BLD_sd", c(1,2,3), "_M", sep=""), paste("CRFVOL_sd", c(1,2,3), "_M", sep=""))

## now run everything in a loop...
for(i in 1:length(var.name)){
 t1 <- newXMLNode("WCS_GDAL")
 t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
 t1.l <- newXMLNode("CoverageName", var.name[i], parent=t1)
 xml.out <- paste(var.name[i], ".xml", sep="")
 saveXML(t1, file=xml.out)
 f.name <- paste(var.name[i], "_Ghana.tif", sep="")
 if(!file.exists(f.name)){
   gdal_translate(src_dataset=xml.out, dst_dataset=f.name, of="GTiff", srcwin=paste(c(o.x, o.y, d.x, d.y), collapse=" "))
 }
}

## read all layers:
library(raster)
ghana1km <- stack(paste(var.name, "_Ghana.tif", sep=""))
ghana1km <- as(ghana1km, "SpatialGridDataFrame")
str(ghana1km@data)
## Formula from: http://www.eoearth.org/view/article/156087/
soc <- function(ORCDRC, BLD, CRFVOL, HSIZE){ ORCDRC/1000 * HSIZE * BLD * (100-CRFVOL)/100 }
## test it:
soc(50, 1200, 12, .3)

sds <- get("stsize", envir = GSIF.opts)
sds
for(j in c(1,2,3)){
 ghana1km@data[,paste("SOC_sd", j, "_M", sep="")] <- soc(ghana1km@data[,paste("ORCDRC_sd", j, "_M_Ghana",sep="")], ghana1km@data[,paste("BLD_sd", j, "_M_Ghana", sep="")], ghana1km@data[,paste("CRFVOL_sd", j, "_M_Ghana", sep="")], HSIZE=sds[j])
}
## aggregate values:
ghana1km$SOC <- rowSums(ghana1km@data[,paste("SOC_sd", c(1,2,3), "_M", sep="")], na.rm=TRUE)
## mask out missing pixels:
ghana1km$SOC <- ifelse(ghana1km$ORCDRC_sd1_M_Ghana==-9999, NA, ghana1km$SOC)
rg <- quantile(ghana1km$SOC, c(0.01, 0.99), na.rm=TRUE)
at <- expm1(seq(log1p(rg[1]), log1p(rg[2]), length.out=20))
library(plotKML)
data(R_pal)
ghana1km$SOCf <- ifelse(ghana1km$SOC<rg[1], rg[1], ifelse(ghana1km$SOC>rg[2], rg[2], ghana1km$SOC))
spplot(ghana1km["SOCf"], at=at, col.regions=R_pal[["soc_pal"]], sp.layout=bnd, scales=list(draw=TRUE), main="Total soil carbon stock (0\22630 cm) in kg per square meter")

##-----------------------------------
## REST service
##-----------------------------------

## Not run: 
library(rjson)
library(sp)
## 2 points:
pnts <- data.frame(lon=c(10.65,5.36), lat=c(51.81,52.48), id=c("p1","p2"))
coordinates(pnts) <- ~lon+lat
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
pnts

sel <- sites@coords[,1]<ghana@bbox[1,2] & sites@coords[,1]>ghana@bbox[1,1] & sites@coords[,2]<ghana@bbox[2,2] & sites@coords[,2]>ghana@bbox[1,2]
pnts <- sites[sel,][400:500,]
## REST example:
soilgrids.r <- REST.SoilGrids(c("ORCDRC","PHIHOX"))
ov <- over(soilgrids.r, pnts)
str(ov)

kml(pnts, colour=SOURCEID, file="PHIHOX_depth.kml", 
    shape=paste("PHIHOX_depth_", 1:nrow(ov), ".png", sep=""), 
    size=6, points_names=pnts$SOURCEID, 
    colour_scale=rep("#FFFFFF", 2))

## plot soil depth curve:
ORCDRC.pnt1 <- data.frame(
  top=unlist(ov[1,grep("depthCodesMeters", names(ov))])*-100, 
  M=unlist(ov[1,grep("ORCDRC.M", names(ov))]), 
  L=unlist(ov[1,grep("ORCDRC.L", names(ov))]), 
  U=unlist(ov[1,grep("ORCDRC.U", names(ov))]))
ORCDRC.pnt1$variable <- "ORCDRC"
## plot the result:
library(lattice)
library(aqp)
data(soil.legends)
## Soil organic carbon:
ORCDRC.range = range(soil.legends[["ORCDRC"]]$MIN, soil.legends[["ORCDRC"]]$MAX)
dev.new(width=5, height=6)
xyplot(top ~ M | variable, data=ORCDRC.pnt1, ylab='Depth in cm',
  xlab='5th and 95th percentiles', xlim=ORCDRC.range,
  lower=ORCDRC.pnt1$L, upper=ORCDRC.pnt1$U, ylim=c(150,0),
  panel=panel.depth_function,
  alpha=0.25, sync.colors=TRUE,
  par.settings=list(superpose.line=list(col='RoyalBlue', lwd=3)),
  strip=strip.custom(bg=grey(0.8))
)

##-----------------------------------
## plot texture classes 15-30 cm
##----------------------------------- 

library(GSIF)
library(XML)
library(gdalUtils)
library(raster)
te.as <- c(17, -10, 19, -9)
o.x2 = 43200 + round(43200*(te.as[1]-bb[1,2])/(bb[1,2]-bb[1,1]))
o.y2 = round(21600*(bb[2,2]-te.as[4])/(bb[2,2]-bb[2,1]))
d.y2 = round(21600*(te.as[4]-te.as[2])/(bb[2,2]-bb[2,1]))
d.x2 = round(43200*(te.as[3]-te.as[1])/(bb[1,2]-bb[1,1]))
tex.name <- c("CLYPPT_sd3_M", "SNDPPT_sd3_M", "SLTPPT_sd3_M")

## run everything in a loop...
for(i in 1:length(tex.name)){
t1 <- newXMLNode("WCS_GDAL")
t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
t1.l <- newXMLNode("CoverageName", tex.name[i], parent=t1)
xml.out <- paste(tex.name[i], ".xml", sep="")
saveXML(t1, file=xml.out)
f.name <- paste(tex.name[i], "_angola.tif", sep="")
if(!file.exists(f.name)){
 gdal_translate(src_dataset=xml.out, dst_dataset=f.name, of="GTiff", srcwin=paste(c(o.x2, o.y2, d.x2, d.y2), collapse=" "))
}
}
angola1km <- stack(paste(tex.name, "_angola.tif", sep=""))
angola1km <- as(angola1km, "SpatialGridDataFrame")
library(soiltexture)
## example:
TT.points.in.classes(tri.data=data.frame(CLAY=20, SILT=40, SAND=40), class.sys="USDA.TT", PiC.type = "t")

## now run this per row:
names(angola1km@data) <- c("CLAY", "SAND", "SILT")
TT <- TT.points.in.classes(tri.data=angola1km@data, class.sys="USDA.TT", PiC.type="t", tri.sum.tst=FALSE)
## many transitional classes! Needs to be filtered...
no.TT <- sapply(TT, function(x){length(strsplit(x, ",")[[1]])})
sel.TT <- which(no.TT > 1)
trim <- function (x){ gsub("^\\s+|\\s+$", "", x) }
angola1km$TEX_sd3_M <- TT
## replace classes (transitional) with a random class:
for(i in sel.TT){
  angola1km$TEX_sd3_M[[i]] <- trim(strsplit(TT[i], ",")[[1]][ceiling(runif(1)*no.TT[i])])
}
angola1km$TEX_sd3_M <- as.factor(angola1km$TEX_sd3_M)
summary(angola1km$TEX_sd3_M)
#plot(raster(angola1km["TEX_sd3_M"]))
levels(angola1km$TEX_sd3_M)
## colors from http://gsif.isric.org/KML/soil_triangle.png
leg.tex <- c(rgb(240/255,190/255,153/255), rgb(195/255,161/255,136/255), rgb(202/255,195/255,169/255), rgb(195/255,188/255,160/255), rgb(228/255,226/255,203/255), rgb(212/255,190/255,161/255), rgb(209/255,203/255,162/255))
library(plotKML)
kml_open("angola1km_tex.kml")
kml_layer(angola1km, colour=TEX_sd3_M, colour_scale=leg.tex)
kml_screen("http://gsif.isric.org/KML/soil_triangle.png", position = "LR", sname = "triangle")
kml_close("angola1km_tex.kml")

## end of script;