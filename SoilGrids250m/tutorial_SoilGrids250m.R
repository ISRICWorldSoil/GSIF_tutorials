# SoilGrids tutorial
# Reference: [http://gsif.isric.org/doku.php?id=wiki:tutorial_soilgrids]
# Tom.Hengl@isric.org

library(RCurl)
library(rgdal)
library(GSIF)
library(raster)
library(plotKML)
library(XML)
library(lattice)
library(aqp)
library(soiltexture)

## GDAL paths:
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
  gdalinfo <- paste0(gdal.dir, "/gdalinfo.exe")
} else {
  gdal_translate = "gdal_translate"
  gdalwarp = "gdalwarp"
  gdalinfo = "gdalinfo"
}

##-----------------------------------
## Accessing data
##-----------------------------------

## (a) FTP download:
## location of soilgrids:
sg.ftp <- "ftp://ftp.soilgrids.org/data/recent/"
filenames = getURL(sg.ftp, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames = strsplit(filenames, "\r*\n")[[1]]
filenames[1:5]

## download to a local directory:
ORC.name <- filenames[grep(filenames, pattern=glob2rx("ORCDRC_M_sl1_250m_ll.tif$"))]
ORC.name
try(download.file(paste(sg.ftp, ORC.name, sep=""), ORC.name))
## 2.8GB Geotiff!!

## check that everything is OK:
GDALinfo(ORC.name)

## We focus on Ghana
wg.url <- url("http://gsif.isric.org/lib/exe/fetch.php?media=admin.af.rda")
load(wg.url)
proj4string(admin.af) <- "+proj=longlat +datum=WGS84"
country.af <- as(admin.af, "SpatialLines")
## Ghana bounding box:
ghana <- admin.af[admin.af$FORMAL_EN=="Republic of Ghana",]
ghana@bbox

## load soil Africa Soil Profile DB:
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
te = as.vector(ghana@bbox)
unlink("ORC_sl1_Ghana.tif")
system(paste0(gdalwarp, ' ', ORC.name, ' ORC_sl1_Ghana.tif -te ', paste(te, collapse=" ")))
ORCDRC_sl1_ghana <- readGDAL("ORC_sl1_Ghana.tif")
plot(log1p(raster(ORCDRC_sl1_ghana)), col=SAGA_pal[[1]])

## (b) Web Coverage Service
## location of service:
wcs = "http://webservices.isric.org/geoserver/wcs?"
## create an XML file:
l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", "orcdrc_m_sl1_250m", parent=l1)
l1
xml.out = "ORCDRC_M_sl1.xml"
saveXML(l1, file=xml.out)

## check if the layer exists:
system(paste(gdalinfo, xml.out))

## Alternative: calculate offset and region dims:
extent(raster(ORC.name))
bb <- matrix(nrow=2, c(-180,-62.00081,179.9999,87.37))
o.x = 172800 + round(172800*(te[1]-bb[1,2])/(bb[1,2]-bb[1,1]))
o.y = round(71698*(bb[2,2]-te[4])/(bb[2,2]-bb[2,1]))
d.y = round(71698*(te[4]-te[2])/(bb[2,2]-bb[2,1]))
d.x = round(172800*(te[3]-te[1])/(bb[1,2]-bb[1,1]))
o.x; o.y; d.x; d.y
## read only the piece of TIF file to R:
#ORCDRC_sl1_ghana <- readGDAL(ORC.name, offset=c(o.y, o.x), region.dim=c(d.y, d.x), silent=TRUE)
## getCoverage using GDAL translate:
system(paste0(gdal_translate, ' ', xml.out, ' ORC_sl1_Ghana.tif -tr ', 1/120, ' ', 1/120, ' -co \"COMPRESS=DEFLATE\" -srcwin ', paste(c(o.x, o.y, d.x, d.y), collapse=" ")))
## This will only fetch pixels for the bounding box
GDALinfo("ORC_sl1_Ghana.tif")

## plot the map:
ORCDRC_sl1_ghana <- readGDAL("ORC_sl1_Ghana.tif")
data(soil.legends)
class.labels = make.unique(as.character(round((soil.legends[["ORCDRC"]]$MAX-soil.legends[["ORCDRC"]]$MIN)/2, 1)))
ORCDRC_sl1_ghana$val <- cut(ORCDRC_sl1_ghana$band1, breaks=c(soil.legends[["ORCDRC"]]$MIN[1], soil.legends[["ORCDRC"]]$MAX), labels = class.labels)
bnd <- list(list("sp.points", sites, pch="+", col="black"), list("sp.lines", country.af))
spplot(ORCDRC_sl1_ghana["val"], col.regions=soil.legends[["ORCDRC"]]$COLOR, sp.layout=bnd, scales=list(draw=TRUE), main="Organic carbon in permilles (0 cm)")
## Visualization of maps in Google Earth:
plotKML(ORCDRC_sl1_ghana["val"], colour_scale=soil.legends[["ORCDRC"]]$COLOR)

## plot soil types:
t1 <- newXMLNode("WCS_GDAL")
t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
t1.l <- newXMLNode("CoverageName", "taxousda_250m", parent=t1)
saveXML(t1, file="TAXOUSDA.xml")

system(paste0(gdal_translate, ' TAXOUSDA.xml TAXOUSDA_Ghana.tif -tr ', 1/120, ' ', 1/120, ' -r \"nearest\" -co \"COMPRESS=DEFLATE\" -srcwin ', paste(c(o.x, o.y, d.x, d.y), collapse=" ")))
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
var.name <- c(paste0("orcdrc_m_sl", c(1,2,3,4), "_250m"), paste0("bldfie_m_sl", c(1,2,3,4), "_250m"), paste0("crfvol_m_sl", c(1,2,3,4), "_250m"))

## download soil maps in a loop...
for(i in 1:length(var.name)){
 t1 <- newXMLNode("WCS_GDAL")
 t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
 t1.l <- newXMLNode("CoverageName", var.name[i], parent=t1)
 xml.out <- paste(var.name[i], ".xml", sep="")
 saveXML(t1, file=xml.out)
 f.name <- paste(var.name[i], "_Ghana.tif", sep="")
 if(!file.exists(f.name)){
   system(paste0(gdal_translate, ' ', xml.out, ' ', f.name, ' -tr ', 1/120, ' ', 1/120, ' -co \"COMPRESS=DEFLATE\" -srcwin ', paste(c(o.x, o.y, d.x, d.y), collapse=" ")))
 }
}

## read all layers:
library(raster)
ghana1km <- stack(paste(var.name, "_Ghana.tif", sep=""))
ghana1km <- as(as(ghana1km, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
## mask out missing pixels:
sel.NA <- !ghana1km$orcdrc_m_sl1_250m_Ghana < 0 & !ghana1km$crfvol_m_sl1_250m_Ghana < 0 summary(sel.NA)
ghana1km <- ghana1km[sel.NA,]
str(ghana1km@data)

## Standard depths:
std <- c(0,5,15,30)
## Soil Organic Carbon stock formula: http://gsif.r-forge.r-project.org/OCSKGM.html
n.lst=c("orcdrc","bldfie","crfvol")
for(d in c(1,2,3)){
  ## Upper/lower horizons:
  U.lst <- paste0(n.lst, "_m_sl", d, "_250m_Ghana")
  L.lst <- paste0(n.lst, "_m_sl", d+1, "_250m_Ghana")
  ghana1km$ORC <- rowMeans(ghana1km@data[,sapply(c(U.lst[1],L.lst[1]), function(x){grep(x, names(ghana1km))})], na.rm = TRUE)
  ghana1km$BLD <- rowMeans(ghana1km@data[,sapply(c(U.lst[2],L.lst[2]), function(x){grep(x, names(ghana1km))})], na.rm = TRUE)
  ghana1km$CRF <- rowMeans(ghana1km@data[,sapply(c(U.lst[3],L.lst[3]), function(x){grep(x, names(ghana1km))})], na.rm = TRUE)
  ## Predict organic carbon stock (in tones / ha):
  ghana1km@data[,paste0("OCS_sd", d)] <- round(as.vector(OCSKGM(ORCDRC=ghana1km$ORC, BLD=ghana1km$BLD, CRFVOL=ghana1km$CRF, HSIZE=get("stsize", envir = GSIF.opts)[d]*100, ORCDRC.sd=8, BLD.sd=120, CRFVOL.sd=10)*10))
}
#plot(log1p(raster(ghana1km["OCS_sd1"])), col=SAGA_pal[[1]])

## aggregate values and optimize for plotting:
ghana1km$OCS_30cm <- rowSums(ghana1km@data[,paste0("OCS_sd", c(1,2,3))], na.rm=TRUE)

rg <- quantile(ghana1km$OCS_30cm, c(0.01, 0.99), na.rm=TRUE)
at <- expm1(seq(log1p(rg[1]), log1p(rg[2]), length.out=20))
ghana1km$OCS_30cmf <- ifelse(ghana1km$OCS_30cm<rg[1], rg[1], ifelse(ghana1km$OCS_30cm>rg[2], rg[2], ghana1km$OCS_30cm))
spplot(ghana1km["OCS_30cmf"], at=at, col.regions=R_pal[["soc_pal"]], sp.layout=bnd, scales=list(draw=TRUE), main="Total soil carbon stock (0--30 cm) in tonnes per ha")

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

## plot points in Google Earth
kml(pnts, colour=SOURCEID, file="PHIHOX_depth.kml", 
    shape=paste("PHIHOX_depth_", 1:nrow(ov), ".png", sep=""), 
    size=6, points_names=pnts$SOURCEID, 
    colour_scale=rep("#FFFFFF", 2))

## plot soil depth curve:
ORCDRC.pnt1 <- data.frame(
  top=unlist(ov[1,grep("depthCodesMeters.sl", names(ov))])*-100, 
  M=unlist(ov[1,grep("ORCDRC.M", names(ov))]))
ORCDRC.pnt1$variable <- "ORCDRC"
## plot the result:
data(soil.legends)
## Soil organic carbon:
ORCDRC.range = range(soil.legends[["ORCDRC"]]$MIN, soil.legends[["ORCDRC"]]$MAX)
dev.new(width=5, height=6)
xyplot(top ~ M | variable, data=ORCDRC.pnt1, ylab='Depth in cm',
  xlab='SoilGrids predictions', xlim=ORCDRC.range,
  ylim=c(200,0), type="l", lwd=3, alpha=0.25, strip=strip.custom(bg=grey(0.8))
)

##-----------------------------------
## plot texture classes 15-30 cm
##----------------------------------- 

## Bounding box:
te.as <- c(17, -10, 19, -9)
o.x2 = 172800 + round(172800*(te.as[1]-bb[1,2])/(bb[1,2]-bb[1,1]))
o.y2 = round(71698*(bb[2,2]-te.as[4])/(bb[2,2]-bb[2,1]))
d.y2 = round(71698*(te.as[4]-te.as[2])/(bb[2,2]-bb[2,1]))
d.x2 = round(172800*(te.as[3]-te.as[1])/(bb[1,2]-bb[1,1]))
tex.name <- c("clyppt_m_sl3_250m", "sndppt_m_sl3_250m", "sltppt_m_sl3_250m")

## run everything in a loop...
for(i in 1:length(tex.name)){
  t1 <- newXMLNode("WCS_GDAL")
  t1.s <- newXMLNode("ServiceURL", wcs, parent=t1)
  t1.l <- newXMLNode("CoverageName", tex.name[i], parent=t1)
  xml.out <- paste(tex.name[i], ".xml", sep="")
  saveXML(t1, file=xml.out)
  f.name <- paste(tex.name[i], "_angola.tif", sep="")
  if(!file.exists(f.name)){
    system(paste0(gdal_translate, ' ', xml.out, ' ', tex.name[i], '_angola.tif -tr ', 1/120, ' ', 1/120, ' -co \"COMPRESS=DEFLATE\" -srcwin ', paste(c(o.x2, o.y2, d.x2, d.y2), collapse=" ")))
  }
}
angola1km <- stack(paste(tex.name, "_angola.tif", sep=""))
angola1km <- as(angola1km, "SpatialGridDataFrame")
plot(stack(angola1km[1:3]))
## example of texture classification:
TT.points.in.classes(tri.data=data.frame(CLAY=20, SILT=40, SAND=40), class.sys="USDA.TT", PiC.type = "t")

## now run this per row:
names(angola1km@data) <- c("CLAY", "SAND", "SILT")
angola1km$TEX_sl3_M <- TT.points.in.classes(tri.data=angola1km@data[,1:3], class.sys="USDA.TT", PiC.type="t", tri.sum.tst=FALSE)
angola1km$TEX_sl3_M <- as.factor(angola1km$TEX_sl3_M)
summary(angola1km$TEX_sl3_M)
plot(raster(angola1km["TEX_sl3_M"]))
levels(angola1km$TEX_sl3_M)
## colors from http://gsif.isric.org/KML/soil_triangle.png
leg.tex <- c(rgb(240/255,190/255,153/255), rgb(195/255,161/255,136/255), rgb(202/255,195/255,169/255), rgb(195/255,188/255,160/255), rgb(228/255,226/255,203/255), rgb(212/255,190/255,161/255), rgb(209/255,203/255,162/255))
library(plotKML)
kml_open("angola1km_tex.kml")
kml_layer(angola1km, colour=TEX_sd3_M, colour_scale=leg.tex)
kml_screen("http://gsif.isric.org/KML/soil_triangle.png", position = "LR", sname = "triangle")
kml_close("angola1km_tex.kml")

## end of script;