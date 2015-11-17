# title         : tutorial_WorldGrids.R
# purpose       : Accessing and using WorldGrids;
# reference     : [http://worldgrids.org/doku.php?id=wiki:functions]
# producer      : Prepared by T. Hengl and J. Mendes de Jesus
# address       : In Wageningen, NL.
# inputs        : WorldGrids GeoTiffs [http://worldgrids.org]
# outputs       : subsets of WorldGrids;

## load all required libraries
library(XML)
library(sp)
library(rgdal)
library(gdalUtils)
gdal_setInstallation() ## CAN TAKE 1-2 MINUTES!

## layers of interest:
wg.lst <- c("worldgrids:CNTGAD3a", "worldgrids:EVMMOD3a", "worldgrids:G12IGB3a", "worldgrids:TDMMOD3a", "worldgrids:TWISRE3a")

## Web Coverage Service
## location of service:
wcs = "http://wms3.isric.org/geoserver/worldgrids/wcs?"
## create an XML file:
l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", wg.lst[1], parent=l1)
l1
xml.out = paste0(gsub(":", "_", wg.lst[1]), ".xml")
saveXML(l1, file=xml.out)

## test it:
gdalinfo(xml.out)
## works ok!

## subset to some bounding box:
te <- c(20,40,22,42)
bb <- matrix(nrow=2, c(-180,-90,180,90))
o.x = 43200 + round(43200*(te[1]-bb[1,2])/(bb[1,2]-bb[1,1]))
o.y = round(21600*(bb[2,2]-te[4])/(bb[2,2]-bb[2,1]))
d.y = round(21600*(te[4]-te[2])/(bb[2,2]-bb[2,1]))
d.x = round(43200*(te[3]-te[1])/(bb[1,2]-bb[1,1]))
## resample and read to R:
tif.out <- paste0(gsub(":", "_", wg.lst[1]), ".tif")
gdal_translate(src_dataset=xml.out, dst_dataset=tif.out, of="GTiff", srcwin=paste(c(o.x, o.y, d.x, d.y), collapse=" "))
wg1km <- readGDAL(tif.out)
library(raster)
plot(raster(wg1km))

## DO NOT RUN:
library(GSIF)
URI = "http://wps.worldgrids.org/pywps.cgi"
server <- list(URI=URI, request="execute",
    version="version=1.0.0", service.name="service=wps",
    identifier="identifier=sampler_local1pt_nogml")
glcesa3.wps <- new("WPS", server=server, inRastername="glcesa3a")
# show(biocl15.wps)
prl <- getProcess(glcesa3.wps)
prl[7]
describe(glcesa3.wps, identifier="overlay")
p1 <- data.frame(lon=15, lat=15)
coordinates(p1) <- ~lon+lat
proj4string(p1) <- CRS("+proj=longlat +datum=WGS84")
p1
over(glcesa3.wps, p1)
# fetch grids and load the to R:
glcesa3 <- subset(glcesa3.wps, bbox=matrix(c(20,40,22,42), nrow=2))
image(glcesa3)