## "Data sharing" examples
## tom.hengl@isric.org

list.of.packages = c("sf", "plotKML", "rgdal", "raster", "RSQLite", "XML", "htmlwidgets", "leaflet", "GSIF", "")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## conversion from raster to table:
library(rgdal)
library(plotKML)
?readGDAL
spnad83 <- readGDAL(system.file("pictures/erdas_spnad83.tif", package = "rgdal")[1])
spnad83.tbl <- as.data.frame(spnad83)
str(spnad83.tbl)

## conversion from vector to table:
library(sf)
data(eberg_zones)
class(eberg_zones)
eberg_zones.tbl <- as(eberg_zones, "sf")
str(eberg_zones.tbl)
write.csv(eberg_zones.tbl, "eberg_zones.csv")
## Simple conversion would not work since any spatial geometry is lost
#eberg_zones.tbl <- as.data.frame(eberg_zones)

## GeoPackage data examples
library(RSQLite)
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
writeOGR(eberg, "eberg.gpkg", "eberg", "GPKG")
## now connect to the DB:
con <- dbConnect(RSQLite::SQLite(), dbname = "eberg.gpkg")
dbListTables(con)
df <- dbGetQuery(con, 'select "soiltype" from eberg')
summary(as.factor(df$soiltype))
dbGetQuery(con, 'select * from gpkg_spatial_ref_sys')[3,"description"]

## points only
eberg_ll = reproject(eberg)
write.csv(as.data.frame(eberg_ll[c("ID","TAXGRSC")]), "eberg_ll.csv")

## connecting to a WCS:
library(XML)
wcs = "http://webservices.isric.org/geoserver/wcs?"
l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", "orcdrc_m_sl1_250m", parent=l1)
l1
xml.out = "ORCDRC_M_sl1.xml"
saveXML(l1, file=xml.out)
system(paste0('gdalinfo ', xml.out))

## spatial query:
spnad83.file = system.file("pictures/erdas_spnad83.tif", package = "rgdal")[1]
spnad83.file
system(paste0('gdallocationinfo ', spnad83.file, ' 100 100'))

## adding metadata:
data("eberg_grid")
gridded(eberg_grid) = ~ x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
writeGDAL(eberg_grid["DEMSRT6"], "eberg_DEM.tif", options="COMPRESS=DEFLATE")
system(paste0('gdal_edit.py -mo \"DESCRIPTION=elevation values from the SRTM DEM\" -mo \"DOWNLOAD_URL=http://geomorphometry.org/content/ebergotzen\" eberg_DEM.tif'))
system('gdalinfo eberg_DEM.tif')

## ploting maps (samples and predictions in the same map):
library(leaflet)
library(htmlwidgets)
library(GSIF)
library(raster)
demo(meuse, echo=FALSE)
omm <- autopredict(meuse["om"], meuse.grid[c("dist","soil","ffreq")], method="ranger", auto.plot=FALSE, rvgm=NULL)
spplot(omm$predicted["om"])
meuse.ll <- reproject(meuse["om"])
m = leaflet() %>% addTiles() %>% addRasterImage(raster(omm$predicted["om"]), colors = SAGA_pal[[1]][4:20]) %>% addCircles(lng = meuse.ll@coords[,1], lat = meuse.ll@coords[,2], color = c('black'), radius=meuse.ll$om)  
saveWidget(m, file="organicmater_predicted.html")
