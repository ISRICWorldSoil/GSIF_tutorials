# title         : examples_sp_st.R
# purpose       : Demo sp and space-time classes in R
# reference     : See further plotKML package
# producer      : T. Hengl
# address       : In Wageningen
# inputs        : various;
# outputs       : KML files / GIS files;

## Example: Google headquarters in Mountain View, CA
library(sp)
lat = 37.423156
lon = -122.084917
name = "Google headquarters"
pnt = data.frame(name, lat, lon)
coordinates(pnt) <- ~lon+lat
proj4string(pnt) <- CRS("+proj=longlat +datum=WGS84")
pnt
str(pnt, max.level=2)

library(XML)
pnt.kml <- newXMLNode("kml")
h2 <- newXMLNode("Document", parent = pnt.kml)
h3 <- newXMLNode("name", "Google headquarters", parent = h2)
h4 <- newXMLNode("Folder", parent=pnt.kml[["Document"]])
txtc <- sprintf('<Placemark><Point><coordinates>
    %.5f,%.5f,%.0f</coordinates></Point></Placemark>',
    coordinates(pnt)[,1], coordinates(pnt)[,2], rep(0, length(pnt)))
parseXMLAndAdd(txtc, parent = h4)
pnt.kml
saveXML(pnt.kml, "pnt.kml")
library(rgdal)
writeOGR(pnt, "pnt.shp", ".", "ESRI Shapefile")

## Example: creating polygons
lat = c(-122.08634, -122.08648, -122.08120, -122.08119, -122.08634)
lon = c(37.42324, 37.42118, 37.42084, 37.42309, 37.42324)
coords <- as.matrix(data.frame(lat, lon))
# 5 points sorted clock-wise; ending at the begin point:
coords
pol <- Polygons(list(Polygon(coords, hole=as.logical(NA))), ID=name)
pol <- SpatialPolygons(list(pol), proj4string = CRS("+proj=longlat +datum=WGS84"))
pol
str(pol, max.level=2)

pol.kml <- newXMLNode("kml")
h2 <- newXMLNode("Document", parent = pol.kml)
h3 <- newXMLNode("name", "Google headquarters", parent = h2)
h4 <- newXMLNode("Folder", parent=pol.kml[["Document"]])
# prepare coordinates for KML format:
coords_kml <- paste(coords[,1], ',', coords[,2], ',', rep(0, nrow(coords)),
    collapse='\n ', sep = "")
txtc <- sprintf('<Placemark><Polygon><tessellate>1</tessellate>
    <outerBoundaryIs><LinearRing>
    <coordinates>%s</coordinates>
    </LinearRing></outerBoundaryIs></Polygon></Placemark>',  coords_kml)
parseXMLAndAdd(txtc, parent = h4)
pol.kml
saveXML(pol.kml, "pol.kml")

## Example: raster class
library(plotKML)
library(raster)
data(LST)
coordinates(LST) <- ~lon+lat
gridded(LST) <- TRUE
proj4string(LST) <- CRS("+proj=longlat +datum=WGS84")
## write to SAGA GIS:
for(j in 1:ncol(LST)){ writeGDAL(LST[j], paste(names(LST)[j], ".sdat", sep=""), driver="SAGA", mvFlag=-99999) }
## list files:
grd.lst <- list.files(pattern=glob2rx("LST*.sgrd"))
grd.lst
## calculate sd for a stack of rasters (RasterBrick in raster package):
library(RSAGA)
rsaga.get.modules("geostatistics_grid")
rsaga.get.usage("geostatistics_grid", 4)
rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(grd.lst, collapse=";", sep=""), STDDEV="LST_sd.sgrd"))
LST$s.d. <- readGDAL("LST_sd.sdat")$band1[LST@grid.index]
LST@grid@cells.dim
data(SAGA_pal)
kml(LST, colour=s.d., colour_scale=SAGA_pal[[1]], png.width=LST@grid@cells.dim[1]*5, png.height=LST@grid@cells.dim[1]*5)

## the same thing (derive s.d. for rasters) in the the raster package:
LSTb <- brick(LST)
## per layer:
signif(cellStats(LSTb, 'sd'), 3)
## per cell:
LST.sd <- calc(LSTb, sd)
image(LST.sd) ## TH: missing pixels make a problem

## "RasterBrick" is different than SpatialGrid -> it provides much more possibilities
class(LSTb)
slotNames(LSTb)
# get the dates from the file names:
dates <- sapply(strsplit(names(LST), "LST"), function(x){x[[2]]})
LSTb@title = "Time series of MODIS LST (8-day mosaics) images"
LSTb <- setZ(LSTb, format(as.Date(dates, "%Y_%m_%d"), "%Y-%m-%dT%H:%M:%SZ"))
names(LSTb)


## Example: space-time points
data(HRtemp08)
names(HRtemp08)
HRtemp08[1,]

p1 = newXMLNode("Placemark")
# temporal support half day plus/minus:
begin <- format(HRtemp08[1,"DATE"]-.5, "%Y-%m-%dT%H:%M:%SZ")
end <- format(HRtemp08[1,"DATE"]+.5, "%Y-%m-%dT%H:%M:%SZ")
txt <- sprintf('<name>%s</name><TimeStamp><begin>%s</begin>
    <end>%s</end></TimeStamp><Point>
    <coordinates>%.4f,%.4f,%.0f</coordinates></Point>',
    HRtemp08[1,"NAME"], begin, end, HRtemp08[1,"Lon"],
    HRtemp08[1,"Lat"], 0)
parseXMLAndAdd(txt, parent=p1)
p1
saveXML(p1, "p1.kml")

# format the time column:
library(spacetime)
HRtemp08$ctime <- as.POSIXct(HRtemp08$DATE, format="%Y-%m-%dT%H:%M:%SZ")
sp <- SpatialPoints(HRtemp08[,c("Lon","Lat")])
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84")
HRtemp08.st <- STIDF(sp, time = HRtemp08$ctime, data = HRtemp08[,c("NAME","TEMP")])
# subset to first 500 records:
HRtemp08_jan <- HRtemp08.st[1:500]
str(HRtemp08_jan)
# plot in Google Earth: 
plotKML(HRtemp08_jan[,,"TEMP"], dtime = 24*3600, LabelScale = .4)

## Example: trajectories
library(fossil)
library(spacetime)
library(adehabitat)
data(gpxbtour)
# format the time column:
gpxbtour$ctime <- as.POSIXct(gpxbtour$time, format="%Y-%m-%dT%H:%M:%SZ")
coordinates(gpxbtour) <- ~lon+lat
proj4string(gpxbtour) <- CRS("+proj=longlat +datum=WGS84")
xy <- as.list(data.frame(t(coordinates(gpxbtour))))
gpxbtour$dist.km <- sapply(xy, function(x) { 
  deg.dist(long1=x[1], lat1=x[2], long2=xy[[1]][1], lat2=xy[[1]][2]) 
} )
# convert to a STTDF class:
gpx.ltraj <- as.ltraj(coordinates(gpxbtour), gpxbtour$ctime, id = "th")
gpx.st <- as(gpx.ltraj, "STTDF")
gpx.st$speed <- gpxbtour$speed
gpx.st@sp@proj4string <- CRS("+proj=longlat +datum=WGS84")
str(gpx.st)
plotKML(gpx.st, colour="speed", kmz=TRUE)

## Example: foot and mouth disease
library(stpp)
data(fmd)
fmd0  <- data.frame(fmd)
coordinates(fmd0) <- c("X", "Y")
proj4string(fmd0) <- CRS("+init=epsg:27700")
fmd_sp <- as(fmd0, "SpatialPoints")
dates <- as.Date("2001-02-18")+fmd0$ReportedDay
library(spacetime)
fmd_ST <- STIDF(fmd_sp, dates, data.frame(ReportedDay=fmd0$ReportedDay))
data(SAGA_pal)
plotKML(fmd_ST, kmz=TRUE, colour_scale=SAGA_pal[[1]])

# end of script;

