## Prediction of heavy metal concentrations for the Geul river valley in NL;
## Prepared by T. Hengl and G.B.M. Heuvelink

library(GSIF)
library(plotKML)
library(RCurl)
library(sp)
library(rgdal)
library(randomForest)
library(ranger)
library(raster)
library(quantregForest)
library(plotGoogleMaps)
library(leaflet)
nl.rd <- getURL("http://spatialreference.org/ref/sr-org/6781/proj4/")

## Geul data:
geul <- read.table("geul.dat", header = TRUE)
coordinates(geul) <- ~x+y
proj4string(geul) <- CRS(nl.rd) 
grd25 <- readGDAL("dem25.txt")
names(grd25) = "dem"
#grd25$river <- factor(readGDAL("river.txt")$band1, labels = "River")
grd25$dis <- readGDAL("riverdist.txt")$band1
## mask out missing pixels:
grd25 <- as(grd25, "SpatialPixelsDataFrame")
proj4string(grd25) <- CRS(nl.rd) 
# spplot(grd25)

## fit a model:
grd25.spc <- spc(grd25, ~ dem+dis)
m = log1p(pb) ~ PC1+PC2
#pbm <- fit.gstatModel(m, observations=geul, grd25.spc@predicted, family=gaussian(link=log))
pbm <- fit.gstatModel(m, observations=geul, grd25.spc@predicted)
## predict values:
pb.rk0 <- predict(pbm, grd25.spc@predicted)
## block kriging by default!
## back-transform:
pb.rk0@predicted$pb <- expm1(pb.rk0@predicted$pb + pb.rk0@predicted$var1.var/2)
pb.rk0
## 57% of variability explained by the model;
plot(pb.rk0)

## compare to randomForest-kriging:
pb.rk <- autopredict(geul["pb"], grd25)
plot(pb.rk)

## geostat simulations:
pb.rks <- predict(pbm, grd25.spc@predicted, nsim=20, block = c(0,0))
plotKML(pb.rks)

## extend the model using a new covariate:
saga_cmd = shortPathName("C:/SAGA-GIS/saga_cmd.exe")
writeGDAL(grd25["dem"], "dem.sdat", "SAGA", mvFlag=-99999)
system(paste0(saga_cmd, ' ta_hydrology 15 -DEM \"dem.sgrd\" -TWI \"swi.sgrd\"'))
grd25$swi <- readGDAL("swi.sdat")$band1[grd25@grid.index]
plot(stack(grd25))

grd25.spc2 <- spc(grd25, ~ dem+dis+swi)
m2 = pb ~ PC1+PC2+PC3
pbm2 <- fit.gstatModel(m2, observations=geul, grd25.spc2@predicted, method="quantregForest")
pb.rk2 <- predict(pbm2, grd25.spc2@predicted)
show(pb.rk2)
plot(pb.rk2)

## plot in Google Earth:
pb.pol <- grid2poly(pb.rk2@predicted["pb"])
kml(pb.pol, file.name="pb_predicted.kml", colour=pb, colour_scale = SAGA_pal[[1]], kmz=TRUE)
kml(geul, file.name="pb_points.kml", colour=pb, labels=geul$pb)
kml_View("pb_points.kml")
kml_View("pb_predicted.kml")

## plot in GoogleMaps:
str(pb.rk2@predicted)
mp <- plotGoogleMaps(pb.rk2@predicted, filename='geul.html', zcol='pb', add=TRUE, colPalette=SAGA_pal[[1]])
ms <- plotGoogleMaps(geul, filename='geul.html', zcol='pb', add=FALSE, previousMap=mp)

## plot using Leaflet:
r = raster(pb.rk2@predicted["pb"])
pal <- colorNumeric(SAGA_pal[[1]], values(r), na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(r, colors=pal, opacity=0.6) %>%
  addLegend(pal=pal, values=values(r), title="Pb concentration")

## points only:
geul.xy <- spTransform(geul, CRS("+proj=longlat +datum=WGS84"))
leaflet(geul.xy) %>% addTiles() %>%
  addCircleMarkers(radius=geul.xy$pb/50, color="red", stroke=FALSE, fillOpacity=0.8)
