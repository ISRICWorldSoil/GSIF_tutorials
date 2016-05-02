## Some examples of how to prepare soil covariates for spatial prediction
## By: tom.hengl@isric.org

library(plotKML)
library(sp)
library(gstat)
library(raster)
?eberg
data(eberg_zones)
spplot(eberg_zones[1])

## Rasterize using 'raster' package:
data("eberg_grid25")
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")
r <- raster(eberg_grid25)
r
eberg_zones_r <- rasterize(eberg_zones, r, field="ZONES")
plot(eberg_zones_r)

## Rasterize using SAGA GIS:
library(rgdal)
writeOGR(eberg_zones, "eberg_zones.shp", ".", "ESRI Shapefile")
saga_cmd = "C:/SAGA-GIS/saga_cmd.exe"
system(paste0(saga_cmd, ' grid_gridding 0 -INPUT=\"eberg_zones.shp\" -FIELD=\"ZONES\" -OUTPUT=\"eberg_zones.sgrd\" -GRID_TYPE=0 -TARGET_DEFINITION=1 -TARGET_USER_SIZE=25 -TARGET_USER_XMIN=', extent(r)[1],' -TARGET_USER_XMAX=', extent(r)[2], ' -TARGET_USER_YMIN=', extent(r)[3],' -TARGET_USER_YMAX=', extent(r)[4]))

