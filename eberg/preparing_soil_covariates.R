## Some examples of how to prepare soil covariates for spatial prediction
## By: tom.hengl@isric.org and bas.kempen@wur.nl

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
## TAKES TIME! Raster package is not maybe most suitable for larger data sets
plot(eberg_zones_r)

## Rasterize using SAGA GIS:
library(rgdal)
eberg_zones$ZONES_int <- as.integer(eberg_zones$ZONES)
writeOGR(eberg_zones["ZONES_int"], "eberg_zones.shp", ".", "ESRI Shapefile")
## specify location of SAGA GIS:
saga_cmd = "C:/Progra~1/SAGA-GIS/saga_cmd.exe"
pix = 25
system(paste0(saga_cmd, ' grid_gridding 0 -INPUT \"eberg_zones.shp\" -FIELD \"ZONES_int\" -GRID \"eberg_zones.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', pix, ' -TARGET_USER_XMIN ', extent(r)[1]+pix/2,' -TARGET_USER_XMAX ', extent(r)[2]-pix/2, ' -TARGET_USER_YMIN ', extent(r)[3]+pix/2,' -TARGET_USER_YMAX ', extent(r)[4]-pix/2))
eberg_zones_r2 <- readGDAL("eberg_zones.sdat")
levels(eberg_zones$ZONES)
eberg_zones_r2$ZONES <- as.factor(eberg_zones_r2$band1)
levels(eberg_zones_r2$ZONES) = levels(eberg_zones$ZONES)
summary(eberg_zones_r2$ZONES)
spplot(eberg_zones_r2["ZONES"])

## Downscale:
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
names(eberg_grid)

writeGDAL(eberg_grid["TWISRT6"], "eberg_grid_TWISRT6.tif")
gdalwarp = "C:/Progra~1/GDAL/gdalwarp.exe"
system(paste0(gdalwarp, ' eberg_grid_TWISRT6.tif eberg_grid_TWISRT6_25m.tif -r \"cubicspline\" -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "),' -tr ', pix, ' ', pix, ' -overwrite'))
## compare maps:
par(mfrow=c(1,2))
zlim = range(eberg_grid$TWISRT6, na.rm=TRUE)
image(raster(eberg_grid["TWISRT6"]), col=SAGA_pal[[1]], zlim=zlim, main="Original", asp=1)
image(raster("eberg_grid_TWISRT6_25m.tif"), col=SAGA_pal[[1]], zlim=zlim, main="Downscaled", asp=1)
## Aggregate:
system(paste0(gdalwarp, ' eberg_grid_TWISRT6.tif eberg_grid_TWISRT6_250m.tif -r \"average\" -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "),' -tr 250 250 -overwrite'))
par(mfrow=c(1,2))
image(raster(eberg_grid["TWISRT6"]), col=SAGA_pal[[1]], zlim=zlim, main="Original", asp=1)
image(raster("eberg_grid_TWISRT6_250m.tif"), col=SAGA_pal[[1]], zlim=zlim, main="Aggregated", asp=1)
grid250m <- readGDAL("eberg_grid_TWISRT6_250m.tif")
grid250m <- grid2poly(grid250m)
## plot in Google Earth
kml(grid250m, colour=band1, colour_scale=SAGA_pal[[1]], kmz=TRUE)

## Derive some standard DEM variables of interest for soil mapping:
saga_DEM_derivatives <- function(INPUT, MASK=NULL, sel=c("SLP","TWI","CRV","VBF","VDP","OPN","DVM")){
  if(!is.null(MASK)){
    ## Fill in missing DEM pixels:
    suppressWarnings( system(paste0(saga_cmd, ' grid_tools 25 -GRID=\"', INPUT, '\" -MASK=\"', MASK, '\" -CLOSED=\"', INPUT, '\"')) )
  }
  ## Slope:
  if(any(sel %in% "SLP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 0 -ELEVATION=\"', INPUT, '\" -SLOPE=\"', gsub(".sgrd", "_slope.sgrd", INPUT), '\" -C_PROF=\"', gsub(".sgrd", "_cprof.sgrd", INPUT), '\"') ) ) )
  }
  ## TWI:
  if(any(sel %in% "TWI")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_hydrology 15 -DEM=\"', INPUT, '\" -TWI=\"', gsub(".sgrd", "_twi.sgrd", INPUT), '\"') ) ) )
  }
  ## MrVBF:
  if(any(sel %in% "VBF")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 8 -DEM=\"', INPUT, '\" -MRVBF=\"', gsub(".sgrd", "_vbf.sgrd", INPUT), '\" -T_SLOPE=10 -P_SLOPE=3') ) ) )
  }
  ## Valley depth:
  if(any(sel %in% "VDP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_channels 7 -ELEVATION=\"', INPUT, '\" -VALLEY_DEPTH=\"', gsub(".sgrd", "_vdepth.sgrd", INPUT), '\"') ) ) )
  }
  ## Openess:
  if(any(sel %in% "OPN")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_lighting 5 -DEM=\"', INPUT, '\" -POS=\"', gsub(".sgrd", "_openp.sgrd", INPUT), '\" -NEG=\"', gsub(".sgrd", "_openn.sgrd", INPUT), '\" -METHOD=0' ) ) ) )
  }
  ## Deviation from Mean Value:
  if(any(sel %in% "DVM")){
    suppressWarnings( system(paste0(saga_cmd, ' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean.sgrd", INPUT), '\" -RADIUS=11' ) ) )
  }
}

## test it and plot all results:
writeGDAL(eberg_grid["DEMSRT6"], "DEMSRT6.sdat", "SAGA")
saga_DEM_derivatives("DEMSRT6.sgrd")
dem.lst <- list.files(pattern=glob2rx("^DEMSRT6_*.sdat"))
plot(stack(dem.lst), col=SAGA_pal[[1]])

## Gap filling:
eberg_grid$test <- eberg_grid$TWISRT6
eberg_grid@data[sample.int(nrow(eberg_grid), 1000),"test"] <- NA
writeGDAL(eberg_grid["test"], "test.sdat", "SAGA")
system(paste0(saga_cmd, ' grid_tools 7 -INPUT \"test.sgrd\" -RESULT \"test.sgrd\"'))
par(mfrow=c(1,2))
image(raster(eberg_grid["test"]), col=SAGA_pal[[1]], zlim=zlim, main="Original", asp=1)
image(raster("test.sdat"), col=SAGA_pal[[1]], zlim=zlim, main="Filtered", asp=1)

## SPC:
library(GSIF)
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
eberg_spc <- spc(eberg_grid, formulaString)
names(eberg_spc@predicted) # 11 components on the end;
## Not run: # plot maps:
rd = range(eberg_spc@predicted@data[,1], na.rm=TRUE)
plot(stack(eberg_spc@predicted[1:11]), zlim=rd, col=rev(rainbow(65)[1:48]))

## Read files to R using 'stack':
writeGDAL(eberg_grid25["DEMTOPx"], "DEMTOPx_25m.tif")
writeGDAL(eberg_grid25["NVILANx"], "NVILANx_25m.tif")
writeGDAL(eberg_zones_r2[1], "eberg_zones_25m.tif")
grd.lst <- list.files(pattern="25m")
grd.lst
grid25m <- stack(grd.lst)
grid25m <- as(grid25m, "SpatialGridDataFrame")
str(grid25m)

### save covariate stack to RDATA file for future use
save(grid25m, file = "covariates25m.RDATA")

## Overlay using 'extract'
library(sp)
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
ov = as.data.frame(extract(stack(grd.lst), eberg))
str(ov[complete.cases(ov),])

## Overlay is implemented in the GSIF package by default:
demo(meuse, echo=FALSE)
## simple model:
omm <- fit.gstatModel(meuse, om~dist+ffreq, meuse.grid, family = gaussian(log))
