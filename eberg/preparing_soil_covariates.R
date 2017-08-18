## Some examples of how to prepare soil covariates for spatial prediction
## http://gsif.isric.org/doku.php/wiki:soil_covariates#working_with_larger_rasters
## By: tom.hengl@isric.org and bas.kempen@wur.nl

library(plotKML)
library(rgdal)
library(gstat)
library(raster)
?eberg
data(eberg_zones)
spplot(eberg_zones)
#plotKML(eberg_zones)

## Rasterize using 'raster' package ----
data("eberg_grid25")
gridded(eberg_grid25) <- ~x+y
proj4string(eberg_grid25) <- CRS("+init=epsg:31467")
r <- raster(eberg_grid25)
r
system.time( eberg_zones_r <- rasterize(eberg_zones, r, field="ZONES") )
## TAKES TIME! Raster package is not maybe most suitable for larger data sets
plot(eberg_zones_r)

## Rasterize using SAGA GIS ----
library(rgdal)
eberg_zones$ZONES_int <- as.integer(eberg_zones$ZONES)
writeOGR(eberg_zones["ZONES_int"], "eberg_zones.shp", ".", "ESRI Shapefile")
## specify location of SAGA GIS:
if(.Platform$OS.type=="windows"){
  saga_cmd = shortPathName("C:/SAGA-GIS/saga_cmd.exe")
} else {
  saga_cmd = "saga_cmd"
}
pix = 25
## http://saga-gis.org/saga_module_doc/2.2.7/grid_gridding_0.html
system(paste0(saga_cmd, ' grid_gridding 0 -INPUT \"eberg_zones.shp\" -FIELD \"ZONES_int\" -GRID \"eberg_zones.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', pix, ' -TARGET_USER_XMIN ', extent(r)[1]+pix/2,' -TARGET_USER_XMAX ', extent(r)[2]-pix/2, ' -TARGET_USER_YMIN ', extent(r)[3]+pix/2,' -TARGET_USER_YMAX ', extent(r)[4]-pix/2))
eberg_zones_r2 <- readGDAL("eberg_zones.sdat")
levels(eberg_zones$ZONES)
eberg_zones_r2$ZONES <- as.factor(eberg_zones_r2$band1)
levels(eberg_zones_r2$ZONES) = levels(eberg_zones$ZONES)
summary(eberg_zones_r2$ZONES)
spplot(eberg_zones_r2["ZONES"])

## Downscale ----
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
## Aggregate ----
system(paste0(gdalwarp, ' eberg_grid_TWISRT6.tif eberg_grid_TWISRT6_250m.tif -r \"average\" -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "),' -tr 250 250 -overwrite'))
par(mfrow=c(1,2))
image(raster(eberg_grid["TWISRT6"]), col=SAGA_pal[[1]], zlim=zlim, main="Original", asp=1)
image(raster("eberg_grid_TWISRT6_250m.tif"), col=SAGA_pal[[1]], zlim=zlim, main="Aggregated", asp=1)
grid250m <- readGDAL("eberg_grid_TWISRT6_250m.tif")
grid250m <- grid2poly(grid250m)
## plot in Google Earth
kml(grid250m, colour=band1, colour_scale=SAGA_pal[[1]], kmz=TRUE)

## Derive DEM variables of interest for soil mapping ----
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
r = readGDAL(dem.lst[6])
plotKML(r, colour_scale=SAGA_pal[[1]])

## Gap filling ----
eberg_grid$test <- eberg_grid$TWISRT6
eberg_grid@data[sample.int(nrow(eberg_grid), 1000),"test"] <- NA
writeGDAL(eberg_grid["test"], "test.sdat", "SAGA")
system(paste0(saga_cmd, ' grid_tools 7 -INPUT \"test.sgrd\" -RESULT \"test.sgrd\"'))
par(mfrow=c(1,2))
image(raster(eberg_grid["test"]), col=SAGA_pal[[1]], zlim=zlim, main="Original", asp=1)
image(raster("test.sdat"), col=SAGA_pal[[1]], zlim=zlim, main="Filtered", asp=1)

## SPC ----
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

## Read files to R using 'stack' ----
writeGDAL(eberg_grid25["DEMTOPx"], "DEMTOPx_25mf.tif", options=c("COMPRESS=DEFLATE"))
writeGDAL(eberg_grid25["NVILANx"], "NVILANx_25mf.tif", options=c("COMPRESS=DEFLATE"))
writeGDAL(eberg_zones_r2[1], "eberg_zones_25mf.tif", options=c("COMPRESS=DEFLATE"))
grd.lst <- list.files(pattern="_25mf.tif")
grd.lst
grid25m <- stack(grd.lst)
grid25m <- as(grid25m, "SpatialGridDataFrame")
str(grid25m)

### save covariate stack to RDATA file for future use
save(grid25m, file = "covariates25m.RDATA")

## Overlay using 'extract' ----
library(sp)
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
ov = as.data.frame(raster::extract(stack(grd.lst), eberg, na.rm=FALSE))
str(ov[complete.cases(ov),])

## If the raster layers can not be stacked and if each layer is available in a different projection system, we can make a more generic function:
overlay.fun = function(i, y){raster::extract(raster(i), na.rm=FALSE, spTransform(y, proj4string(raster(i))))}
library(parallel)
ov = data.frame(mclapply(grd.lst, FUN=overlay.fun, y=eberg))
names(ov) = basename(grd.lst)
str(ov[complete.cases(ov),])

## Overlay is implemented in the GSIF package by default:
demo(meuse, echo=FALSE)
## simple model:
omm <- fit.gstatModel(meuse, om~dist+ffreq, meuse.grid, family = gaussian(log))
plot(omm)

## Tiling - working with larger rasters ----
fn = system.file("pictures/SP27GTIF.TIF", package = "rgdal")
obj <- GDALinfo(fn)
tiles <- GSIF::getSpatialTiles(obj, block.x=2000, return.SpatialPolygons = FALSE)
tiles.pol <- GSIF::getSpatialTiles(obj, block.x=2000, return.SpatialPolygons = TRUE)
tile.pol = SpatialPolygonsDataFrame(tiles.pol, tiles)
plot(raster(fn), col=bpy.colors(20))
lines(tile.pol, lwd=2)

## function to parallelize:
fun_mask <- function(i, tiles, dir="./tiled/", threshold=190){
  out.tif = paste0(dir, "T", i, ".tif")
  if(!file.exists(out.tif)){
    x = readGDAL(fn, offset=unlist(tiles[i,c("offset.y","offset.x")]), region.dim=unlist(tiles[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tiles[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    x$mask = ifelse(x$band1>threshold, 1, 0)
    writeGDAL(x["mask"], type="Byte", mvFlag = 255, out.tif, options=c("COMPRESS=DEFLATE"))
  }
}

## Run in parallel:
dir.create("./tiled")
x0 = mclapply(1:nrow(tiles), FUN=fun_mask, tiles=tiles)
## Mosaick back results of computing:
t.lst <- list.files(path="./tiled", pattern=glob2rx("^T*.tif$"), full.names=TRUE, recursive=TRUE)
cat(t.lst, sep="\n", file="SP27GTIF_tiles.txt")
system('gdalbuildvrt -input_file_list SP27GTIF_tiles.txt SP27GTIF.vrt')
system('gdalwarp SP27GTIF.vrt SP27GTIF_mask.tif -ot \"Byte\" -dstnodata 255 -co \"BIGTIFF=YES\" -r \"near\" -overwrite -co \"COMPRESS=DEFLATE\"')
plot(raster("SP27GTIF_mask.tif"))

## Parallization of ploting of the map from the plotKML package ----
library(parallel)
library(snowfall)
plotKML.GDALobj(obj, tiles=tiles, z.lim=c(0,185))
