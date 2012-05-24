
# title         : rw_riodoce.R
# purpose       : preparation of the Rio Doce dataset for DSM;
# reference     : [http://gsif.r-forge.r-project.org/rw_riodoce.R];
# producer      : Prepared by Eliana De Souza & T. Hengl
# last update   : In Wageningen, NL, June 2012.
# inputs        : MS Excel sheet with soil profile data and various GIS data (shape files and GeoTIFFs);
# outputs       : rda files to be used in the GSIF package;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org]; 

# -------------------------------------
# Initial settings
# -------------------------------------

library(rgdal)
library(RSAGA)
library(plyr)
library(aqp)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
si <- Sys.info()
utm.csy <- "+proj=utm +zone=23 +south +ellps=GRS67 +units=m +no_defs"

# -------------------------------------
# Soil profile data
# -------------------------------------

library(gdata)
perl <- gdata:::findPerl("perl")
# perl = "C:/Perl64/bin/perl.exe"
perfis <- read.xls("soil_profile_sad69_20_maio.xls", perl=perl)
plyr:::nunique(perfis$identif)  ## 461 profiles
perfis$SOURCEID <- as.factor(paste("UFV", perfis$identif, sep=""))
# rename columns:
names(perfis)[which(names(perfis)=="Horiz")] <- "HRDLOC" # horizon designation
names(perfis)[which(names(perfis)=="legend")] <- "TAXBRC" # Taxonomy Brasilian Soil classification
names(perfis)[which(names(perfis)=="x")] <- "LONWGS84" # longitude
names(perfis)[which(names(perfis)=="y")] <- "LATWGS84" # latitude
names(perfis)[which(names(perfis)=="Prof_top")] <- "UHDICM"
names(perfis)[which(names(perfis)=="Prof_boton")] <- "LHDICM"
names(perfis)[which(names(perfis)=="PH_H2O")] <- "PHIHO5"
summary(perfis$PHIHO5)
perfis$PHIHO5 <- ifelse(perfis$PHIHO5 < 2, NA, perfis$PHIHO5)
names(perfis)[which(names(perfis)=="PH_KCL")] <- "PHIKCL"
perfis$PHIKCL <- ifelse(perfis$PHIKCL < 2, NA, perfis$PHIKCL)
names(perfis)[which(names(perfis)=="SB")] <- "BSABCL" # Base saturation
perfis$BSABCL <- ifelse(perfis$BSABCL == 0, NA, perfis$BSABCL)
names(perfis)[which(names(perfis)=="CTC_T1")] <- "CEXBCL" # effective Cations Exchange capacity;
perfis$CEXBCL <- ifelse(perfis$CEXBCL == 0, NA, perfis$CEXBCL)
names(perfis)[which(names(perfis)=="MO")] <- "ORCDRC" # Organic carbon in permilles;
perfis$ORCDRC <- ifelse(perfis$ORCDRC == 0, NA, perfis$ORCDRC)
names(perfis)[which(names(perfis)=="Total_sand")] <- "SNDPPT" # total sand;
names(perfis)[which(names(perfis)=="Silt")] <- "SLTPPT" # total silt;
names(perfis)[which(names(perfis)=="clay")] <- "CLYPPT" # total clay;
names(perfis)[which(names(perfis)=="bulk_densi")] <- "BLDVOL" # bulk density;
perfis$BLDVOL <- ifelse(perfis$BLDVOL==0, NA, perfis$BLDVOL)
names(perfis)[which(names(perfis)=="water_capa")] <- "AWAIMM" # water capacity;
perfis$AWAIMM <- ifelse(perfis$AWAIMM==0, NA, perfis$AWAIMM)

depths(perfis) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(perfis) <- ~ LONWGS84 + LATWGS84 + TAXBRC  
# generate SpatialPoints
# coordinates(perfis) <- ~ LONWGS84 + LATWGS84
# assign CRS data
# proj4string(perfis) <- "+proj=latlong +datum=WGS84"
## some can be duplicated:
# perfis@site$SOURCEID[duplicated(perfis@site$SOURCEID)]
str(perfis)
perfis@sp@bbox

riodoce <- list(sites=perfis@site, horizons=perfis@horizons[,c("SOURCEID","HRDLOC","UHDICM","LHDICM","PHIHO5","PHIKCL","ORCDRC","SNDPPT","SLTPPT","CLYPPT","BSABCL","CEXBCL","BLDVOL","AWAIMM")])
str(riodoce)
save(riodoce, file="riodoce.rda", compress="xz")

# -------------------------------------
# Grids
# -------------------------------------

# read the soil map (study area):
mapas <- readShapePoly("mapa_solos_wgs.shp")
bbox <- mapas@bbox
# extend for ca 20 km or 20 pixels:
bbox[,"min"] <- bbox[,"min"]-.2
bbox[,"max"] <- bbox[,"max"]+.2
proj4string(mapas) <- "+proj=longlat +datum=WGS84"
mapas.xy <- spTransform(mapas, CRS(utm.csy))
mapas.xy$legendm <- as.integer(mapas.xy$Legend_RD)
writePolyShape(mapas.xy["legendm"], "mapa_solos_xy.shp")

# read the geology shape:
geol <- readShapePoly("geologia_sad69.shp")
proj4string(geol) <- utm.csy
geol$legendm <- as.integer(geol$Litologia)
writePolyShape(geol["legendm"], "geologia_sad69_xy.shp")

#### 1 km resolutions

# SRTM DEM:
GDALinfo("srtm90m.tif")
system(paste(gdalwarp, "srtm90m.tif -t_srs \"+proj=longlat +datum=WGS84\" DEMSRT3a.sdat -of \"SAGA\" -r bilinear -te", bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], "-tr", 1/120, 1/120))
grid.ll <- readGDAL("DEMSRT3a.sdat")
proj4string(grid.ll) <- "+proj=longlat +datum=WGS84"
grid.xy <- spTransform(grid.ll, CRS(utm.csy))
grid.xy@bbox

# convert soil polygon map to a raster:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="TAXBRC3a_utm.sgrd", GRID_TYPE=0, INPUT="mapa_solos_xy.shp", FIELD=1, TARGET=0, LINE_TYPE=1, USER_SIZE=1000, USER_XMIN=grid.xy@bbox[1,1]+1000/2, USER_XMAX=grid.xy@bbox[1,2]-1000/2, USER_YMIN=grid.xy@bbox[2,1]+1000/2, USER_YMAX=grid.xy@bbox[2,2]-1000/2))
# convert geological map to a raster:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="GEOGLS3a_utm.sgrd", GRID_TYPE=0, INPUT="geologia_sad69_xy.shp", FIELD=1, TARGET=0, LINE_TYPE=1, USER_SIZE=1000, USER_XMIN=grid.xy@bbox[1,1]+1000/2, USER_XMAX=grid.xy@bbox[1,2]-1000/2, USER_YMIN=grid.xy@bbox[2,1]+1000/2, USER_YMAX=grid.xy@bbox[2,2]-1000/2))


# metric system:
system(paste(gdalwarp, "srtm90m.tif -t_srs \"", utm.csy, "\" DEMSRT3a_utm.sdat -of \"SAGA\" -r bilinear -te", grid.xy@bbox[1,1], grid.xy@bbox[2,1], grid.xy@bbox[1,2], grid.xy@bbox[2,2], "-tr", 1e3, 1e3))
gridxy <- readGDAL("DEMSRT3a_utm.sdat")
proj4string(gridxy) <- utm.csy
grid.ll <- spTransform(gridxy, CRS("+proj=longlat +datum=WGS84"))
grid.ll <- as.data.frame(grid.ll)
gridxy$lat <- grid.ll$y
gridxy$long <- grid.ll$x
writeGDAL(gridxy["lat"], "lat_utm.sdat", "SAGA")
writeGDAL(gridxy["long"], "long_utm.sdat", "SAGA")

# Slope map:
rsaga.slope(in.dem="DEMSRT3a_utm.sgrd", out.slope="SLPSRT3a_utm.sgrd")
# Topographic Wetness Index:
rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM="DEMSRT3a_utm.sgrd", C="catharea.sgrd", GN="catchslope.sgrd", CS="modcatharea.sgrd", SB="TWISRT3a_utm.sgrd", T=10))
# valley depth:
rsaga.geoprocessor(lib="ta_morphometry", module=14, param=list(DEM="DEMSRT3a_utm.sgrd", HO="tmp.sgrd", HU="VDPSRT3a_utm.sgrd", NH="tmp.sgrd", SH="tmp.sgrd", MS="tmp.sgrd", W=12, T=120, E=4))
# incoming solar radiation:
## TAKES > 30 mins!
rsaga.geoprocessor(lib="ta_lighting", module=2, param=list(GRD_DEM="DEMSRT3a_utm.sgrd", GRD_DIRECT="tmp.sgrd", GRD_DIFFUS="tmp.sgrd", GRD_TOTAL="INSSRT3a_utm.sgrd", DURATION="durat.sgrd", GRD_LAT="lat_utm.sgrd", GRD_LON="long_utm.sgrd", PERIOD=2, DHOUR=2, DAY_A=1, MON_B=1, DAY_B=31, MON_B=12, DDAYS=5))  # time-consuming!
# convergence index:
rsaga.geoprocessor(lib="ta_morphometry", module=2, param=list(ELEVATION="DEMSRT3a_utm.sgrd", CONVERGENCE="CNVSRT3a_utm.sgrd", RADIUS=3, DISTANCE_WEIGHTING_WEIGHTING=0, SLOPE=FALSE))

## MODIS LST and EVI images (PCA):
mevi.lst <- list.files(path = "E:/MODIS", pattern="EVI", full.names = TRUE)
mlstd.lst <- list.files(path = "E:/MODIS", pattern="LST_Day", full.names = TRUE)
mlstn.lst <- list.files(path = "E:/MODIS", pattern="LST_Night", full.names = TRUE)
# resample to Rio doce grid:
for(i in 1:length(mevi.lst)){
    outname <- set.file.extension(paste("EVI2011", i, sep="_"), ".tif")
    system(paste(gdalwarp, mevi.lst[i], '-t_srs \"', utm.csy, '\" ', outname,' -te', grid.xy@bbox[1,1], grid.xy@bbox[2,1], grid.xy@bbox[1,2], grid.xy@bbox[2,2], '-r bilinear -tr', 1000, 1000))
}
# LST
for(i in 1:length(mlstd.lst)){
    outname <- set.file.extension(paste("LSTD2011", i, sep="_"), ".tif")
    system(paste(gdalwarp, mlstd.lst[i], '-t_srs \"', utm.csy, '\" ', outname,' -te', grid.xy@bbox[1,1], grid.xy@bbox[2,1], grid.xy@bbox[1,2], grid.xy@bbox[2,2], ' -dstnodata \"-32768\" -r bilinear -tr', 1000, 1000))
}
for(i in 1:length(mlstn.lst)){
    outname <- set.file.extension(paste("LSTN2011", i, sep="_"), ".tif")
    system(paste(gdalwarp, mlstn.lst[i], '-t_srs \"', utm.csy, '\" ', outname,' -te', grid.xy@bbox[1,1], grid.xy@bbox[2,1], grid.xy@bbox[1,2], grid.xy@bbox[2,2], ' -dstnodata \"-32768\" -r bilinear -tr', 1000, 1000))
}

## Derive principal components:
mevi.lst <- list.files(pattern="EVI2011")
evigrid <- stack(mevi.lst)
evigrid <- brick(evigrid)
f <- as.formula(paste("~", paste(evigrid@layernames, collapse="+")))
pc.EVI <- prcomp(f, scale=TRUE, as.data.frame(getValues(evigrid)))
evigrid@data@values <- pc.EVI$x
evigrid@layernames <- attr(pc.EVI$x, "dimnames")[2][[1]]
plot(evigrid)
writeRaster(evigrid[["PC1"]], file="EV1MODD.tif", format="GTiff")
writeRaster(evigrid[["PC2"]], file="EV2MODD.tif", format="GTiff")

mlstd.lst <- list.files(pattern="LSTD2011")
lstdgrid <- stack(mlstd.lst)
lstdgrid <- brick(lstdgrid)
# filter the missing values:
x <- scale(getValues(lstdgrid)) 
x[is.na(x)] <- 0 
# PCA:
f2 <- as.formula(paste("~", paste(lstdgrid@layernames, collapse="+")))
pc.LSTD <- prcomp(f2, as.data.frame(x))
lstdgrid@data@values <- pc.LSTD$x
lstdgrid@layernames <- attr(pc.LSTD$x, "dimnames")[2][[1]]
# plot(lstdgrid)
writeRaster(lstdgrid[["PC1"]], file="ST1MODD.tif", format="GTiff")
writeRaster(lstdgrid[["PC2"]], file="ST2MODD.tif", format="GTiff")

mlstn.lst <- list.files(pattern="LSTN2011")
lstngrid <- stack(mlstn.lst)
lstngrid <- brick(lstngrid)
# filter the missing values:
x <- scale(getValues(lstngrid)) 
x[is.na(x)] <- 0 
# PCA:
f2 <- as.formula(paste("~", paste(lstngrid@layernames, collapse="+")))
pc.LSTN <- prcomp(f2, as.data.frame(x))
lstngrid@data@values <- pc.LSTN$x
lstngrid@layernames <- attr(pc.LSTN$x, "dimnames")[2][[1]]
# plot(lstngrid)
writeRaster(lstngrid[["PC1"]], file="ST1MODN.tif", format="GTiff")
writeRaster(lstngrid[["PC2"]], file="ST2MODN.tif", format="GTiff")

# -------------------------------------
# Read grids to R and save as rda
# -------------------------------------

riodoce.grid <- readGDAL("TAXBRC3a_utm.sdat")
names(riodoce.grid) <- "TAXBRC3"
riodoce.grid$GEOGLS3 <- readGDAL("GEOGLS3a_utm.sdat")$band1
riodoce.grid$TAXBRC3 <- ifelse(riodoce.grid$TAXBRC3==129, NA, riodoce.grid$TAXBRC3)
riodoce.grid$GEOGLS3 <- ifelse(riodoce.grid$GEOGLS3==129, NA, riodoce.grid$GEOGLS3)
riodoce.grid$TAXBRC3 <- as.factor(riodoce.grid$TAXBRC3)
riodoce.grid$GEOGLS3 <- as.factor(riodoce.grid$GEOGLS3)
# attach the original names:
levels(riodoce.grid$TAXBRC3) <- levels(mapas$Legend_RD)
# iconv(levels(geol$Litologia), "latin1", to="ASCII", "?") ## problems with non-standard names
levels(riodoce.grid$GEOGLS3) <- c("AGL","ANF","ARC","ARM","DIO","END","FER","GNA","GRA","KIN","MTB","MTR","QUR","ROC","ROP","SOL","TON","TUR")
riodoce.grid$DEMSRT3 <- round(readGDAL("DEMSRT3a_utm.sdat")$band1, 0)
riodoce.grid$SLPSRT3 <- round(readGDAL("SLPSRT3a_utm.sdat")$band1, 2)
riodoce.grid$TWISRT3 <- round(readGDAL("TWISRT3a_utm.sdat")$band1, 1)
# riodoce.grid$VDPSRT3 <- round(readGDAL("VDPSRT3a_utm.sdat")$band1, 0)
riodoce.grid$INSSRT3 <- round(readGDAL("INSSRT3a_utm.sdat")$band1, 0)
# riodoce.grid$CNVSRT3 <- round(readGDAL("CNVSRT3a_utm.sdat")$band1, 1)
riodoce.grid$EV1MOD3 <- round(readGDAL("EV1MODD.tif")$band1, 1)
riodoce.grid$EV2MOD3 <- round(readGDAL("EV2MODD.tif")$band1, 1)
riodoce.grid$TD1MOD3 <- round(readGDAL("ST1MODD.tif")$band1, 1)
# riodoce.grid$TD2MOD3 <- round(readGDAL("ST2MODD.tif")$band1, 1)
riodoce.grid$TN1MOD3a <- round(readGDAL("ST1MODN.tif")$band1, 1)
# riodoce.grid$TN2MOD3 <- round(readGDAL("ST2MODN.tif")$band1, 1)
riodoce.grid@bbox <- gridxy@bbox
# spplot(riodoce.grid[8], col.regions=bpy.colors(30))
# spplot(riodoce.grid[2], col.regions=bpy.colors(18))

riodoce.grids <- as.data.frame(riodoce.grid)
riodoce.grids <- riodoce.grids[!is.na(riodoce.grids$TAXBRC3),] 
str(riodoce.grids)
save(riodoce.grids, file="riodoce.grids.rda", compress="xz")


# end of script;