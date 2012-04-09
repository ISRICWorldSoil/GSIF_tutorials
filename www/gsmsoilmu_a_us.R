# title         : gsmsoilmu_a_us.R
# purpose       : preparing U.S. General Soil Map data for DSM;
# reference     : http://soildatamart.sc.egov.usda.gov
# producer      : Prepared by T. Hengl and Dylan Beaudette
# last update   : In Wageningen, NL, March 2012.
# inputs        : geotiff with dominant soil orders, sub-groups; 
# outputs       : Geographic Coordinate System (NAD83), original file is about 380 MB;

# ------------------------------------------------------------
# Initial settings and data loading
# ------------------------------------------------------------

require(aqp)
require(gstat)
require(rgdal)
fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
gdalwarp <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
workd <- normalizePath(getwd())
aea.prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# bounding box whole of USA:
minX = -2405000
maxX = 2295000
minY = 260000
maxY = 3220000
pixsize = 5000

## download the General soil map of USA:
# download.file("http://soildatamart.sc.egov.usda.gov/download/export/e_1819206/gsmsoil_us.zip", "gsmsoil_us.zip")
download.file("http://globalsoilmap.net/data/gsmsoil_us.zip", "gsmsoil_us.zip") # 347 MB!!
unzip("gsmsoil_us.zip")

# ------------------------------------------------------------
# Derive dominant soil taxa
# ------------------------------------------------------------
 
# reproject shape:
rsaga.geoprocessor(lib="pj_proj4", module=5, param=list(SOURCE="gsmsoil_us/spatial/gsmsoilmu_a_us.shp", TARGET="gsmsoilmu_a_us_aea.shp", SOURCE_PROJ=paste('\"', "+proj=latlong +datum=NAD83", '\"', sep=""), TARGET_PROJ=paste('\"', aea.prj, '\"', sep="")))
# gsmsoil <- readShapePoly("gsmsoilmu_a_us_aea.shp")
# str(gsmsoil@data)
# THE DATABASE TABLE CAN NOT BE READ TO R! fix the DBF file manually:
dbf <- foreign::read.dbf("gsmsoilmu_a_us_aea.dbf", as.is=TRUE)
dbf$MUKEY <- as.numeric(sapply(dbf$MUSYM, function(x){substr(x, start=nchar(x)-5, stop=nchar(x))}))
dbf$MUSYM <- sapply(dbf$MUSYM, function(x){substr(x, start=1, stop=nchar(x)-6)})
unlink("gsmsoilmu_a_us_aea.dbf")
foreign::write.dbf(dbf, file="gsmsoilmu_a_us_aea.dbf")

# read the table data (taxonomy):
comp <- read.table("gsmsoil_us/tabular/comp.txt", sep="|")
dim(comp)
str(comp) # V85 = "suborder",  # V108 = MUKEY 
comp$V108 <- as.factor(comp$V108)
length(levels(comp$V108))  
length(levels(as.factor(comp$V85)))  # 69
length(levels(as.factor(comp$V86)))  # 251
# re-name for simplicity
names(comp)[2] <- 'comppct'
names(comp)[85] <- 'suborder'
names(comp)[86] <- 'grtgroup'
names(comp)[108] <- 'mukey'


# ------------------------------------------------------------
# Illinois and Indiana only
# ------------------------------------------------------------

# rasterize:
pixsize = 2000
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(INPUT="gsmsoilmu_a_us_aea.shp", FIELD=3, MULTIPLE=1, LINE_TYPE=0, GRID_TYPE=2, TARGET=0, USER_XMIN=379466.3+pixsize/2, USER_YMIN=1569558.9+pixsize/2, USER_XMAX=960405.3-pixsize/2, USER_YMAX=2192919.3-pixsize/2, USER_SIZE=pixsize, USER_GRID="gsmsoil_mukey_ilin.sgrd"))
gsmsoil.ilin <- readGDAL("gsmsoil_mukey_ilin.sdat")
names(gsmsoil.ilin) <- c("mukey")
gsmsoil.ilin$mukey <- as.factor(gsmsoil.ilin$mukey)
length(levels(gsmsoil.ilin$mukey))

# function for computing total component percentages
f.sum <- function(i) {
	n <- nrow(i)
	s <- sum(i$comppct)
	return(data.frame(pct=s, n=n))
}

# tabulate percentage of suborder-level taxa within each map unit (mukey)
comp.suborder.sums <- ddply(comp, .(mukey, suborder), f.sum, .progress='text')
# takes 3-4 mins
str(comp.suborder.sums)

# function for picking the largest suborder from within each map unit (mukey)
f.largest <- function(i) {
	i.sorted <- i[order(i$pct, decreasing=TRUE), ]
	top.suborder <- i.sorted$suborder[1]
	top.suborder.pct <- i.sorted$pct[1]
	return(data.frame(suborder=top.suborder, pct=top.suborder.pct))
}

# keep the largest suborder, and its associated total percentage
comp.suborder <- ddply(comp.suborder.sums, .(mukey), f.largest,
.progress='text')

## check, notice that picking the most common suborder from the raw data gives the wrong result:
# comp.check <- comp[comp$mukey == 660908, c('comppct', 'suborder', 'mukey')]
# comp.check[order(comp.check$comppct, decreasing=TRUE), ]

# we get the correct answer from the data that has been tabulated by suborder within each map unit key
comp.suborder[comp.suborder$mukey == 660908, ]
# ... the difference between orthents and xererts is huge, however with nearly equal percentages picking either one will result in a poor description of the map unit as a whole

## important notes:
# 1. there are cases where the 'largest' | 'dominant' | 'most common' suborder is < 30% of the map unit.... which means that picking a single suborer is not likely all that informative.
# 2. this process is much more efficient in SQL via RDBMS

# merge dominant suborder and the coordiantes:
taxsuborder.ilin <- merge(comp.suborder, as.data.frame(gsmsoil.ilin), by="mukey", all.y=TRUE)
coordinates(taxsuborder.ilin) <- ~x+y
gridded(taxsuborder.ilin) <- TRUE
fullgrid(taxsuborder.ilin) <- TRUE
# to write to SAGA GIF format requires indicators
taxsuborder.ilin$suborder.i <- as.integer(taxsuborder.ilin$suborder)
proj4string(taxsuborder.ilin) <- CRS(aea.prj)
writeGDAL(taxsuborder.ilin["suborder.i"], "taxsuborder.tif", "GTiff", type = "Byte", mvFlag = 255)
writeGDAL(taxsuborder.ilin["pct"], "taxsuborder_pct.sdat", "SAGA", type = "Byte", mvFlag = 255)
write.table(data.frame(code=1:length(levels(taxsuborder.ilin$suborder)), levels=levels(taxsuborder.ilin$suborder)), "taxsuborder.txt")

# ------------------------------------------------------------
# Whole of USA
# ------------------------------------------------------------

# rasterize:
pixsize = 5000
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(INPUT="gsmsoilmu_a_us_aea.shp", FIELD=3, MULTIPLE=1, LINE_TYPE=0, GRID_TYPE=2, TARGET=0, USER_XMIN=minX+pixsize/2, USER_YMIN=minY+pixsize/2, USER_XMAX=maxX-pixsize/2, USER_YMAX=maxY-pixsize/2, USER_SIZE=pixsize, USER_GRID="gsmsoil_mukey.sgrd"))
gsmsoil <- readGDAL("gsmsoil_mukey.sdat")
names(gsmsoil) <- c("mukey")
gsmsoil$mukey <- as.factor(gsmsoil$mukey)
length(levels(gsmsoil$mukey))
gsmsoil <- merge(comp.suborder, as.data.frame(gsmsoil), by="mukey", all.y=TRUE)
coordinates(gsmsoil) <- ~x+y
gridded(gsmsoil) <- TRUE
fullgrid(gsmsoil) <- TRUE
# to write to SAGA GIF format requires indicators
gsmsoil$suborder.i <- as.integer(gsmsoil$suborder)
proj4string(gsmsoil) <- CRS(aea.prj)
writeGDAL(gsmsoil["suborder.i"], "taxsuborder.tif", "GTiff", type = "Byte", mvFlag = 255)
writeGDAL(gsmsoil["pct"], "taxsuborder_pct.tif", "GTiff", type = "Byte", mvFlag = 255)
write.table(data.frame(code=1:length(levels(gsmsoil$suborder)), levels=levels(gsmsoil$suborder)), "taxsuborder.txt")

# ------------------------------------------------------------
# Add extra soil predictors from Worldmaps (5 km)
# ------------------------------------------------------------

URL <- "http://spatial-analyst.net/worldmaps/"
map.list <- c("biocl5", "biocl6", "PRECm", "biocl15", "globedem", "slope", "globcov", "LAIm", "PCEVI1", "PCEVI2", "PCEVI3", "PCEVI4", "wildness", "LSTDm", "LSTDs", "LSTNm", "LSTNs", "anthroms")
# download the zipped maps one by one and resample to 5 km resolution:
for(i in 1:length(map.list)) {
  if(is.na(file.info(paste(map.list[i], "_us.tif", sep=""))$size)){
  download.file(paste(URL, map.list[i], ".zip", sep=""), destfile=paste(getwd(), "/", map.list[i], ".zip", sep=""))
  unzip(paste(getwd(), "/", map.list[i], ".zip", sep=""))
  # resample to USA grid:
  if(map.list[i]=="globcov"|map.list[i]=="anthroms"){
  system(paste(gdalwarp, ' ', set.file.extension(map.list[i], ".tif"), ' -t_srs \"', aea.prj, '\" ', paste(map.list[i], "_us.tif", sep=""),' -te ', gsmsoil@bbox[1,1], ' ', gsmsoil@bbox[2,1], ' ', gsmsoil@bbox[1,2], ' ', gsmsoil@bbox[2,2], ' -r near -tr ', pixsize, ' ', pixsize, sep=""))
  }
  else {
  system(paste(gdalwarp, ' ', set.file.extension(map.list[i], ".tif"), ' -t_srs \"', aea.prj, '\" ', paste(map.list[i], "_us.tif", sep=""),' -te ', gsmsoil@bbox[1,1], ' ', gsmsoil@bbox[2,1], ' ', gsmsoil@bbox[1,2], ' ', gsmsoil@bbox[2,2], ' -r bilinear -tr ', pixsize, ' ', pixsize, sep=""))
  }
  # Delete temporary file:
  unlink(paste(map.list[i], ".zip", sep=""))
  unlink(paste(map.list[i], ".tif", sep=""))  
  }
}

tif.lst <- list.files(pattern="*_us.tif")
# read all resampled maps to R:
for(i in 1:length(tif.lst)){
  gsmsoil@data[,strsplit(tif.lst[i], "_")[[1]][1]] <- readGDAL(tif.lst[i])$band1
}
str(gsmsoil@data)
gsmsoil$globcov <- as.factor(gsmsoil$globcov)
gsmsoil$anthroms <- as.factor(gsmsoil$anthroms)
# recyclye some columns not needed any more:
# gsmsoil$mukey <- NULL
gsmsoil$suborder.i <- NULL


# ------------------------------------------------------------
# Reformat and export maps
# ------------------------------------------------------------

USgrids <- as.data.frame(gsmsoil)
USgrids <- USgrids[!is.na(USgrids$suborder),]
str(USgrids)
# 312,366 pixels
# round up the numbers:
for(j in c("LSTDm", "LSTDs", "LSTNm", "LSTNs", "PCEVI1", "PCEVI2", "PCEVI3", "PCEVI4")){
  USgrids[,j] <- round(USgrids[,j], 1)
}
save(USgrids, file="USgrids.rda", compress="xz")

# end of script;