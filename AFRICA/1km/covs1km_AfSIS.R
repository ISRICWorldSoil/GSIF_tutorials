# title         : covs1km_AfSIS.R
# purpose       : Preparing soil covariate layers for Africa;
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl (tom.hengl@wur.nl)
# address       : In Wageningen, NL.
# inputs        : WorldGrids maps
# outputs       : Rda file with all covariates;
# remarks 1     : 7376 x 7877 pixels (18,394,090 available);

library(rgdal)
library(plotKML)
data(SAGA_pal)

fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
if(is.na(file.info("7za.exe")$size)){
  download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
  unzip("7za920.zip")
  unlink("7za920.zip")
}
si <- Sys.info()
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

######################################################
#  DATA IMPORT AND FILTERING
######################################################

## soil mask map:
afgrid1km <- readGDAL("../soilmasks/SMKMOD3x.tif")
names(afgrid1km) <- "SMKMOD3x"
summary(afgrid1km$SMKMOD3x)
bbox = afgrid1km@bbox
## convert to spatialpixels:
afgrid1km <- as(afgrid1km, "SpatialPixelsDataFrame") 
## takes >10 minutes and >8GB RAM!!!
gc()
str(afgrid1km)

# obtain the 1 km resolution grids:
g.lst <- c("SLPSRT3", "IFLGRE3", "EVMMOD3", "EVSMOD3", "TDMMOD3", "TDSMOD3", "LAMMOD3", "TDHMOD3", "TNMMOD3", "TNSMOD3", "PREGSM1", "LMTGSH3", paste("G0", 1:9, "ESA3", sep=""), paste("G", 10:22, "ESA3", sep=""), "L3POBI3") ## "DEMSRE3" not on the server!
for(j in 1:length(g.lst)){
  outname = paste(g.lst[j], "a.tif.gz", sep="")
  outname.tif = paste(g.lst[j], "a.tif", sep="")
  af_outname.tif = paste("af_", g.lst[j], "a.tif", sep="")
  if(is.na(file.info(af_outname.tif)$size)){
    download.file(paste("http://worldgrids.org/lib/exe/fetch.php?media=", outname, sep=""), outname)
    system(paste("7za e", outname))
    unlink(outname)
  }
  if(is.na(file.info(af_outname.tif)$size)){
    ## resample to the bounding box:
    if(outname.tif %in% c("IFLGRE3a.tif", "LMTGSH3a.tif", "PREGSM1a.tif", "L3POBI3a.tif")){
      system(paste(gdalwarp, outname.tif, af_outname.tif, '-t_srs', paste('\"', af.csy, '\"', sep=""), '-r near -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 1000, 1000))
    }
    if(outname.tif %in% c("TDMMOD3a.tif", "TDHMOD3a.tif", "TNMMOD3a.tif", "EVMMOD3a.tif", "EVSMOD3a.tif")){
      system(paste(gdalwarp, outname.tif, af_outname.tif, '-dstnodata \"-32767\" -t_srs', paste('\"', af.csy, '\"', sep=""), '-r bilinear -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 1000, 1000))
    }
    if(outname.tif %in% c("DEMSRE3a.tif", "SLPSRT3a.tif", "TNSMOD3a.tif", "TDSMOD3a.tif", "LAMMOD3a.tif", paste("G0", 1:9, "ESA3a.tif", sep=""), paste("G", 10:22, "ESA3a.tif", sep=""))){
       system(paste(gdalwarp, outname.tif, af_outname.tif, '-t_srs', paste('\"', af.csy, '\"', sep=""), '-r bilinear -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 1000, 1000))
    }    
  }
  ## remove geotif once done
  unlink(outname.tif)  
}

# import to R (18M pixels per band!!):
tif3.lst <- list.files(pattern=glob2rx("af_*a.tif$"))
for(i in 1:length(tif3.lst)){
    afgrid1km@data[,strsplit(strsplit(tif3.lst[i], ".tif")[[1]][1], "af_")[[1]][2]] <- readGDAL(tif3.lst[i])$band1[afgrid1km@grid.index]
} ## Takes >5 mins!
gc()

## derive SAGA TWI:
#system(paste(gdal_translate, 'af_DEMSRE3a.tif -of \"SAGA\"', 'af_DEMSRE3a.sdat'))
#rsaga.wetness.index("af_DEMSRE3a.sgrd", out.wetness.index="af_TWISRE3a.sgrd")
afgrid1km$TWISRE3a <- readGDAL("af_TWISRE3a.sdat")$band1[afgrid1km@grid.index]

## Africa geology:
#af_geo <- readShapePoly("../../WORLDGRIDS/geology/geo7_2ag.shp", proj4string=CRS("+proj=longlat +datum=WGS84"))
#af_geo$GLG_i <- as.integer(af_geo$GLG)
#af_geo_lead <- spTransform(af_geo["GLG_i"], CRS(af.csy))
#writePolyShape(af_geo_lead, "af_SGEUSG.shp")
#rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="af_SGEUSG3x.sgrd", INPUT="af_SGEUSG.shp", FIELD=1, TARGET=0, LINE_TYPE=1, USER_SIZE=1000, USER_XMIN=bbox[1,1]+500, USER_XMAX=bbox[1,2]-500, USER_YMIN=bbox[2,1]+500, USER_YMAX=bbox[2,2]-500))
# write.table(data.frame(code=1:length(levels(af_geo$GLG)), NAME=levels(af_geo$GLG)), file="af_SGEUSG3x.txt")
afgrid1km$SGEUSG3x <- as.factor(readGDAL("af_SGEUSG3x.sdat")$band1[afgrid1km@grid.index])
xg = read.table("af_SGEUSG3x.txt", header=TRUE)
levels(afgrid1km$SGEUSG3x) <- paste(xg$NAME)
gc()

## filter some images:
afgrid1km$PREGSM1a <- ifelse(afgrid1km$PREGSM1a<0, NA, afgrid1km$PREGSM1a)
xl = read.table("http://worldgrids.org/lib/exe/fetch.php?media=l3pobi.txt", sep="\t", header = TRUE)
afgrid1km$L3POBI3a <- as.factor(afgrid1km$L3POBI3a)
levels(afgrid1km$L3POBI3a) <- c(NA, paste(xl$NAME))
#str(afgrid1km@data)
summary(afgrid1km$L3POBI3a)

save(afgrid1km, file="afgrid1km.rda", compress="xz")
## takes >10 minutes to save!!
#system("xcopy afgrid1km.rda E:\\DropBox\\AFSIS-ISRIC\\SP_predictions\\1km\\ /Y")

# end of script;