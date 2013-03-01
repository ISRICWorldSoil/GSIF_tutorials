# title         : covs20km_AfSIS.R
# purpose       : Preparing soil covariate layers for Africa;
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Dec 9 2012.
# inputs        : WorldGrids maps
# outputs       : Rda file with all covariates;
# remarks 1     : ;

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
afgrid20 <- readGDAL("../soilmasks/SMKMOD0x.tif")
names(afgrid20) <- "SMKMOD0x"
bbox = afgrid20@bbox
# obtain the 20 km resolution grids:
g.lst <- c("DEMSRE0", "SLPSRT0", "IFLGRE0", "EVMMOD0", "EVSMOD0", "TDMMOD0", "TDSMOD0", "TDHMOD0", "TNMMOD0", "TNSMOD0", "PREGSM0", "LMTGSH0", paste("G0", 1:9, "ESA0", sep=""), paste("G", 10:22, "ESA0", sep=""), "L3POBI0")
for(j in 1:length(g.lst)){
  outname = paste(g.lst[j], "a.tif.gz", sep="")
  outname.tif = paste(g.lst[j], "a.tif", sep="")
  af_outname.tif = paste("af_", g.lst[j], "a.tif", sep="")
  if(is.na(file.info(outname.tif)$size)){
    download.file(paste("http://worldgrids.org/lib/exe/fetch.php?media=", outname, sep=""), outname)
    system(paste("7za e", outname))
    unlink(outname)
  }
  if(is.na(file.info(af_outname.tif)$size)){
    ## resample to the bounding box:
    if(outname.tif %in% c("IFLGRE0a.tif", "LMTGSH0a.tif", "PREGSM0a.tif", "L3POBI0a.tif", "EVMMOD0a.tif")){
      system(paste(gdalwarp, outname.tif, af_outname.tif, '-t_srs', paste('\"', af.csy, '\"', sep=""), '-r near -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 20000, 20000))
    }
    if(outname.tif %in% c("TDMMOD0a.tif", "TDHMOD0a.tif", "TNMMOD0a.tif")){
      system(paste(gdalwarp, outname.tif, af_outname.tif, '-dstnodata \"-32767\" -t_srs', paste('\"', af.csy, '\"', sep=""), '-r bilinear -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 20000, 20000))
    }
    if(outname.tif %in% c("DEMSRE0a.tif", "SLPSRT0a.tif", "TNSMOD0a.tif", "TDSMOD0a.tif", "EVSMOD0a.tif", paste("G0", 1:9, "ESA0a.tif", sep=""), paste("G", 10:22, "ESA0a.tif", sep=""))){
       system(paste(gdalwarp, outname.tif, af_outname.tif, '-t_srs', paste('\"', af.csy, '\"', sep=""), '-r bilinear -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 20000, 20000))
    }    
  }  
}

# import to R (initial grids are large, but the landmask is only .5Mi pixels):
tif0.lst <- list.files(pattern=glob2rx("af_*0a.tif$"))
for(i in 1:length(tif0.lst)){
    afgrid20@data[,strsplit(strsplit(tif0.lst[i], ".tif")[[1]][1], "af_")[[1]][2]] <- readGDAL(tif0.lst[i], silent = TRUE)$band1
} ## Takes 2 mins!

## additional grids (not YET on WorldGrids.org):
afgrid20$TWISRE0x <- readGDAL("af_TWISRE0x.tif", silent = TRUE)$band1
## geology:
afgrid20$SGEUSG0x <- as.factor(readGDAL("af_SGEUSG0x.tif", silent = TRUE)$band1)
xg = read.table("af_SGEUSG3x.txt", header=TRUE)
levels(afgrid20$SGEUSG0x) <- paste(xg$NAME)
xl = read.table("http://worldgrids.org/lib/exe/fetch.php?media=l3pobi.txt", sep="\t", header = TRUE)
afgrid20$L3POBI0a <- as.factor(afgrid20$L3POBI0a)
levels(afgrid20$L3POBI0a) <- c(NA, paste(xl$NAME))

## Convert to spatial pixels:
afgrid20 <- as(afgrid20, "SpatialPixelsDataFrame")
afgrid20 <- afgrid20[!is.na(afgrid20$SMKMOD0x),]
## filter some images:
afgrid20$PREGSM0a <- ifelse(afgrid20$PREGSM0a<0, NA, afgrid20$PREGSM0a)
afgrid20$EVMMOD0a <- ifelse(afgrid20$EVMMOD0a<0, 0, afgrid20$EVMMOD0a)

## check that everything is fine:
str(afgrid20@data)
#spplot(afgrid20["DEMSRE0a"], col.regions=SAGA_pal[[1]])
#spplot(afgrid20["PREGSM0a"], col.regions=SAGA_pal[[1]])
#spplot(afgrid20["TDMMOD0a"], col.regions=SAGA_pal[[1]])
#spplot(afgrid20["EVMMOD0a"], col.regions=SAGA_pal[[1]])
#spplot(afgrid20["L3POBI0a"])
#summary(afgrid20$L3POBI0a)
#spplot(afgrid20["SGEUSG0x"])

save(afgrid20, file="afgrid20.rda", compress="xz")
#system("xcopy afgrid20.rda E:\\DropBox\\AFSIS-ISRIC\\SP_predictions\\20km\\ /Y")

# end of script;