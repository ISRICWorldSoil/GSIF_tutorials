# title         : soilmask_Africa.R
# purpose       : preparation of a soil mask for the Sub-Saharan Africa;
# reference     : [http://africasoils.net]
# producer      : Prepared by T. Hengl (tom.hengl@wur.nl)
# address       : In Wageningen, NL.
# inputs        : list of countries and the mean Leaf Area index [http://worldgrids.org/doku.php?id=wiki:lammod3]; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : Running on a HP Z420 station with 32GB RAM (64-bit Windows 7 pro);
# remarks 2     : According to the MODIS LAI specifications [https://lpdaac.usgs.gov/products/modis_products_table/mod15a2], perennial salt or inland fresh water, barren, sparse vegetation (rock, tundra, desert), perennial snow, ice, "permanent" wetlands/inundated marshlands, urban/built-up land areas have already been removed from the LAI images;

library(raster)
library(rgdal)
library(maptools)
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

# ------------------------------------------------------------
# Vector maps (GADM):
# ------------------------------------------------------------

## Donwload  shape files:
download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", "ne_10m_admin_0_countries.zip") 
unzip("ne_10m_admin_0_countries.zip")

## Sub-Saharan Africa list of countries of interest:
cnt.lst <- c("CMR", "ETH", "GHA", "KEN", "MLI", "MWI", "NGA", "RWA", "SEN", "TZA", "UGA", "ERI", "DJI", "SSD", "SDN", "CAF", "BEN", "TGO", "GIN", "GNB", "LBR", "SLE", "GNQ", "GAB", "COG","COD", "BDI", "ZWE", "MOZ", "NAM", "SWZ", "BWA", "LSO", "CIV", "GMB", "NER", "TCD", "CPV", "BFA", "ZMB", "MDG", "COM", "SOM", "DJI", "STP", "AGO", "ZAF", "SOL")
cnt2.lst <- c("EGY", "SDN", "LBY", "TUN", "DZA", "MRT", "SAH", "MAR")
admin.wrld = readOGR("ne_10m_admin_0_countries.shp", "ne_10m_admin_0_countries")
## one level is missing! should be "SSD" but it is "SDS" in the "ne_10m_admin_0_countries.shp"!
x = which(!(cnt.lst %in% levels(admin.wrld$sov_a3)))
cnt.lst[x]
## subset to SS Africa:
admin = admin.wrld[which(admin.wrld$SOV_A3 %in% c(cnt.lst, "SDS")),]
admin.af = admin.wrld[which(admin.wrld$SOV_A3 %in% c(cnt.lst, cnt2.lst, "SDS")),]
spplot(admin["SOV_A3"])
spplot(admin.af["SOV_A3"])
save(admin, file="admin.rda", compress="xz")
save(admin.af, file="admin.af.rda")

bbox = admin@bbox
## TH: cout out this one South African island:
bbox[2,1] = -35

## get ID for each country:
load("ISO.country.rda")
str(ISO.country)
cnt_ID.lst <- merge(x=data.frame(ISO=cnt.lst), y=ISO.country, all.y=FALSE)


# ------------------------------------------------------------
# Soil productive areas for Sub-Saharan Africa:
# ------------------------------------------------------------

## download LAMMOD (Leaf Area) and CNTGAD (countries) images:
g.lst <- c("EVMMOD3", "LAMMOD3", "CNTGAD3")
for(j in 1:length(g.lst)){
  outname = paste(g.lst[j], "a.tif.gz", sep="")
  outname.tif = paste(g.lst[j], "a.tif", sep="")
  afname = paste("af_", g.lst[j], "a.tif", sep="")
  if(is.na(file.info(afname)$size)){
    download.file(paste("http://worldgrids.org/lib/exe/fetch.php?media=", outname, sep=""), outname)
    system(paste("7za e", outname))
    unlink(outname)
    ## resample to the bounding box:
    system(paste(gdalwarp, outname.tif, afname, '-r near -te', bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2], '-tr', 1/120, 1/120))  
    unlink(outname.tif)
  }
}

## read to R and subset:
go.lst <- list.files(pattern=glob2rx('af_*.tif$'))
grid1km <- readGDAL(go.lst[1])
for(j in 2:length(go.lst)){ grid1km@data[,go.lst[j]] = readGDAL(go.lst[j])$band1 }
names(grid1km) = go.lst
## subset to countries of interest with positive biomass or LAI!
grid1km$mask <- ifelse(grid1km$af_CNTGAD3a.tif %in% cnt_ID.lst$id & grid1km$af_LAMMOD3a.tif>0, 1, NA)
## total number of pixels:
summary(grid1km$mask)

writeGDAL(grid1km["mask"], "afmask.tif", type="Byte", mvFlag=0)
system(paste("7za a", "-tgzip afmask.tif.gz afmask.tif"))

## reproject to cartesian coordinates:
system(paste(gdalwarp, 'afmask.tif afmask_laea.tif -r near -t_srs \"', af.csy, '\" -dstnodata \"0\" -tr', 1000, 1000))
GDALinfo("afmask_laea.tif")
system(paste("7za a", "-tgzip afmask_laea.tif.gz afmask_laea.tif"))

## convert to a data frame:
# afmask <- as.data.frame(grid1km["mask"])
# save(afmask, "afmask.rda", compress="xz")

## soil mask and predictors prepared by the AfSIS team:
library(RSAGA)
tdir <- "G:/AFSIS/JCdata"
pdir <- "E:/DropBox/AFSIS-ISRIC/SP_predictions/1km"

## read each grid in a loop and save and compress with the default name:
tnames = read.table("layer_names1km.txt", header = TRUE)
## two times mask!
for(j in 1:length(tnames$WorldGridsName)){
  outname = set.file.extension(tnames$WorldGridsName[j], "tif")
  tgzname = paste(pdir, set.file.extension(tnames$WorldGridsName[j], "tif.gz"), sep="/")
  if(is.na(file.info(tgzname)$size)&!(j==5|j==11)){
    x = readGDAL(paste(tdir, "predictgrid.tif", sep="/"), band=j)
    names(x) = paste(tnames$WorldGridsName[j])
    proj4string(x) = af.csy
    if(outname %in% c("BI1WCL1x.tif", "WTISRE1x.tif", "SLPSRE1x.tif")){
      writeGDAL(x, outname, type = "Byte", mvFlag="0")
    }
    if(outname %in% c("BI2WCL1x.tif", "DEMSRE1x.tif", "EVMMOD1x.tif", "R01MOD1x.tif", "R02MOD1x.tif", "R03MOD1x.tif", "R07MOD1x.tif", "NPMMOD1x.tif", "RLFSRE1x.tif", "TDMMOD1x.tif", "TNMMOD1x.tif")){
      writeGDAL(x, outname, type = "Int16", mvFlag="-9999")  
    }
    # compress maps:
    system(paste("7za a", "-tgzip", tgzname, outname))
    unlink(outname)
  }
}

## create a mask at 20 km:
x = readGDAL(paste(tdir, "predictgrid.tif", sep="/"), band=11)
names(x) = "SMKMOD3x"
proj4string(x) = af.csy
x$SMKMOD3x <- ifelse(is.na(x$SMKMOD3x), 0, x$SMKMOD3x)
writeGDAL(x, "SMKMOD3x.tif", type = "Byte")
unlink("SMKMOD0x.tif")
system(paste(gdalwarp, 'SMKMOD3x.tif SMKMOD0x.tif -r near -dstnodata \"0\" -tr', 20000, 20000))
unlink("SMKMOD1x.tif")
system(paste(gdal_translate, 'SMKMOD3x.tif SMKMOD3x_b.tif -ot \"Float32\"'))
system(paste(gdalwarp, 'SMKMOD3x_b.tif SMKMOD1x.tif -r cubicspline -tr', 5000, 5000))
afgrid5km <- readGDAL("SMKMOD1x.tif")
names(afgrid5km) = "SMKMOD1x"
afgrid5km@data[,1] <- ifelse(afgrid5km@data[,1]>0, 1, NA)
afgrid5km <- as(afgrid5km, "SpatialPixelsDataFrame")
save(afgrid5km, file="afgrid5km.rda", compress="xz")
system("7za a -tgzip SMKMOD1x.tif.gz SMKMOD1x.tif")
system("7za a -tgzip SMKMOD3x.tif.gz SMKMOD3x.tif")
unlink("SMKMOD3x_b.tif")
rm(x)

# mask map:
afmask20km = readGDAL("SMKMOD0x.tif")
afmask20km.txt = as.data.frame(afmask20km)
str(afmask20km.txt)
## 46,000 pixels!
library(plotKML)
pl = grid2poly(afmask20km)
kml(pl, file="afmask20km.kml", kmz=TRUE)

# end of script; 
