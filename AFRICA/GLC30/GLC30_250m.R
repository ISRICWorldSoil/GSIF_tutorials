# title         : GLC30_250m.R
# purpose       : Processing GLC30 for Africa;
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl (tom.hengl@wur.nl)
# address       : In Wageningen, NL.
# inputs        : GLC30 zipped files download from http://www.globallandcover.com/GLC30Download/download.aspx
# outputs       : 250 m resolution imagery;
# remarks 1     : ;

library(RSAGA)
library(gdalUtils)
library(rgdal)
library(utils)
library(snowfall)
library(raster)
csize = 250
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
classes = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
zip.lst <- list.files(pattern=glob2rx("*.zip$"), full.names=FALSE)
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)

resample.glc30 <- function(i, prj=af.csy, csize=csize, classes=classes){
  lst <- unzip(zip.lst[i], list = TRUE)
  inname <- lst$Name[grep(lst$Name, pattern=glob2rx("*.tif$"))]
  if(nchar(inname)>0){
    unzip(zip.lst[i], inname)
    ## resample to local coordinates:
    gdalwarp(srcfile=inname, dstfile=gsub(".tif", ".sdat", inname), t_srs=prj, of="SAGA", r="near", tr=c(30,30), srcnodata=0, ot="Byte", dstnodata=255, overwrite=TRUE)
    gdal_translate(gsub(".tif", ".sdat", inname), gsub(".tif", paste0("_30m.tif"), inname), co="COMPRESS=DEFLATE", ot="Byte")
    ## mask values per each class:
    for(k in 1:length(classes)){
      outn <- gsub(".tif", paste0("_", classes[k] ,"_250m.tif"), inname)
      if(!file.exists(outn)){
        rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=gsub(".tif", ".sgrd", inname), RESULT=gsub(".tif", paste0("_", classes[k] ,".sgrd"), inname), FORMULA=paste("ifelse(a=",classes[k],",100,0)", sep="")), show.output.on.console = FALSE)
        ## resample:
        gdalwarp(srcfile=gsub(".tif", paste0("_", classes[k] ,".sdat"), inname), dstfile=outn, r="average", ot="Byte", tr=c(250,250), dstnodata=255, overwrite=TRUE)
        ## clean up:
        unlink(gsub(".tif", paste0("_", classes[k], ".sgrd"), inname))
        unlink(gsub(".tif", paste0("_", classes[k], ".mgrd"), inname))
        unlink(gsub(".tif", paste0("_", classes[k], ".sdat"), inname))
        unlink(gsub(".tif", paste0("_", classes[k], ".prj"), inname))
      }
    }
    unlink(gsub(".tif", ".sdat", inname))
    unlink(gsub(".tif", ".sgrd", inname))
    unlink(gsub(".tif", ".prj", inname))
    unlink(gsub(".tif", ".sdat.aux.xml", inname))
    unlink(inname)
  }
}

lapply(1:length(zip.lst), resample.glc30, classes=classes)

#sfInit(parallel=TRUE, cpus=6)
#sfLibrary(RSAGA)
#sfLibrary(utils)
#sfLibrary(gdalUtils)
#sfExport("zip.lst", "classes", "csize", "af.csy")
#t <- sfLapply(1:length(zip.lst), resample.glc30, classes=classes)
#sfStop()

## Create mosaick:
outdir <- "H:\\AFSIS\\africagrids"  
## list all geotifs, create mosaicks and compress to output dir:
tifout.lst <- paste0("*_", (1:10)*10, "_250m.tif$")
for(j in 1:length(tifout.lst)){
  vrt.name <- paste0("af_GLC_", j*10, "_250m")
  outn <- paste0(vrt.name, ".tif")
  if(!file.exists(paste0(outdir, "\\", outn, ".gz"))){
    tmp.lst <- list.files(pattern=glob2rx(tifout.lst[j]), recursive = TRUE)
    unlink("my_liste.txt")
    cat(tmp.lst, sep="\n", file="my_liste.txt")
    gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0(vrt.name, ".vrt"))
    gdalwarp(paste0(vrt.name, ".vrt"), dstfile=outn, tr=c(250,250), r="near", dstnodata=255, te=c(-3977625,-4321625,3397625,3554625), ot="Byte")
    system(paste("7za a", "-tgzip", paste0(outn, ".gz"), outn)) 
    system(paste("xcopy", paste0(outn, ".gz"), shortPathName(normalizePath(outdir))))
    unlink(paste0(vrt.name, ".vrt"))
    unlink(paste0(outn, ".gz"))
    #unlink(outn)
  }
}


## end of script;