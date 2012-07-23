# Purpose        : Prepare prediction locations following the GlobalSoilMap specs;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : GlobalSoilMap specifications are constantly updated at [http://globalsoilmap.net/specifications];


## make prediction locations in WGS84 (from point to grid):
setMethod("make.3Dgrid", signature(obj = "SpatialPixelsDataFrame"),  function(obj, proj4s = get("ref_CRS", envir = GSIF.opts), pixsize = get("cellsize", envir = GSIF.opts)[2], resample = "bilinear", NAflag = get("NAflag", envir = GSIF.opts), stdepths = get("stdepths", envir = GSIF.opts), tmp.file = TRUE, show.output.on.console = FALSE, ...){   
  
  # look for FWTools path:  
  require(plotKML)
  gdalwarp <- get("gdalwarp", envir = plotKML.opts)
  if(nchar(gdalwarp)==0){
      plotKML.env(silent = FALSE)
      gdalwarp <- get("gdalwarp", envir = plotKML.opts)
  }
  
  if(!nchar(gdalwarp)==0){
    message(paste("Resampling", length(names(obj)), "layers to", proj4s, "with grid cell size of", pixsize, "..."))
    pb <- txtProgressBar(min=0, max=ncol(obj), style=3)
    for(i in 1:ncol(obj)){
  
    if(tmp.file==TRUE){
        tf <- tempfile() 
        }
    else { 
        tf <- paste(normalizeFilename(deparse(substitute(obj, env = parent.frame()))), names(obj)[i], sep="_")
       }

        # write SPDF to a file:
        if(is.factor(obj@data[,i])){
          x <- obj[i]
          x@data[,1] <- as.integer(x@data[,1])
          writeGDAL(x, paste(tf, ".tif", sep=""), "GTiff")
        }        
        else {
          writeGDAL(obj[i], paste(tf, ".tif", sep=""), "GTiff")
        }
        
        # resample to WGS84 system:
        if(is.factor(obj@data[,i])){
          system(paste(gdalwarp, ' ', tf, '.tif', ' -t_srs \"', proj4s, '\" ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r near', ' -tr ', pixsize, ' ', pixsize, sep=""), show.output.on.console = show.output.on.console)
        }
        else {
          system(paste(gdalwarp, ' ', tf, '.tif', ' -t_srs \"', proj4s, '\" ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r ', resample, ' -tr ', pixsize, ' ', pixsize, sep=""), show.output.on.console = show.output.on.console)
        }
        # read images back to R:
        if(i==1){
          res <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = TRUE)
          names(res) <- names(obj)[i]
        }
        else{
          res@data[names(obj)[i]] <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = TRUE)$band1
        }
        
        # reformat to the original factors:
        if(is.factor(obj@data[,i])){
          res@data[,i] <- factor(res@data[,i], levels=levels(obj@data[,i]))
        }
    setTxtProgressBar(pb, i)          
    }
  close(pb)
  } else { stop("Could not locate FWTools. First install and test FWTools. See 'plotKML.env()' for more info.") }
  
  # make a list of grids with standard depths 
  out <- sp3D(obj=res, proj4s = proj4s, stdepths = stdepths)
  
  return(out)

})


## make prediction locations in WGS84 (from point to grid):
setMethod("make.3Dgrid", signature(obj = "RasterBrick"),  function(obj, proj4s = get("ref_CRS", envir = GSIF.opts), pixsize = get("cellsize", envir = GSIF.opts)[1], resample = "bilinear", NAflag = get("NAflag", envir = GSIF.opts), stdepths = get("stdepths", envir = GSIF.opts), tmp.file = TRUE, show.output.on.console = FALSE, ...){
    
    library(plotKML)
    # for each layer layers:
    if (ncol(obj) > 1) {

      r <- stack(obj)
      r <- stack(lapply(r@layers, plotKML::reproject, CRS = CRS, method = resample, ...))
      res <- as(r, "SpatialGridDataFrame")
      names(res) <- names(obj)
    }

    # single layer:
    else {
      r <- raster(obj)
      res <- as(plotKML::reproject(r, CRS = CRS, method = resample, ...), "SpatialGridDataFrame")
      names(res) <- names(obj)
    }
    
  # make a list of grids with standard depths 
  out <- sp3D(obj=res, proj4s = proj4s, stdepths = stdepths)
  
  return(out)

})


## convert to 3D spatial pixels;
sp3D <- function(obj, proj4s = proj4string(obj), stdepths = get("stdepths", envir = GSIF.opts), stsize = get("stsize", envir = GSIF.opts)){
  if(!(class(obj)=="SpatialPixelsDataFrame"|class(obj)=="SpatialGridDataFrame"))
    stop("Object of class 'SpatialPixelsDataFrame' or 'RasterBrick' required")
  
  # convert to a data frame:
  x <- as(obj, "SpatialPixelsDataFrame")
  x <- as.data.frame(x)
  # rename the column names so they correspond to the geosamples class:
  names(x)[c(length(names(x))-1, length(names(x)))] <- c("longitude", "latitude")

  out <- list(NULL)
  for(j in 1:length(stdepths)){
    XYD <- x
    XYD$altitude <- stdepths[j] 
    # sp complains by default, so better mask out the warnings:
    suppressWarnings(gridded(XYD) <- ~ longitude + latitude + altitude)
    proj4string(XYD) <- proj4s
    XYD@grid@cellsize[3] <- stsize[j]
    out[[j]] <- XYD
  }
  
  return(out)
}

# make GlobalSoilMap class:
GlobalSoilMap <- function(obj, varname, period = c(Sys.Date()-1, Sys.Date())){
  if(!class(obj)=="list"){
    stop("Object of class 'list' required")
  }
  out = new("GlobalSoilMap", varname = varname, TimeSpan.begin = as.POSIXct(period[1]), TimeSpan.end = as.POSIXct(period[2]), sd1=as(obj[[1]], "SpatialPixelsDataFrame"), sd2=as(obj[[2]], "SpatialPixelsDataFrame"), sd3=as(obj[[3]], "SpatialPixelsDataFrame"), sd4=as(obj[[4]], "SpatialPixelsDataFrame"), sd5=as(obj[[5]], "SpatialPixelsDataFrame"), sd6=as(obj[[6]], "SpatialPixelsDataFrame"))
  return(out)
}


# end of script;