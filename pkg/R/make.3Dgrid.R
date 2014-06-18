# Purpose        : Prepare prediction locations following the GlobalSoilMap specs;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : GlobalSoilMap specifications are constantly updated at [http://globalsoilmap.net/specifications];


## make prediction locations in WGS84 (from point to grid):
setMethod("make.3Dgrid", signature(obj = "SpatialPixelsDataFrame"), function(obj, proj4s = get("ref_CRS", envir = GSIF.opts), pixsize = get("cellsize", envir = GSIF.opts)[2], resampling_method = "bilinear", NAflag = get("NAflag", envir = GSIF.opts), stdepths = get("stdepths", envir = GSIF.opts), tmp.file = TRUE, show.output.on.console = TRUE, ...){   
  
  res <- warp(obj, proj4s = proj4s, pixsize = pixsize, resampling_method = resampling_method, NAflag = NAflag, tmp.file = tmp.file, show.output.on.console = show.output.on.console)    
  # make a list of grids with standard depths:
  res <- as(res, "SpatialPixelsDataFrame") 
  out <- sp3D(obj=res, proj4s = proj4s, stdepths = stdepths)
  return(out)
  
})


## make prediction locations in WGS84 (from point to grid):
setMethod("make.3Dgrid", signature(obj = "RasterBrick"),  function(obj, proj4s = get("ref_CRS", envir = GSIF.opts), pixsize = get("cellsize", envir = GSIF.opts)[2], resampling_method = "bilinear", NAflag = get("NAflag", envir = GSIF.opts), stdepths = get("stdepths", envir = GSIF.opts), tmp.file = TRUE, show.output.on.console = TRUE, ...){
    
    # for each layer layers:
    if (ncol(obj) > 1) {

      r <- stack(obj)
      r <- stack(lapply(r@layers, plotKML::reproject, CRS = CRS, method = resampling_method, ...))
      res <- as(r, "SpatialGridDataFrame")
      names(res) <- names(obj)
    }

    # single layer:
    else {
      r <- raster(obj)
      res <- as(plotKML::reproject(r, CRS = CRS, method = resampling_method, ...), "SpatialGridDataFrame")
      names(res) <- names(obj)
    }
    
  # make a list of grids with standard depths 
  res <- as(res, "SpatialPixelsDataFrame")
  out <- sp3D(obj=res, proj4s = proj4s, stdepths = stdepths)
  
  return(out)

})


## convert to 3D spatial pixels;
setMethod("sp3D", signature(obj = "SpatialPixelsDataFrame"), function(obj, proj4s = proj4string(obj), stdepths = get("stdepths", envir = GSIF.opts), stsize = get("stsize", envir = GSIF.opts), make.grid = TRUE){
    
  ## convert to a data frame:
  x <- data.frame(obj)
  ## rename the column names so they correspond to the geosamples class:
  sel <- names(x) %in% attr(obj@coords, "dimnames")[[2]]
  ## check if these are 2D or 3D grids:
  if(sum(sel)==3){ 
      names(x)[sel] <- c("longitude", "latitude", "altitude")
  } else {
    if(sum(sel)==2) { 
      names(x)[sel] <- c("longitude", "latitude") 
    }
  }

  out <- list(NULL)
  for(j in 1:length(stdepths)){
    XYD <- x
    XYD$altitude <- rep(stdepths[j], nrow(XYD))
    coordinates(XYD) <- ~ longitude + latitude + altitude
    if(make.grid == TRUE){
      ## sp complains by default, so better mask out the warnings:
      try( suppressWarnings( gridded(XYD) <- TRUE ) )
      if(!class(.Last.value)[1]=="try-error"){
        ## fix the cell size and cellcentre.offset:
        XYD@grid@cellsize[3] <- stsize[j]
      }
    }
    XYD@bbox[3,1] <- stdepths[j]-stsize[j]/2
    XYD@bbox[3,2] <- stdepths[j]+stsize[j]/2    
    proj4string(XYD) <- proj4s
    out[[j]] <- XYD
  }
  
  return(out)
})


## downsample and convert to 3D spatial pixels;
setMethod("sp3D", signature(obj = "list"), function(obj, proj4s = proj4string(obj[[1]]), stdepths = get("stdepths", envir = GSIF.opts), stsize = get("stsize", envir = GSIF.opts), tmp.file = TRUE, ...){

  # check if the objects are valid:
  if(length(unique(sapply(obj, FUN=function(x){proj4string(x)})))>1){
    stop("Not all coordinate systems in the 'obj' list match")
  }
  if(!any(sapply(obj, class)=="SpatialPixelsDataFrame")){
    stop("List of covariates of class 'SpatialPixelsDataFrame' expected")
  }
  if(any(!(sapply(obj, FUN=function(x){x@grid@cellsize[1]})==sapply(obj, FUN=function(x){x@grid@cellsize[2]})))){
    stop("Covariate layers with equal grid cell size expected")
  }
  
  # pick the most detailed scale:
  cellsize.l <- sapply(obj, FUN=function(x){x@grid@cellsize[1]})
  tc <- which(cellsize.l == min(cellsize.l))
  x <- obj[[tc[1]]]
  fullgrid(x) <- TRUE
  obj[[tc[1]]] <- NULL
  
  # resample grids not available in the finest resolution:
  ret <- list(NULL)
  for(j in 1:length(obj)){
    # check if it is available in the finest resolution or if it does not overlap:
    if(cellsize.l[j] > min(cellsize.l)|!identical(x@bbox, obj[[j]]@bbox)){
      ret[[j]] <- warp(obj[[j]], proj4s = proj4string(x), pixsize = min(cellsize.l), GridTopology = x@grid, resampling_method = "cubicspline", tmp.file = tmp.file)
    }
    else {
    ret[[j]] <- obj[[j]]
    }
  }
  
  ret <- lapply(ret, FUN=function(x){slot(x, "data")})
  x@data <- cbind(x@data, do.call(cbind, ret))
  x <- as(x, "SpatialPixelsDataFrame")
 
  out <- sp3D(obj = x, stdepths = stdepths, stsize = stsize)      
 
  return(out)
})



## make GlobalSoilMap object:
GlobalSoilMap <- function(obj, varname, TimeSpan = list(begin=Sys.Date()-1, end=Sys.Date())){
  if(!class(obj)=="list"){
    stop("Object of class 'list' required")
  }
  
  out = new("GlobalSoilMap", varname = varname, TimeSpan = TimeSpan, sd1=as(obj[[1]], "SpatialPixelsDataFrame"), sd2=as(obj[[2]], "SpatialPixelsDataFrame"), sd3=as(obj[[3]], "SpatialPixelsDataFrame"), sd4=as(obj[[4]], "SpatialPixelsDataFrame"), sd5=as(obj[[5]], "SpatialPixelsDataFrame"), sd6=as(obj[[6]], "SpatialPixelsDataFrame"))
  
  return(out)
}

## make SoilGrids object:
SoilGrids <- function(obj, varname, TimeSpan = get("TimeSpan", envir = GSIF.opts)){
  if(!class(obj)=="list"){
    stop("Object of class 'list' required")
  }
  
  out = new("SoilGrids", varname = varname, TimeSpan = TimeSpan, sd1=as(obj[[1]], "SpatialPixelsDataFrame"), sd2=as(obj[[2]], "SpatialPixelsDataFrame"), sd3=as(obj[[3]], "SpatialPixelsDataFrame"), sd4=as(obj[[4]], "SpatialPixelsDataFrame"), sd5=as(obj[[5]], "SpatialPixelsDataFrame"), sd6=as(obj[[6]], "SpatialPixelsDataFrame"))
  
  return(out)
}


# end of script;