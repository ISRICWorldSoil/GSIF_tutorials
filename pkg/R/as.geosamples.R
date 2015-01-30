# Purpose        : Converts a SoilProfileCollection to loose records (KML placemarks);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Hannes I. Reuter;
# Status         : pre-alpha
# Note           : see also "as.data.frame" operation;

## coerce SoilProfileCollection to "geosamples":
setMethod("as.geosamples", signature(obj = "SoilProfileCollection"), 
  function(obj, registry = as.character(NA), sample.area = 1, mxd = 2, TimeSpan.begin, TimeSpan.end) 
  {

  # reproject if necessary:
  if(requireNamespace("aqp", quietly = TRUE)){
    if(!check_projection(obj@sp)){
       obj@sp <- reproject(obj@sp)
    }
  
    # check for duplicates:
    dp <- duplicated(obj@site[,obj@idcol])
    if(sum(dp)>0){
      warning(paste("Duplicated IDs detected in the 'site' slot and will be removed:", paste(obj@site[dp, obj@idcol], collapse=", ", sep="")))
      obj@site <- obj@site[!dp,]
      obj@sp <- obj@sp[!dp,]
    } 
  
    # estimate thickness in m and depths:
    sampleThickness <- abs(obj@horizons[,obj@depthcols[2]] - obj@horizons[,obj@depthcols[1]])/100 
    depths <- - (obj@horizons[,obj@depthcols[1]] + (obj@horizons[,obj@depthcols[2]] - obj@horizons[,obj@depthcols[1]])/2)/100
  
    # add the time coordinate if missing:
    if(!missing(TimeSpan.begin)&!missing(TimeSpan.end)){
      if(!length(TimeSpan.begin)==length(TimeSpan.end)){ stop("Arguments 'TimeSpan.begin' and 'TimeSpan.end' of the same length required") }
      if(!any(class(TimeSpan.begin)=="POSIXct")&!any(class(TimeSpan.end)=="POSIXct")){ stop("Arguments 'TimeSpan.begin' and 'TimeSpan.end' of class 'POSIXct' required") }
      XYT <- data.frame(cbind(obj@sp@coords, time=(unclass(TimeSpan.begin)+unclass(TimeSpan.end))/2), dtime=(unclass(TimeSpan.end) - unclass(TimeSpan.begin))/2)
    } else {
      XYT <- data.frame(cbind(obj@sp@coords, time=rep(NA, nrow(obj@sp@coords))), dtime=0)
    }
  
    names(XYT)[1:2] <- c("x", "y")
    XYT$ID <- obj@site[,obj@idcol]
    ## TH: this assumes that the 'sp' slot and the 'site' slot have the same order;
    
    # convert site data to geosamples:
    x <- NULL
    site <- obj@site
    # remove columns that are not of interest:
    site[,obj@idcol] <- NULL
  
    ## if empty skip this step:
    if(ncol(site)>0){
    # add the location error:
    locationError = attr(obj@sp@coords, "locationError")
    if(is.null(locationError)) { locationError = rep(as.character(NA), nrow(obj@sp@coords)) }
    XYT$locationError <- locationError
    
    # for each soil variable
    for(j in 1:length(names(site))){
      ll <- length(site[,names(site)[j]])
      observationid = attr(site[,names(site)[j]], "IGSN")
      if(is.null(observationid)) { observationid = rep(as.character(NA), ll) } 
      measurementError = attr(site[,names(site)[j]], "measurementError")
      if(is.null(measurementError)) { measurementError = rep(as.character(NA), ll) }
      sampleArea = attr(site[,names(site)[j]], "sampleArea")
      if(is.null(sampleArea)) { sampleArea = rep(sample.area, ll) }    
      x[[j]] <- data.frame(observationid = as.character(observationid), sampleid = obj@site[,obj@idcol], longitude = XYT[,1], latitude = XYT[,2], locationError = as.numeric(locationError), TimeSpan.begin = as.POSIXct(XYT[,3]-XYT[,"dtime"]/2, origin="1970-01-01"), TimeSpan.end = as.POSIXct(XYT[,3]+XYT[,"dtime"]/2, origin="1970-01-01"), altitude = as.numeric(rep(0, ll)), altitudeMode = rep("relativeToGround", ll), sampleArea = sampleArea, sampleThickness = rep(mxd*sample.area, ll), observedValue = as.character(site[,names(site)[j]]), methodid = rep(names(site)[j], ll), measurementError = as.numeric(measurementError), stringsAsFactors = FALSE) 
    }
      rx <- do.call(rbind, x)
    } else {
      rx <- NULL
    }
        
    # convert horizon data to geosamples:
    y <- NULL
    hors <- obj@horizons
    # remove columns that are not of interest:
    hors[,obj@idcol] <- NULL
    hors[,obj@depthcols[1]] <- NULL
    hors[,obj@depthcols[2]] <- NULL
    
    ## if empty skip this step:
    if(ncol(hors)>0){
      # add coordinates:
      XYTh <- merge(data.frame(ID=obj@horizons[,obj@idcol]), XYT, by="ID", all.x=TRUE)
    
    # for each soil variable
      for(j in 1:length(names(hors))){
        ll <- length(hors[,names(hors)[j]])
        observationid = attr(hors[,names(hors)[j]], "IGSN")
        if(is.null(observationid)) { observationid = rep(as.character(NA), ll) } 
        measurementError = attr(hors[,names(hors)[j]], "measurementError")
        if(is.null(measurementError)) { measurementError = rep(as.character(NA), ll) }
        sampleArea = attr(hors[,names(hors)[j]], "sampleArea")
        if(is.null(sampleArea)) { sampleArea = rep(sample.area, ll) }
          y[[j]] <- data.frame(observationid = as.character(observationid), sampleid = XYTh$ID, longitude = XYTh$x, latitude = XYTh$y, locationError = as.numeric(XYTh$locationError), TimeSpan.begin = as.POSIXct(XYTh$time-XYTh$dtime/2, origin="1970-01-01"), TimeSpan.end = as.POSIXct(XYTh$time+XYTh$dtime/2, origin="1970-01-01"), altitude = as.numeric(depths), altitudeMode = rep("relativeToGround", ll), sampleArea = sampleArea, sampleThickness = sampleThickness, observedValue = as.character(hors[,names(hors)[j]]), methodid = rep(names(hors)[j], ll), measurementError = as.numeric(measurementError), stringsAsFactors = FALSE) 
      }
      ry <- do.call(rbind, y)
    } else {
      ry <- NULL
    }
    
    # merge the sites and horizons tables:
    tb <- rbind(rx, ry)
    tb$methodid <- as.factor(tb$methodid)
  
    # check if the metadata comply with the geosamples standard:
    mnames = c("methodid", "description", "units", "detectionLimit")
    if(any(!(names(obj@metadata) %in% mnames))){ 
      # warning(paste("Missing column names in the 'metadata' slot:", paste(mnames, collapse=", "))) 
      tnames <- levels(tb$methodid)
      metadata <- data.frame(tnames, rep(NA, length(tnames)), rep(NA, length(tnames)), rep(NA, length(tnames)))
      names(metadata) <- mnames
    } else {
      metadata = obj@metadata
    }
   
    # make geosamples:
    gs <- new("geosamples", registry = registry, methods = metadata, data = tb)
      
    return(gs)
  }
})


setMethod("as.geosamples", signature(obj = "SpatialPointsDataFrame"), 
  function(obj, registry = as.character(NA), sample.area = 1, mxd = 2, TimeSpan.begin, TimeSpan.end) 
  {
  
  ## SpatialPoints should be in the WGS84 projection system:
  if(is.na(proj4string(obj))|!check_projection(obj)){ 
    stop(paste("proj4 string", get("ref_CRS", envir = plotKML.opts), "required")) 
  }
  ## this assumes some standard coordinates column names!
  xyn = attr(obj@bbox, "dimnames")[[1]]
  if(!any(xyn %in% c('longitude', 'latitude', 'altitude', 'time'))){
    warning("'longitude', 'latitude', 'altitude', and 'time' coordinates names expected")
  }

  ## space-time coordinates:
  longitude = obj@coords[,1]
  latitude = obj@coords[,2]
  if(!length(xyn)>2){ 
    altitude = rep(0, nrow(obj@coords))
  } else {
    altitude = obj@coords[,3]
  }

  ## altitudeMode
  altitudeMode = attr(obj@coords, "altitudeMode")
    if(is.null(altitudeMode)) { altitudeMode = "relativeToGround" }  
  
  ## look for location error:
  locationError = attr(obj@coords, "locationError")
    if(is.null(locationError)) { locationError = as.character(NA) }
  ## sampleArea/thickness:
  sampleArea = attr(obj@coords, "sampleArea")
    if(is.null(sampleArea)) { sampleArea = sample.area }
  sampleThickness = attr(obj@coords, "sampleThickness")
    if(is.null(sampleThickness)) { sampleThickness = 0 }
  
  ## look for observation ID:
  if("observationid" %in% names(obj@data)){
    attr(obj@data, "observationid") <- obj@data[,"observationid"]
    obj@data[,"observationid"] <- NULL
  }

  ## look for sample ID:
  if("sampleid" %in% names(obj@data)){
    sampleid = obj@data[,"sampleid"]
    obj@data[,"sampleid"] <- NULL
  } else {
    sampleid = attr(obj@data, "sampleid")
    if(is.null(sampleid)){
      sampleid = as.character(1:length(obj))
  }}
  ## look for measurement error:
  if("measurementError" %in% names(obj@data) & length(names(obj@data))==2){
    ## TH: Only works for a single column data frames!
    attr(obj@data, "measurementError") <- obj@data[,"measurementError"]
    obj@data[,"measurementError"] <- NULL
  }

  ## Try to get the times:
  if(!missing(TimeSpan.begin)&!missing(TimeSpan.end)){
    if(!length(TimeSpan.begin)==length(TimeSpan.end)){ stop("Arguments 'TimeSpan.begin' and 'TimeSpan.end' of the same length required") }
    if(!any(class(TimeSpan.begin)=="POSIXct")&!any(class(TimeSpan.end)=="POSIXct")){ stop("Arguments 'TimeSpan.begin' and 'TimeSpan.end' of class 'POSIXct' required") }
  } else {
      ## try to get the times from coordinates:
      if("time" %in% xyn){
        TimeSpan.begin = obj@coords[,"time"]
        TimeSpan.end = obj@coords[,"time"]
      } else {
        TimeSpan.begin = NA
        TimeSpan.end = NA
      }
  } 

  ## first table:
  rx <- data.frame(sampleid, longitude, latitude, locationError = as.numeric(locationError), TimeSpan.begin = as.POSIXct(TimeSpan.begin, origin="1970-01-01"), TimeSpan.end = as.POSIXct(TimeSpan.end, origin="1970-01-01"), altitude, altitudeMode, sampleArea, sampleThickness = as.numeric(sampleThickness), stringsAsFactors = FALSE)

  ## observed values:
  observedValue = as.vector(sapply(as.list(obj@data), function(x){paste(x)}))
  methodid = as.vector(sapply(names(obj), function(x){rep(x, nrow(obj))}))
  observationid = rep(attr(obj@data, "observationid"), length(names(obj)))
    if(is.null(observationid)) { observationid = rep(as.character(NA), length(observedValue)) }
  if(ncol(obj@data)==1 & !is.null(attr(obj@data, "measurementError"))){
    measurementError = attr(obj@data, "measurementError")
  } else { 
    measurementError = rep(as.character(NA), length(observedValue))
  }  
  sampleid2 <- rep(sampleid, ncol(obj@data))
  
  ## second table:
  ry <- data.frame(sampleid = sampleid2, observationid, observedValue, methodid, measurementError, stringsAsFactors = FALSE)

  ## merge the two tables:
  tb <- merge(rx, ry, by="sampleid")
  
  ## look for description of the methods:
  if(ncol(obj@data)==1 & !is.null(attr(obj@data, "description"))){
    description = attr(obj@data, "description")
  } else {
    description = rep(as.character(NA), length(names(obj))) 
  }
  if(ncol(obj@data)==1 & !is.null(attr(obj@data, "units"))){
    units = attr(obj@data, "units")
  } else {   
    units = rep(as.character(NA), length(names(obj)))
  }
  if(ncol(obj@data)==1 & !is.null(attr(obj@data, "detectionLimit"))){
    detectionLimit = attr(obj@data, "detectionLimit")
  } else {
    detectionLimit = rep(as.character(NA), length(names(obj)))
  }
  
  ## try to generate the metadata:
  metadata = data.frame(methodid = names(obj), description, units, detectionLimit, stringsAsFactors = FALSE)

  # make geosamples:
  gs <- new("geosamples", registry = registry, methods = metadata, data = tb[,c("observationid", "sampleid", "longitude", "latitude", "locationError", "TimeSpan.begin", "TimeSpan.end", "altitude", "altitudeMode", "sampleArea", "sampleThickness", "observedValue", "methodid", "measurementError")])  
  
})


## subsetting geosamples:
.subset.geosamples <- function(x, method){
  ret <- x@data[x@data$methodid==method,]
  if(nrow(ret)==0){ warning("Empty object. Methodid possibly not available") }
  attr(ret$methodid, "description") <- x@methods[x@methods$methodid==method,"description"]
  attr(ret$methodid, "units") <- x@methods[x@methods$methodid==method,"units"]
  attr(ret$methodid, "detectionLimit") <- x@methods[x@methods$methodid==method,"detectionLimit"]
  return(ret)  
}

setMethod("subset", signature(x = "geosamples"), .subset.geosamples)

.stack.geosamples <- function(x, lst=levels(x@methods$methodid), geo.columns=c("sampleid","longitude","latitude","altitude")){
  if(requireNamespace("reshape", quietly = TRUE)){
    x.df.lst <- list(NULL)
    for(j in 1:length(lst)){
      x.df.lst[[j]] <- subset(x, method=lst[j])
      x.df.lst[[j]][,lst[j]] <- as.numeric(x.df.lst[[j]][,"observedValue"])
      x.df.lst[[j]] <- x.df.lst[[j]][!is.na(x.df.lst[[j]][,lst[j]]),c(geo.columns,lst[j])]
    }
    x.df <- reshape::merge_recurse(x.df.lst)
    return(x.df)
  }
}

setMethod("stack", signature(x = "geosamples"), .stack.geosamples)

## summary values:
setMethod("show", signature(object = "geosamples"), 
  function(object) 
  {
  cat("  Registry            :", object@registry, "\n")
  cat("  Variables           :", paste(levels(object@data$methodid), collapse=", "), "\n")
  cat("  Total samples       :", nrow(object@data), "\n")
  # create sp object:
  sp <- object@data
  coordinates(sp) <- ~longitude+latitude
  proj4string(sp) <- get("ref_CRS", envir = GSIF.opts)
  cat("  Unique locations    :", length(unique(sp@coords))/2, "\n")
  le = mean(object@data$locationError, na.rm=TRUE)
  cat("  Mean location error :", signif(ifelse(is.nan(le), NA, le), 5), "\n")
  cat("  Min longitude       :", range(object@data$longitude)[1], "\n")  
  cat("  Max longitude       :", range(object@data$longitude)[2], "\n")  
  cat("  Min latitude        :", range(object@data$latitude)[1], "\n")  
  cat("  Max latitude        :", range(object@data$latitude)[2], "\n")  
  # estimate the total area covered by the samples:
  sp.gc <- spTransform(sp, CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))
  Tarea <- signif(diff(sp.gc@bbox[1,])*diff(sp.gc@bbox[2,])/1e6, 4)
  cat("  Total area (app.)   :", paste(Tarea), "(square-km)", "\n")
})


## Extract regression matrix:

.over.geosamplesSP <- function(x, y, methodid, var.type = "numeric"){
  
  if(length(attr(x@coords, "dimnames")[[2]])>2){
    warning("'SpatialPixelsDataFrame' object with two dimensions expected")
  }
  if(!any(y@data$altitudeMode == "relativeToGround")){
    warning("AltitudeMode accepts only 'relativeToGround' values")
  }
  
  pnts = .subset.geosamples(y, method=methodid)
  ## check if it results in an empty set:
  if(nrow(pnts)==0|is.null(pnts)){
    stop(paste("Subsetting the geosamples based on method", methodid, "results in an empy set."))
  }
  
  ## reformat observed values:
  if(var.type=="numeric"){
    pnts$observedValue = as.numeric(pnts$observedValue)
  } else { 
      if(var.type=="factor"){
      pnts$observedValue = as.factor(pnts$observedValue)
      }
  }
  
  coordinates(pnts) <- ~longitude+latitude
  proj4string(pnts) <- get("ref_CRS", envir = GSIF.opts) 
  pnts.t <- spTransform(pnts, x@proj4string)
  attr(pnts.t@coords, "dimnames")[[2]] = attr(x@coords, "dimnames")[[2]]
  attr(pnts.t@bbox, "dimnames")[[1]] = attr(x@bbox, "dimnames")[[2]]
  
  ov <- sp::over(pnts.t, x)                    
  out <- cbind(data.frame(pnts), ov, data.frame(pnts.t@coords))
  if(nrow(out)==0){ 
    warning("Overlay resulted in an empty table") 
  }
  
  return(out)
  
}

setMethod("over", signature(x = "SpatialPixelsDataFrame", y = "geosamples"), .over.geosamplesSP)


.over.geosamplesRaster <- function(x, y, methodid, var.type = "numeric"){
  
  if(!any(y@data$altitudeMode == "relativeToGround")){
    warning("AltitudeMode accepts only 'relativeToGround' values")
  }
  
  pnts = .subset.geosamples(y, method=methodid)
  ## check if it results in an empty set:
  if(nrow(pnts)==0|is.null(pnts)){
    stop(paste("Subsetting the geosamples based on method", methodid, "results in an empy set."))
  }
  
  ## reformat observed values:
  if(var.type=="numeric"){
    pnts$observedValue = as.numeric(pnts$observedValue)
  } else { 
      if(var.type=="factor"){
      pnts$observedValue = as.factor(pnts$observedValue)
      }
  }
  
  coordinates(pnts) <- ~longitude+latitude
  proj4string(pnts) <- get("ref_CRS", envir = GSIF.opts) 
  pnts.t <- spTransform(pnts, CRS(proj4string(x)))
 
  ov <- as.data.frame(raster::extract(x, pnts.t))                    
  out <- cbind(data.frame(pnts), ov, data.frame(pnts.t@coords))
  if(nrow(out)==0){ 
    warning("Overlay resulted in an empty table") 
  }
  
  return(out)
  
}

setMethod("over", signature(x = "RasterStack", y = "geosamples"), .over.geosamplesRaster)


# end of script;