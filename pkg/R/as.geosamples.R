# Purpose        : Converts a SoilProfileCollection to loose records (KML placemarks);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Hannes I. Reuter;
# Status         : pre-alpha
# Note           : see also "as.data.frame" operation;

## coerce SoilProfileCollection to "geosamples":
setMethod("as.geosamples", signature(obj = "SoilProfileCollection"), 
  function(obj, registry = as.character(NA), sample.area = 1, mxd = 2, dtime = 3600) 
  {
  # reproject if necessary:
  require(plotKML)
  if(!check_projection(obj@sp)){
     obj@sp <- reproject(obj@sp)
  }
 
  # estimate volume in m^3 and depths:
  vols <- abs(obj@horizons[,obj@depthcols[2]] - obj@horizons[,obj@depthcols[1]])/100 * sample.area
  depths <- - (obj@horizons[,obj@depthcols[1]] + (obj@horizons[,obj@depthcols[2]] - obj@horizons[,obj@depthcols[1]])/2)/100

  # add the time coordinate if missing:
  if(ncol(obj@sp@coords)==2){
    XYT <- data.frame(cbind(obj@sp@coords, time=rep(NA, nrow(obj@sp@coords))))
  }
  else {
    XYT <- data.frame(obj@sp@coords)
  }
  XYT$ID <- profile_id(obj)

  # convert site data to geosamples:
  x <- NULL
  site <- obj@site
  # remove columns that are not of interest:
  site[,obj@idcol] <- NULL
  locationError = attr(obj@sp@coords, "locationError")
  if(is.null(locationError)) { locationError = rep(as.character(NA), nrow(obj@sp@coords)) }
  
  # for each soil variable
  for(j in 1:length(names(site))){
    ll <- length(site[,names(site)[j]])
    sampleid = attr(site[,names(site)[j]], "IGSN")
    if(is.null(sampleid)) { sampleid = rep(as.character(NA), ll) } 
    measurementError = attr(site[,names(site)[j]], "measurementError")
    if(is.null(measurementError)) { measurementError = rep(as.character(NA), ll) }
    x[[j]] <- data.frame(sampleid = sampleid, producerid = profile_id(obj), longitude = XYT[,1], latitude = XYT[,2], locationError = locationError, TimeSpan.begin = as.POSIXct(XYT[,3]-dtime/2, origin="1970-01-01"), TimeSpan.end = as.POSIXct(XYT[,3]+dtime/2, origin="1970-01-01"), altitude = rep("0", ll), altitudeMode = rep("relativeToGround", ll), volume = rep(mxd*sample.area, ll), observedValue = site[,names(site)[j]], methodid = rep(names(site)[j], ll), measurementError = measurementError) 
  }
  rx <- do.call(rbind, x)
  # reformat values:
  rx$sampleid <- as.character(rx$sampleid)
  rx$locationError <- as.numeric(rx$locationError)
  rx$altitude <- as.numeric(rx$altitude)
  rx$observedValue <- as.character(rx$observedValue)
  rx$measurementError <- as.numeric(rx$measurementError)
    
  # convert horizon data to geosamples:
  y <- NULL
  hors <- obj@horizons
  # remove columns that are not of interest:
  hors[,obj@idcol] <- NULL
  hors[,obj@depthcols[1]] <- NULL
  hors[,obj@depthcols[2]] <- NULL
  # add coordinates:
  XYTh <- merge(data.frame(ID=obj@horizons[,obj@idcol], dtime=dtime), XYT, by=obj@idcol, all.x=TRUE)
  
  # for each soil variable
  for(j in 1:length(names(hors))){
    ll <- length(hors[,names(hors)[j]])
    sampleid = attr(hors[,names(hors)[j]], "IGSN")
    if(is.null(sampleid)) { sampleid = rep(as.character(NA), ll) } 
    measurementError = attr(hors[,names(hors)[j]], "measurementError")
    if(is.null(measurementError)) { measurementError = rep(as.character(NA), ll) }
    y[[j]] <- data.frame(sampleid = sampleid, producerid = XYTh[,1], longitude = XYTh[,3], latitude = XYTh[,4], locationError = locationError, TimeSpan.begin = as.POSIXct(XYTh[,5]-XYTh[,2]/2, origin="1970-01-01"), TimeSpan.end = as.POSIXct(XYTh[,5]+XYTh[,2]/2, origin="1970-01-01"), altitude = depths, altitudeMode = rep("relativeToGround", ll), volume = vols, observedValue = hors[,names(hors)[j]], methodid = rep(names(hors)[j], ll), measurementError = measurementError) 
  }
  ry <- do.call(rbind, y)
  # reformat values:
  ry$sampleid <- as.character(ry$sampleid)
  ry$locationError <- as.numeric(ry$locationError)
  ry$altitude <- as.numeric(ry$altitude)
  ry$observedValue <- as.character(ry$observedValue)
  ry$measurementError <- as.numeric(ry$measurementError)
  
  # merge the sites and horizons tables:
  tb <- rbind(rx, ry)
 
  # make geosamples:
  gs <- new("geosamples", registry = registry, methods = obj@metadata, data = tb)
    
  return(gs)
})



## subsetting geosamples:
setMethod("subset", signature(x = "geosamples"), 
  function(x, method) 
  {
  
  ret <- x@data[x@data$methodid==method,]
  attr(ret$methodid, "description") <- x@methods[x@methods$methodid==method,"description"]
  attr(ret$methodid, "units") <- x@methods[x@methods$methodid==method,"units"]
  attr(ret$methodid, "detectionLimit") <- x@methods[x@methods$methodid==method,"detectionLimit"]
  return(ret)  
})


## summary values:
setMethod("show", signature(object = "geosamples"), 
  function(object) 
  {
  cat("  Registry            :", object@registry, "\n")
  cat("  Variables           :", paste(levels(object@data$methodid), collapse=", "), "\n")
  cat("  Total samples       :", nrow(object@data), "\n")
  sp <- object@data
  coordinates(sp) <- ~longitude+latitude
  proj4string(sp) <- get("ref_CRS", envir = plotKML.opts)
  cat("  Unique locations    :", length(unique(sp@coords))/2, "\n")
  cat("  Mean location error :", mean(object@data$locationError, na.rm=TRUE), "\n")
  cat("  Min longitude       :", range(object@data$longitude)[1], "\n")  
  cat("  Max longitude       :", range(object@data$longitude)[2], "\n")  
  cat("  Min latitude        :", range(object@data$latitude)[1], "\n")  
  cat("  Max latitude        :", range(object@data$latitude)[2], "\n")  
  # estimate the total area covered by the samples:
  sp.gc <- spTransform(sp, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"))
  Tarea <- signif(diff(sp.gc@bbox[1,])*diff(sp.gc@bbox[2,])/1e6, 4)
  cat("  Total area          :", paste(Tarea), "(square-km)", "\n")
})


# end of script;